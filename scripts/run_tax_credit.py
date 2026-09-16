#!/usr/bin/env python
"""Tourmaline tax-credit orchestrator: simulation, manifest, evaluation, plotting."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import os
import shutil
import subprocess
import sys
from glob import glob
from os.path import abspath, basename, dirname, expandvars, getmtime, join, exists, splitext
from pathlib import Path

import pandas as pd
import yaml

# Ensure tax-credit is importable
def _ensure_tax_credit(package_dir: str) -> None:
    root = abspath(expandvars(package_dir))
    if root not in sys.path:
        sys.path.insert(0, root)
    try:
        import tax_credit  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            f"tax-credit not installed. Run: pip install -e {root}\n{exc}"
        ) from exc


def load_config(path: str) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    validate_blca_config(cfg)
    return cfg


_BLCA_CUTOFF_KEYS = ("blca_perc_identity", "blca_query_cov")


def validate_blca_config(cfg: dict) -> None:
    """Require bt2-blca's own cutoffs whenever bt2-blca is benchmarked.

    Older configs shared ``perc_identity`` / ``query_cov`` between the consensus
    methods and bt2-blca. Refuse them rather than silently running BLCA with
    different cutoffs than the config appears to request.
    """
    if "bt2-blca" not in classify_methods(cfg):
        return
    missing = [key for key in _BLCA_CUTOFF_KEYS if cfg.get(key) is None]
    if missing:
        raise ValueError(
            f"bt2-blca is in classify_methods but {', '.join(missing)} is not set. "
            "bt2-blca no longer reads perc_identity / query_cov (those are for "
            "consensus-blast / consensus-vsearch only). Add blca_perc_identity and "
            "blca_query_cov to the BT2-BLCA OPTIONS section of the config; see "
            "config_04_tax_credit.yaml."
        )


def run_output_dir(cfg: dict) -> str:
    return join(
        abspath(expandvars(cfg["output_dir"])),
        cfg["run_name"] + "-tax-credit",
    )


def data_dir(cfg: dict) -> str:
    return join(run_output_dir(cfg), cfg.get("data_subdir", "data"))


def results_dir(cfg: dict) -> str:
    """Directory for taxonomy assignment outputs (under data/)."""
    return join(data_dir(cfg), cfg.get("results_subdir", "results"))


_DB_SIMULATION_DEFAULTS = {
    "read_length": 250,
    "min_read_length": 140,
    "trim_primers": True,
    "truncate": True,
}


def _db_simulation_params(db: dict) -> dict:
    """Per-reference simulation settings with defaults."""
    params = {}
    for key, default in _DB_SIMULATION_DEFAULTS.items():
        value = db.get(key, default)
        if key in ("trim_primers", "truncate"):
            params[key] = bool(value)
        else:
            params[key] = value
    if params["truncate"] and params["read_length"] is None:
        raise ValueError(
            f"reference_databases entry {db['id']!r}: read_length is required "
            "when truncate is true"
        )
    return params


def _run_qiime(cmd: str) -> None:
    print(cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True)


def _qiime_export(qza_path: str, dest_path: str, output_format: str) -> None:
    """Export via ``qiime tools export`` to a single output file."""
    if exists(dest_path):
        os.remove(dest_path)
    parent = dirname(dest_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    _run_qiime(
        "qiime tools export "
        f"--input-path '{qza_path}' "
        f"--output-path '{dest_path}' "
        f"--output-format {output_format}"
    )


def _write_headerless_taxonomy_tsv(src_path: str, dest_path: str) -> None:
    with open(src_path, encoding="utf-8", errors="replace") as src, open(
        dest_path, "w", encoding="utf-8"
    ) as dest:
        for line in src:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("Feature ID"):
                continue
            parts = stripped.split("\t")
            if len(parts) < 2:
                continue
            dest.write(f"{parts[0]}\t{parts[1]}\n")


def _ensure_text_reference(path: str, staging_dir: str, kind: str) -> str:
    """Return a plain FASTA/TSV path, exporting QIIME 2 .qza artifacts when needed."""
    if not path.lower().endswith(".qza"):
        return path

    os.makedirs(staging_dir, exist_ok=True)
    stem = splitext(basename(path))[0]
    if kind == "sequences":
        dest = join(staging_dir, stem + ".fasta")
    elif kind == "taxonomy":
        dest = join(staging_dir, stem + ".tsv")
    else:
        raise ValueError(f"Unknown reference kind: {kind}")

    if exists(dest) and getmtime(dest) >= getmtime(path):
        return dest

    if kind == "sequences":
        _qiime_export(path, dest, "DNAFASTAFormat")
    else:
        try:
            _qiime_export(path, dest, "HeaderlessTSVTaxonomyFormat")
        except subprocess.CalledProcessError:
            tmp_tsv = dest + ".qiime_export.tsv"
            if exists(tmp_tsv):
                os.remove(tmp_tsv)
            _qiime_export(path, tmp_tsv, "TSVTaxonomyFormat")
            _write_headerless_taxonomy_tsv(tmp_tsv, dest)
            os.remove(tmp_tsv)

    return dest


def reference_dataframe(cfg: dict) -> pd.DataFrame:
    dbs = cfg.get("reference_databases", [])
    exclude = set(cfg.get("exclude_databases") or [])
    staging_root = join(data_dir(cfg), "ref_dbs")
    records = {}
    for db in dbs:
        db_id = db["id"]
        if db_id in exclude:
            continue
        db_staging = join(staging_root, db_id)
        refseqs = _ensure_text_reference(
            abspath(expandvars(db["refseqs_file"])),
            db_staging,
            "sequences",
        )
        taxa = _ensure_text_reference(
            abspath(expandvars(db["taxa_file"])),
            db_staging,
            "taxonomy",
        )
        records[db_id] = [
            refseqs,
            taxa,
            db_id,
            db["fwd_primer"],
            db["rev_primer"],
            db.get("fwd_primer_id", "F"),
            db.get("rev_primer_id", "R"),
        ]
    if not records:
        raise ValueError("No reference_databases configured (after exclude_databases).")
    df = pd.DataFrame.from_dict(records, orient="index")
    df.columns = [
        "Reference file path",
        "Reference tax path",
        "Reference id",
        "Fwd primer",
        "Rev primer",
        "Fwd primer id",
        "Rev primer id",
    ]
    return df


def simulation_methods_for_eval(evaluation_methods: list) -> list:
    mapping = {
        "cross-validated": "cross-validated-taxa",
        "cross-validated-taxa": "cross-validated-taxa",
        "cross-validated-trad": "cross-validated-trad",
        "novel-taxa": "novel-taxa",
    }
    methods = []
    for m in evaluation_methods:
        if m in mapping:
            methods.append(mapping[m])
        elif m not in ("mock-community", "self-validated"):
            raise ValueError(f"Unknown evaluation_method: {m}")
    return list(dict.fromkeys(methods))


def analysis_data_subdir(evaluation_method: str) -> str:
    if evaluation_method in ("cross-validated", "cross-validated-taxa"):
        return "cross-validated"
    if evaluation_method == "cross-validated-trad":
        return "cross-validated-trad"
    if evaluation_method == "novel-taxa":
        return "novel-taxa-simulations"
    if evaluation_method == "self-validated":
        return "self-validated"
    raise ValueError(evaluation_method)


def _cv_recall_levels(cfg: dict) -> tuple[int, int]:
    """Cross-validated manifest levels: one pass at max_level (min = max - 1)."""
    max_level = cfg.get("cv_recall_max_level", 6)
    if max_level < 1:
        raise ValueError("cv_recall_max_level must be >= 1")
    return max_level, max_level - 1


_SUPPORTED_CLASSIFY_METHODS = frozenset({
    "naive-bayes",
    "consensus-blast",
    "consensus-vsearch",
    "bt2-blca",
})

_SHARED_FIT_METHODS = frozenset({"naive-bayes", "bt2-blca"})


def classify_methods(cfg: dict) -> list[str]:
    """Return ordered list of taxonomy assignment methods to benchmark."""
    methods = cfg.get("classify_methods")
    if methods is None:
        single = cfg.get("classify_method", "naive-bayes")
        methods = [single] if isinstance(single, str) else list(single)
    elif isinstance(methods, str):
        methods = [methods]
    else:
        methods = list(methods)
    if not methods:
        raise ValueError("At least one classify method is required.")
    unknown = set(methods) - _SUPPORTED_CLASSIFY_METHODS
    if unknown:
        raise ValueError(
            f"Unsupported classify method(s): {sorted(unknown)}. "
            f"Supported: {sorted(_SUPPORTED_CLASSIFY_METHODS)}"
        )
    return methods


def _param_value_list(cfg: dict, key: str, default: float | None = None) -> list[float]:
    val = cfg.get(key, default)
    if val is None:
        raise ValueError(f"Config key {key!r} is required and has no value.")
    if isinstance(val, (list, tuple)):
        return [float(v) for v in val]
    return [float(val)]


def _confidence_value_list(cfg: dict, method: str) -> list[float]:
    """Return confidence thresholds to sweep for naive-bayes or bt2-blca."""
    if method == "naive-bayes":
        keys = ("nb_confidence_values", "confidence_values")
        fallback = cfg.get("skl_confidence", 0.7)
    elif method == "bt2-blca":
        keys = ("blca_confidence_values", "confidence_values")
        fallback = cfg.get("confidence_thres", 0.8)
    else:
        return []
    for key in keys:
        val = cfg.get(key)
        if val is None:
            continue
        if isinstance(val, (list, tuple)):
            return [float(v) for v in val]
        return [float(val)]
    return [float(fallback)]


def _manifest_na() -> str:
    """Placeholder for manifest fields not used by a job."""
    return ""


def consensus_param_combinations(cfg: dict) -> list[dict[str, float]]:
    """Cartesian product of consensus classifier parameters from config."""
    perc = _param_value_list(cfg, "perc_identity", 0.8)
    qc = _param_value_list(cfg, "query_cov", 0.8)
    mc = _param_value_list(cfg, "min_consensus", 0.51)
    return [
        {"perc_identity": pi, "query_cov": q, "min_consensus": m}
        for pi, q, m in itertools.product(perc, qc, mc)
    ]


def bt2_param_combinations(cfg: dict) -> list[dict[str, float]]:
    """Cartesian product of bt2-blca BLCA cutoff parameters from config.

    Read from ``blca_perc_identity`` / ``blca_query_cov``, separate from the
    consensus-method ``perc_identity`` / ``query_cov`` keys. Both are required.
    """
    validate_blca_config(cfg)
    perc = _param_value_list(cfg, "blca_perc_identity")
    qc = _param_value_list(cfg, "blca_query_cov")
    return [
        {"perc_identity": pi, "query_cov": q}
        for pi, q in itertools.product(perc, qc)
    ]


def _format_sweep_param(value: float) -> str:
    return format(float(value), "g")


def consensus_param_id(combo: dict[str, float]) -> str:
    return (
        f"pi{_format_sweep_param(combo['perc_identity'])}-"
        f"qc{_format_sweep_param(combo['query_cov'])}-"
        f"mc{_format_sweep_param(combo['min_consensus'])}"
    )


def bt2_param_id(combo: dict[str, float]) -> str:
    return (
        f"pi{_format_sweep_param(combo['perc_identity'])}-"
        f"qc{_format_sweep_param(combo['query_cov'])}"
    )


def method_settings(cfg: dict, method: str) -> dict:
    """Per-method assignment settings derived from config."""
    classify_params = cfg.get("classify_params") or ""
    if method == "naive-bayes":
        classify_params = cfg.get("naive_bayes_classify_params", classify_params)
    elif method in ("consensus-blast", "consensus-vsearch"):
        classify_params = cfg.get("consensus_classify_params", classify_params)
    settings = {
        "classify_method": method,
        "classify_threads": cfg.get("classify_threads", 5),
        "fit_params": cfg.get("fit_params") or "",
        "classify_params": classify_params or "",
    }
    if method == "bt2-blca":
        settings["taxa_ranks"] = cfg.get(
            "taxa_ranks", "kingdom,phylum,class,order,family,genus,species"
        )
        settings["perc_identity"] = cfg["blca_perc_identity"]
        settings["query_cov"] = cfg["blca_query_cov"]
    return settings


def confidences_for_method(cfg: dict, method: str) -> list[float]:
    return _confidence_value_list(cfg, method)


def param_id(
    method: str,
    method_cfg: dict,
    confidence: float,
    consensus_combo: dict[str, float] | None = None,
    bt2_combo: dict[str, float] | None = None,
) -> str:
    if method in ("consensus-blast", "consensus-vsearch"):
        if consensus_combo is None:
            raise ValueError("consensus_combo is required for consensus methods")
        return consensus_param_id(consensus_combo)
    if method == "bt2-blca":
        if bt2_combo is None:
            raise ValueError("bt2_combo is required for bt2-blca sweeps")
        base = bt2_param_id(bt2_combo)
        return f"{base}-conf{confidence}"
    fit_params = (method_cfg.get("fit_params") or "").strip()
    if method == "naive-bayes" and fit_params:
        digest = hashlib.md5(fit_params.encode()).hexdigest()[:8]
        return f"nb-{digest}-conf{confidence}"
    return f"{method}-conf{confidence}"


def fit_param_id(method: str, method_cfg: dict) -> str:
    """Shared classifier directory name (confidence-independent)."""
    fit_params = (method_cfg.get("fit_params") or "").strip()
    if method == "naive-bayes" and fit_params:
        digest = hashlib.md5(fit_params.encode()).hexdigest()[:8]
        return f"nb-{digest}"
    return method.replace("-", "_")


def _manifest_assign_params(
    method: str,
    method_cfg: dict,
    consensus_combo: dict[str, float] | None = None,
    bt2_combo: dict[str, float] | None = None,
) -> dict:
    na = _manifest_na()
    params = {
        "classify_method": method,
        "fit_params": na,
        "classify_params": na,
        "perc_identity": na,
        "query_cov": na,
        "min_consensus": na,
        "taxa_ranks": na,
    }
    if method == "naive-bayes":
        params["fit_params"] = method_cfg.get("fit_params") or na
        params["classify_params"] = method_cfg.get("classify_params") or na
    elif method in ("consensus-blast", "consensus-vsearch"):
        if consensus_combo is None:
            raise ValueError("consensus_combo is required for consensus methods")
        params["classify_params"] = method_cfg.get("classify_params") or na
        params["perc_identity"] = float(consensus_combo["perc_identity"])
        params["query_cov"] = float(consensus_combo["query_cov"])
        params["min_consensus"] = float(consensus_combo["min_consensus"])
    elif method == "bt2-blca":
        if bt2_combo is None:
            return params
        params["taxa_ranks"] = method_cfg.get(
            "taxa_ranks", "kingdom,phylum,class,order,family,genus,species"
        )
        params["perc_identity"] = float(bt2_combo["perc_identity"])
        params["query_cov"] = float(bt2_combo["query_cov"])
    else:
        raise ValueError(f"Unsupported classify method for manifest: {method}")
    return params


def _empty_manifest_artifact_fields() -> dict:
    return {
        "classifier_qza": "",
        "bowtie_index_dir": "",
    }


def _manifest_confidence(method: str, confidence: float | None) -> str | float:
    """Return confidence for manifest rows; blank when not used by the method."""
    if method in _SHARED_FIT_METHODS and confidence is not None:
        return confidence
    return _manifest_na()


def _append_fold_assignment_rows(
    rows: list,
    *,
    cfg: dict,
    results_root: str,
    subdir: str,
    dataset_id: str,
    reference_id: str,
    query: str,
    ref_seqs: str,
    ref_taxa: str,
    eval_method: str,
    method: str,
    method_cfg: dict,
    confidences: list,
    job_id_prefix: str,
) -> None:
    """Add manifest rows for one simulated fold (or mock-community combo)."""
    if method in ("consensus-blast", "consensus-vsearch"):
        for combo in consensus_param_combinations(cfg):
            assign_params = _manifest_assign_params(method, method_cfg, combo)
            p = param_id(method, method_cfg, 0.0, combo)
            assign_dir = join(results_root, subdir, dataset_id, reference_id, method, p)
            rows.append({
                "job_id": f"{job_id_prefix}-{p}",
                "evaluation_method": eval_method,
                "dataset_id": dataset_id,
                "reference_id": reference_id,
                "query_qza": query,
                "ref_seqs": ref_seqs,
                "ref_taxa": ref_taxa,
                "output_dir": assign_dir,
                "confidence": _manifest_na(),
                "skip_fit": False,
                "trad_fit": False,
                "fit_only": False,
                "fit_job_id": "",
                **_empty_manifest_artifact_fields(),
                **assign_params,
            })
        return

    if method == "bt2-blca":
        bt2_combos = bt2_param_combinations(cfg)
        multi_shared = len(confidences) > 1 or len(bt2_combos) > 1
        if multi_shared:
            shared_dir = join(
                results_root, subdir, dataset_id, reference_id, method,
                fit_param_id(method, method_cfg),
            )
            fit_job_id = f"{job_id_prefix}-{fit_param_id(method, method_cfg)}-fit"
            rows.append({
                "job_id": fit_job_id,
                "evaluation_method": eval_method,
                "dataset_id": dataset_id,
                "reference_id": reference_id,
                "query_qza": query,
                "ref_seqs": ref_seqs,
                "ref_taxa": ref_taxa,
                "output_dir": shared_dir,
                "confidence": _manifest_na(),
                "skip_fit": False,
                "trad_fit": False,
                "fit_only": True,
                "fit_job_id": "",
                **_empty_manifest_artifact_fields(),
                **_manifest_assign_params(method, method_cfg),
            })
            shared_artifacts = {
                **_empty_manifest_artifact_fields(),
                "bowtie_index_dir": join(shared_dir, "bowtie2_index"),
            }
            for combo in bt2_combos:
                for conf in confidences:
                    p = param_id(method, method_cfg, conf, bt2_combo=combo)
                    assign_dir = join(
                        results_root, subdir, dataset_id, reference_id, method, p
                    )
                    rows.append({
                        "job_id": f"{job_id_prefix}-{p}",
                        "evaluation_method": eval_method,
                        "dataset_id": dataset_id,
                        "reference_id": reference_id,
                        "query_qza": query,
                        "ref_seqs": ref_seqs,
                        "ref_taxa": ref_taxa,
                        "output_dir": assign_dir,
                        "confidence": _manifest_confidence(method, conf),
                        "skip_fit": True,
                        "trad_fit": False,
                        "fit_only": False,
                        "fit_job_id": fit_job_id,
                        **shared_artifacts,
                        **_manifest_assign_params(method, method_cfg, bt2_combo=combo),
                    })
            return
        combo = bt2_combos[0]
        assign_params = _manifest_assign_params(method, method_cfg, bt2_combo=combo)
        for conf in confidences:
            p = param_id(method, method_cfg, conf, bt2_combo=combo)
            assign_dir = join(results_root, subdir, dataset_id, reference_id, method, p)
            rows.append({
                "job_id": f"{job_id_prefix}-{p}",
                "evaluation_method": eval_method,
                "dataset_id": dataset_id,
                "reference_id": reference_id,
                "query_qza": query,
                "ref_seqs": ref_seqs,
                "ref_taxa": ref_taxa,
                "output_dir": assign_dir,
                "confidence": _manifest_confidence(method, conf),
                "skip_fit": False,
                "trad_fit": False,
                "fit_only": False,
                "fit_job_id": "",
                **_empty_manifest_artifact_fields(),
                **assign_params,
            })
        return

    assign_params = _manifest_assign_params(method, method_cfg)
    multi_conf_shared = method in _SHARED_FIT_METHODS and len(confidences) > 1

    if multi_conf_shared:
        shared_dir = join(
            results_root, subdir, dataset_id, reference_id, method, fit_param_id(method, method_cfg)
        )
        fit_job_id = f"{job_id_prefix}-{fit_param_id(method, method_cfg)}-fit"
        rows.append({
            "job_id": fit_job_id,
            "evaluation_method": eval_method,
            "dataset_id": dataset_id,
            "reference_id": reference_id,
            "query_qza": query,
            "ref_seqs": ref_seqs,
            "ref_taxa": ref_taxa,
            "output_dir": shared_dir,
            "confidence": _manifest_na(),
            "skip_fit": False,
            "trad_fit": False,
            "fit_only": True,
            "fit_job_id": "",
            **_empty_manifest_artifact_fields(),
            **assign_params,
        })
        shared_artifacts = _empty_manifest_artifact_fields()
        if method == "naive-bayes":
            shared_artifacts["classifier_qza"] = join(shared_dir, "classifier.qza")
        else:
            shared_artifacts["bowtie_index_dir"] = join(shared_dir, "bowtie2_index")
        for conf in confidences:
            p = param_id(method, method_cfg, conf)
            assign_dir = join(
                results_root, subdir, dataset_id, reference_id, method, p
            )
            rows.append({
                "job_id": f"{job_id_prefix}-{p}",
                "evaluation_method": eval_method,
                "dataset_id": dataset_id,
                "reference_id": reference_id,
                "query_qza": query,
                "ref_seqs": ref_seqs,
                "ref_taxa": ref_taxa,
                "output_dir": assign_dir,
                "confidence": _manifest_confidence(method, conf),
                "skip_fit": True,
                "trad_fit": False,
                "fit_only": False,
                "fit_job_id": fit_job_id,
                **shared_artifacts,
                **assign_params,
            })
        return

    for conf in confidences:
        p = param_id(method, method_cfg, conf)
        assign_dir = join(results_root, subdir, dataset_id, reference_id, method, p)
        rows.append({
            "job_id": f"{job_id_prefix}-{p}",
            "evaluation_method": eval_method,
            "dataset_id": dataset_id,
            "reference_id": reference_id,
            "query_qza": query,
            "ref_seqs": ref_seqs,
            "ref_taxa": ref_taxa,
            "output_dir": assign_dir,
            "confidence": _manifest_confidence(method, conf),
            "skip_fit": False,
            "trad_fit": False,
            "fit_only": False,
            "fit_job_id": "",
            **_empty_manifest_artifact_fields(),
            **assign_params,
        })


def prepare_datasets(cfg: dict) -> None:
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.framework_functions import (
        generate_self_validated_datasets,
        generate_simulated_datasets,
    )

    out = run_output_dir(cfg)
    os.makedirs(out, exist_ok=True)
    ddir = data_dir(cfg)
    os.makedirs(ddir, exist_ok=True)

    df = reference_dataframe(cfg)
    eval_methods = cfg.get("evaluation_methods", [])
    sim_methods = simulation_methods_for_eval(eval_methods)

    needs_generation = (
        bool(sim_methods)
        or "self-validated" in eval_methods
        or "mock-community" in eval_methods
    )
    if not needs_generation:
        raise ValueError("No evaluation_methods require dataset generation.")

    if sim_methods:
        levelrange = cfg.get("novel_taxa_levels", [6, 5, 4, 3])
        force = cfg.get("force_regenerate", False)
        exclude = set(cfg.get("exclude_databases") or [])
        for db in cfg.get("reference_databases", []):
            db_id = db["id"]
            if db_id in exclude or db_id not in df.index:
                continue
            sim_params = _db_simulation_params(db)
            generate_simulated_datasets(
                df.loc[[db_id]],
                ddir,
                cfg["iterations"],
                levelrange=levelrange,
                force=force,
                simulation_method=sim_methods,
                **sim_params,
            )

    if "self-validated" in eval_methods:
        force = cfg.get("force_regenerate", False)
        exclude = set(cfg.get("exclude_databases") or [])
        for db in cfg.get("reference_databases", []):
            db_id = db["id"]
            if db_id in exclude or db_id not in df.index:
                continue
            sim_params = _db_simulation_params(db)
            generate_self_validated_datasets(
                df.loc[[db_id]],
                ddir,
                force=force,
                **sim_params,
            )

    if "mock-community" in eval_methods:
        stage_mock_communities(cfg)


def stage_mock_communities(cfg: dict) -> None:
    mock_root = cfg.get("mock_dir") or join(data_dir(cfg), "mock-community")
    os.makedirs(mock_root, exist_ok=True)
    for mock in cfg.get("mock_communities") or []:
        mock_id = mock["id"]
        dest = join(mock_root, mock_id)
        os.makedirs(dest, exist_ok=True)
        if mock.get("feature_table_biom"):
            shutil.copy2(
                abspath(expandvars(mock["feature_table_biom"])),
                join(dest, "feature_table.biom"),
            )
        if mock.get("rep_seqs_fasta"):
            rep = join(dest, "rep_seqs.fna")
            shutil.copy2(abspath(expandvars(mock["rep_seqs_fasta"])), rep)
        for ref in mock.get("references") or []:
            ref_dest = join(dest, ref["id"], "expected")
            os.makedirs(ref_dest, exist_ok=True)
            expected = abspath(expandvars(ref["expected_dir"]))
            for name in os.listdir(expected):
                shutil.copy2(join(expected, name), join(ref_dest, name))


def prepare_manifest(cfg: dict) -> str:
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.framework_functions import (
        recall_self_validated_dirs,
        recall_simulated_taxa_dirs,
    )

    out = run_output_dir(cfg)
    manifest_fp = join(out, "assignment_manifest.tsv")
    rows = []
    ddir = data_dir(cfg)
    results_root = results_dir(cfg)
    os.makedirs(results_root, exist_ok=True)
    db_ids = list(reference_dataframe(cfg).index)
    cv_max_level, cv_min_level = _cv_recall_levels(cfg)

    for classify_method in classify_methods(cfg):
        method_cfg = method_settings(cfg, classify_method)
        confidences = confidences_for_method(cfg, classify_method)

        for eval_method in cfg.get("evaluation_methods", []):
            if eval_method == "mock-community":
                mock_root = cfg.get("mock_dir") or join(ddir, "mock-community")
                for mock in cfg.get("mock_communities") or []:
                    mock_id = mock["id"]
                    query = join(mock_root, mock_id, "rep_seqs.qza")
                    for ref in mock.get("references") or []:
                        ref_id = ref["id"]
                        _append_fold_assignment_rows(
                            rows,
                            cfg=cfg,
                            results_root=results_root,
                            subdir="mock-community",
                            dataset_id=mock_id,
                            reference_id=ref_id,
                            query=query,
                            ref_seqs=join(ddir, "ref_dbs", ref_id, "ref_seqs.qza"),
                            ref_taxa=join(ddir, "ref_dbs", ref_id, "ref_taxa.qza"),
                            eval_method=eval_method,
                            method=classify_method,
                            method_cfg=method_cfg,
                            confidences=confidences,
                            job_id_prefix=f"mock-{mock_id}-{ref_id}-{classify_method}",
                        )
                continue

            subdir = analysis_data_subdir(eval_method)
            sim_dir = join(ddir, subdir)
            if eval_method in ("cross-validated", "cross-validated-taxa"):
                combos, ref_dbs = recall_simulated_taxa_dirs(
                    sim_dir, db_ids, cfg["iterations"],
                    ref_seqs="ref_seqs.qza", ref_taxa="ref_taxa.qza",
                    max_level=cv_max_level, min_level=cv_min_level,
                    multilevel=False,
                )
            elif eval_method == "novel-taxa":
                combos, ref_dbs = recall_simulated_taxa_dirs(
                    sim_dir, db_ids, cfg["iterations"],
                    ref_seqs="ref_seqs.qza", ref_taxa="ref_taxa.qza",
                    max_level=6, min_level=cfg.get("novel_recall_min_level", 3),
                    multilevel=True,
                )
            elif eval_method == "cross-validated-trad":
                from tax_credit.framework_functions import trad_cv_shared_reference_qzas

                combos, ref_dbs = recall_simulated_taxa_dirs(
                    sim_dir, db_ids, cfg["iterations"],
                    ref_seqs="ref_seqs.qza", ref_taxa="ref_taxa.qza",
                    max_level=cv_max_level, min_level=cv_min_level,
                    multilevel=False,
                )
                fit_id = fit_param_id(classify_method, method_cfg)

                if classify_method in _SHARED_FIT_METHODS:
                    assign_params = _manifest_assign_params(classify_method, method_cfg)
                    for db_id in db_ids:
                        ref_seqs, ref_taxa = trad_cv_shared_reference_qzas(ddir, db_id)
                        shared_dir = join(
                            results_root, "trad-fit", db_id, classify_method, fit_id
                        )
                        rows.append({
                            "job_id": f"trad-fit-{db_id}-{classify_method}-{fit_id}",
                            "evaluation_method": eval_method,
                            "dataset_id": db_id,
                            "reference_id": db_id,
                            "query_qza": "",
                            "ref_seqs": ref_seqs,
                            "ref_taxa": ref_taxa,
                            "output_dir": shared_dir,
                            "confidence": _manifest_na(),
                            "skip_fit": False,
                            "trad_fit": True,
                            "fit_only": False,
                            "fit_job_id": "",
                            **_empty_manifest_artifact_fields(),
                            **assign_params,
                        })
                    for dataset_id, reference_id in combos:
                        fold_dir = join(sim_dir, dataset_id)
                        query = join(fold_dir, "query.qza")
                        shared_artifacts = _empty_manifest_artifact_fields()
                        if classify_method == "naive-bayes":
                            shared_artifacts["classifier_qza"] = join(
                                results_root,
                                "trad-fit",
                                reference_id,
                                classify_method,
                                fit_id,
                                "classifier.qza",
                            )
                        else:
                            shared_artifacts["bowtie_index_dir"] = join(
                                results_root,
                                "trad-fit",
                                reference_id,
                                classify_method,
                                fit_id,
                                "bowtie2_index",
                            )
                        for conf in confidences:
                            p = param_id(classify_method, method_cfg, conf)
                            assign_dir = join(
                                results_root,
                                subdir,
                                dataset_id,
                                reference_id,
                                classify_method,
                                p,
                            )
                            rows.append({
                                "job_id": f"trad-{dataset_id}-{classify_method}-{p}",
                                "evaluation_method": eval_method,
                                "dataset_id": dataset_id,
                                "reference_id": reference_id,
                                "query_qza": query,
                                "ref_seqs": ref_dbs[dataset_id][0],
                                "ref_taxa": ref_dbs[dataset_id][1],
                                "output_dir": assign_dir,
                                "confidence": _manifest_confidence(classify_method, conf),
                                "skip_fit": True,
                                "trad_fit": False,
                                "fit_only": False,
                                "fit_job_id": "",
                                **shared_artifacts,
                                **assign_params,
                            })
                else:
                    for dataset_id, reference_id in combos:
                        ref_seqs, ref_taxa = ref_dbs[dataset_id]
                        query = join(sim_dir, dataset_id, "query.qza")
                        _append_fold_assignment_rows(
                            rows,
                            cfg=cfg,
                            results_root=results_root,
                            subdir=subdir,
                            dataset_id=dataset_id,
                            reference_id=reference_id,
                            query=query,
                            ref_seqs=ref_seqs,
                            ref_taxa=ref_taxa,
                            eval_method=eval_method,
                            method=classify_method,
                            method_cfg=method_cfg,
                            confidences=confidences,
                            job_id_prefix=f"trad-{dataset_id}-{classify_method}",
                        )
                continue
            elif eval_method == "self-validated":
                combos, ref_dbs = recall_self_validated_dirs(
                    sim_dir, db_ids,
                    ref_seqs="ref_seqs.qza", ref_taxa="ref_taxa.qza",
                )
            else:
                raise ValueError(f"Unknown evaluation_method: {eval_method}")

            for dataset_id, reference_id in combos:
                ref_seqs, ref_taxa = ref_dbs[dataset_id]
                query = join(sim_dir, dataset_id, "query.qza")
                _append_fold_assignment_rows(
                    rows,
                    cfg=cfg,
                    results_root=results_root,
                    subdir=subdir,
                    dataset_id=dataset_id,
                    reference_id=reference_id,
                    query=query,
                    ref_seqs=ref_seqs,
                    ref_taxa=ref_taxa,
                    eval_method=eval_method,
                    method=classify_method,
                    method_cfg=method_cfg,
                    confidences=confidences,
                    job_id_prefix=f"{subdir}-{dataset_id}-{classify_method}",
                )

    manifest = pd.DataFrame(rows)
    manifest.to_csv(manifest_fp, sep="\t", index=False)
    return manifest_fp


def evaluate(cfg: dict) -> None:
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.framework_functions import novel_taxa_classification_evaluation
    from tax_credit.mock_evaluation import evaluate_results
    from tax_credit.paths import list_assignment_result_dirs

    out = run_output_dir(cfg)
    summaries_dir = join(out, "summaries")
    os.makedirs(summaries_dir, exist_ok=True)
    ddir = data_dir(cfg)
    results_root = results_dir(cfg)
    force = cfg.get("force_evaluation", False)
    summary_names = cfg.get("summary_filenames") or {}

    for eval_method in cfg.get("evaluation_methods", []):
        if eval_method == "mock-community":
            mock_root = cfg.get("mock_dir") or join(ddir, "mock-community")
            results_dirs = list_assignment_result_dirs(
                join(results_root, "mock-community")
            )
            summary_fp = join(
                summaries_dir,
                summary_names.get("mock-community", "mock_evaluation_summary.tsv"),
            )
            evaluate_results(
                results_dirs,
                expected_results_dir=mock_root,
                results_fp=summary_fp,
                mock_dir=mock_root,
                taxonomy_level_range=cfg.get("mock_taxonomy_level_range", range(2, 7)),
                per_seq_precision=cfg.get("mock_per_seq_precision", False),
                force=force,
            )
            continue

        sub = analysis_data_subdir(eval_method)
        computed = join(results_root, sub)
        expected = join(ddir, sub)
        results_dirs = list_assignment_result_dirs(computed)
        if not results_dirs:
            print(f"No assignment results under {computed}; skipping evaluation.")
            continue

        default_name = {
            "cross-validated": "evaluate_classification_summary_CV.csv",
            "cross-validated-taxa": "evaluate_classification_summary_CV.csv",
            "novel-taxa": "evaluate_classification_summary_novel.csv",
            "cross-validated-trad": "evaluate_classification_summary_CV_trad.csv",
            "self-validated": "evaluate_classification_summary_self_validated.csv",
        }.get(eval_method, f"evaluate_{sub}.csv")
        summary_fp = join(summaries_dir, summary_names.get(eval_method, default_name))

        test_type = "novel-taxa"
        if eval_method in ("cross-validated", "cross-validated-taxa"):
            test_type = "cross-validated"
        elif eval_method == "cross-validated-trad":
            test_type = "cross-validated-trad"
        elif eval_method == "self-validated":
            test_type = "self-validated"

        novel_taxa_classification_evaluation(
            results_dirs,
            expected,
            summary_fp,
            test_type=test_type,
        )


def _summary_has_per_level_lists(df: pd.DataFrame, column: str = "Precision") -> bool:
    """True when summary metrics are stored as per-rank lists (CV-style)."""
    if column not in df.columns or df.empty:
        return False
    sample = df[column].iloc[0]
    if isinstance(sample, str):
        return sample.strip().startswith("[")
    return isinstance(sample, (list, tuple))


def _rank_name(level: int) -> str:
    from tax_credit.log_analysis import RANK_NAMES

    return RANK_NAMES[int(level)]


def _per_level_plot_table(
    summary_df: pd.DataFrame, assignment_dirs: list[str],
) -> pd.DataFrame:
    """One row per Dataset / iteration / Method / Parameters / level.

    Precision, Recall and F-measure come from the summary's per-level lists; the
    four classification ratios come from the per-read logs. Every metric is
    therefore reported for the same fold and taxonomic level.
    """
    from tax_credit.novel_evaluation import (
        extract_per_level_accuracy,
        extract_per_level_classification_ratios_by_fold,
    )

    keys = ["Dataset", "iteration", "Method", "Parameters", "level"]
    table = extract_per_level_accuracy(
        summary_df, columns=["Precision", "Recall", "F-measure"],
    )
    table["iteration"] = table["iteration"].astype(str)
    table["level"] = table["level"].astype(int)
    ratios = extract_per_level_classification_ratios_by_fold(assignment_dirs)
    if not ratios.empty:
        ratios = ratios.drop(columns="novel_level")
        ratios["iteration"] = ratios["iteration"].astype(str)
        ratios["level"] = ratios["level"].astype(int)
        table = table.merge(ratios, on=keys, how="left")
    table["rank"] = table["level"].map(_rank_name)
    return table


def _novel_plot_table(summary_df: pd.DataFrame) -> pd.DataFrame:
    """One row per novel-taxa fold; ``novel_level`` labels the simulation (L5, L6, ...)."""
    table = summary_df.drop(columns=["mismatch_level_list"], errors="ignore").copy()
    table["novel_level"] = "L" + table["level"].astype(int).astype(str)
    return table


def _plot_ranks(cfg: dict) -> list[str]:
    """Ranks for per-rank boxplots (config ``plot_ranks``)."""
    from tax_credit.log_analysis import RANK_TO_LEVEL

    ranks = cfg.get("plot_ranks") or ["genus", "species"]
    if isinstance(ranks, str):
        ranks = [ranks]
    ranks = [str(rank).strip().lower() for rank in ranks]
    unknown = [rank for rank in ranks if rank not in RANK_TO_LEVEL]
    if unknown:
        raise ValueError(
            f"Unknown plot_ranks entries: {unknown}. Use one of: {list(RANK_TO_LEVEL)}"
        )
    return ranks


def _best_run_rank(cfg: dict) -> str:
    """Rank at which best_run_stacked_bar compares runs (config ``best_run_rank``)."""
    from tax_credit.log_analysis import RANK_TO_LEVEL

    rank = str(cfg.get("best_run_rank") or "species").strip().lower()
    if rank not in RANK_TO_LEVEL:
        raise ValueError(
            f"Unknown best_run_rank: {rank!r}. Use one of: {list(RANK_TO_LEVEL)}"
        )
    return rank


def _metric_label(metric: str) -> str:
    """Readable metric name, e.g. ``misclassification_ratio`` -> ``misclassification``."""
    return metric.replace("_ratio", "").replace("_", " ")


def _plot_best_runs(
    table: pd.DataFrame,
    ratio_df: pd.DataFrame,
    level_col: str,
    metrics: list[str],
    best_run_rank: str,
    plot_label: str,
    base: str,
    plots_out: str,
    csv_out: str,
    rank_labels: dict,
    rank_axis_label: str,
) -> None:
    """Plot classification ratios for the best method + parameters per metric.

    For each database, picks the best run for every metric (averaged over
    folds) and draws its ratios by rank. One figure holds every database in the
    evaluation method: databases (and novel levels) are rows, metrics columns.
    Cross-validated and self-validated runs are compared at *best_run_rank*;
    novel-taxa runs are compared per novel level. Selections are written to
    *csv_out*.
    """
    from tax_credit.log_analysis import RANK_TO_LEVEL
    from tax_credit.novel_evaluation import select_best_runs
    from tax_credit.plotting_functions import (
        stacked_classification_panels_from_data_frames,
    )
    import matplotlib.pyplot as plt

    novel = level_col == "novel_level"
    if novel:
        group_cols = ["Dataset", "novel_level"]
        select_df = table
    else:
        group_cols = ["Dataset"]
        select_df = table[table["level"] == RANK_TO_LEVEL[best_run_rank]]
    best = select_best_runs(select_df, metrics, group_cols=group_cols)
    if best.empty or ratio_df.empty:
        return
    if not novel:
        best.insert(1, "rank", best_run_rank)
    os.makedirs(dirname(csv_out), exist_ok=True)
    best.to_csv(csv_out, index=False)

    # one row per database (per novel level), one column per metric
    panels = []
    row_labels = []
    for group_key, picks in best.groupby(group_cols, sort=True):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        dataset = group_key[0]
        runs = ratio_df[ratio_df["Dataset"] == dataset]
        if novel:
            novel_level = int(group_key[1][1:])
            # only ranks above the novel rank have a meaningful expected taxonomy
            runs = runs[(runs["novel_level"] == novel_level) & (runs["level"] < novel_level)]
            row_labels.append(f"{dataset}\nL{novel_level}")
        else:
            row_labels.append(dataset)

        picked = {pick.metric: pick for pick in picks.itertuples(index=False)}
        for metric in metrics:
            pick = picked.get(metric)
            if pick is None:
                panels.append((f"{_metric_label(metric)}: no data", None))
                continue
            run = runs[(runs["Method"] == pick.Method) & (runs["Parameters"] == pick.Parameters)]
            others = pick.n_tied - 1
            ties = f"\n(tied with {others} other run{'s' if others > 1 else ''})" if others else ""
            panels.append((
                f"best {_metric_label(pick.metric)} ({pick.direction}): {pick.value:.3f}\n"
                f"{pick.Method}\n{pick.Parameters}{ties}",
                run,
            ))

    if not panels:
        return
    selected_at = (
        "each novel level's scores" if novel else f"{best_run_rank} rank"
    )
    fig = stacked_classification_panels_from_data_frames(
        panels,
        level_labels=rank_labels,
        level_axis_label=rank_axis_label,
        ncols=max(len(metrics), 1),
        row_labels=row_labels,
        title=(
            f"{plot_label}: best method and parameters per metric "
            f"(selected at {selected_at})"
        ),
        show=False,
    )
    fig.savefig(
        join(plots_out, f"{base}-best-run-stacked-barplot.pdf"),
        bbox_inches="tight",
    )
    plt.close(fig)


def _plot_type_set(cfg: dict) -> set[str]:
    """Normalize configured plot types (boxplot, pointplot, heatmap, stacked_bar,
    best_run_stacked_bar)."""
    aliases = {
        "boxplot": "boxplot",
        "box": "boxplot",
        "pointplot": "pointplot",
        "point": "pointplot",
        "line": "pointplot",
        "lineplot": "pointplot",
        "heatmap": "heatmap",
        "stacked_bar": "stacked_bar",
        "stacked-bar": "stacked_bar",
        "stackedbar": "stacked_bar",
        "stacked_barplot": "stacked_bar",
        "stacked-barplot": "stacked_bar",
        "best_run_stacked_bar": "best_run_stacked_bar",
        "best-run-stacked-bar": "best_run_stacked_bar",
    }
    configured = cfg.get("plot_types") or ["boxplot", "pointplot", "heatmap"]
    if isinstance(configured, str):
        configured = [configured]
    types = set()
    for item in configured:
        key = str(item).strip().lower()
        if key not in aliases:
            raise ValueError(
                f"Unknown plot_types entry: {item!r}. "
                f"Use one of: {sorted(set(aliases))}"
            )
        types.add(aliases[key])
    return types


def _eval_method_plot_label(eval_method: str) -> str:
    if eval_method in ("cross-validated", "cross-validated-taxa"):
        return "cross-validated"
    if eval_method == "novel-taxa":
        return "novel"
    return eval_method


def _log_analysis_ranks(cfg: dict) -> list[str]:
    ranks = cfg.get("log_analysis_ranks")
    if ranks is None:
        return [cfg.get("log_analysis_rank", "species")]
    if isinstance(ranks, str):
        return [ranks]
    return [str(rank) for rank in ranks]


def analyze_log_results(cfg: dict) -> None:
    if not cfg.get("generate_log_analysis", True):
        return
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.log_analysis import (
        RANK_NAMES,
        RANK_TO_LEVEL,
        filter_sensitivity_to_rank,
        load_classification_accuracy_logs,
        summarize_confusion_pairs,
        summarize_cross_fold_stability,
        summarize_method_parameter_sensitivity,
        summarize_taxon_errors,
        select_top_sensitivity_taxa,
    )
    from tax_credit.paths import list_assignment_result_dirs
    from tax_credit.log_plotting import (
        method_parameter_sensitivity_heatmap_from_data_frame,
    )
    import matplotlib.pyplot as plt

    rank = cfg.get("log_analysis_rank", "species")
    ranks = _log_analysis_ranks(cfg)
    min_obvs = int(
        cfg.get("log_analysis_min_obvs", cfg.get("log_analysis_min_reads", 3))
    )
    top_n = int(cfg.get("log_analysis_top_n", 25))
    out = run_output_dir(cfg)
    summaries_root = join(out, "summaries", "log_analysis")
    plots_dir = join(out, cfg.get("plots_subdir", "plots"))

    for eval_method in cfg.get("evaluation_methods", []):
        if eval_method == "mock-community":
            continue
        plot_label = _eval_method_plot_label(eval_method)
        sub = analysis_data_subdir(eval_method)
        computed = join(results_dir(cfg), sub)
        assignment_dirs = list_assignment_result_dirs(computed)
        log_df = load_classification_accuracy_logs(assignment_dirs)
        if log_df.empty:
            continue

        method_summary_dir = join(summaries_root, eval_method)
        method_plots_dir = join(plots_dir, eval_method)
        os.makedirs(method_summary_dir, exist_ok=True)
        os.makedirs(method_plots_dir, exist_ok=True)

        taxon_summary = summarize_taxon_errors(
            log_df, rank=rank, min_obvs=min_obvs,
        )
        taxon_summary.to_csv(
            join(method_summary_dir, "taxon_error_profiles.csv"), index=False,
        )

        confusion_pairs = summarize_confusion_pairs(log_df, rank=rank)
        confusion_pairs.to_csv(
            join(method_summary_dir, "confusion_pairs.csv"), index=False,
        )

        cross_fold = summarize_cross_fold_stability(
            log_df, rank=rank, min_obvs=min_obvs,
        )
        if not cross_fold.empty:
            cross_fold.to_csv(
                join(method_summary_dir, "cross_fold_stability.csv"), index=False,
            )

        # (rank, log rows, title label, file suffix) for each sensitivity table
        if eval_method == "novel-taxa":
            # Expected taxonomies stop above the novel rank, so each novel level
            # is analysed at the deepest rank its expected taxonomy reaches.
            panels = []
            for novel_level in sorted(log_df["level"].dropna().astype(int).unique()):
                parent_rank = RANK_NAMES[novel_level - 1]
                if parent_rank not in RANK_TO_LEVEL:
                    print(
                        f"WARNING: no sensitivity analysis for {plot_label} "
                        f"L{novel_level}: {parent_rank} rank is not supported",
                        file=sys.stderr,
                    )
                    continue
                panels.append((
                    parent_rank,
                    log_df[log_df["level"] == novel_level],
                    f"{plot_label} L{novel_level}",
                    f"L{novel_level}-{parent_rank}",
                ))
        else:
            panels = [(rank_, log_df, plot_label, rank_) for rank_ in ranks]

        for sensitivity_rank, panel_df, label, suffix in panels:
            sensitivity = summarize_method_parameter_sensitivity(
                panel_df,
                rank=sensitivity_rank,
                min_obvs=min_obvs,
                metric="pct_mis",
            )
            if sensitivity.empty:
                print(
                    f"WARNING: no {sensitivity_rank}-rank sensitivity data for {label}",
                    file=sys.stderr,
                )
                continue
            sensitivity.to_csv(
                join(method_summary_dir, f"method_parameter_sensitivity-{suffix}.csv"),
            )
            at_rank, n_excluded = filter_sensitivity_to_rank(sensitivity, sensitivity_rank)
            if at_rank.empty:
                print(
                    f"WARNING: no expected taxa resolved to {sensitivity_rank} rank "
                    f"for {label}; skipping sensitivity plot",
                    file=sys.stderr,
                )
                continue
            plot_pivot = select_top_sensitivity_taxa(at_rank, top_n=top_n)
            title = (
                f"{label}: fraction of reads misclassified at {sensitivity_rank} "
                f"rank (top {len(plot_pivot)} taxa)"
            )
            if n_excluded:
                title += f"\n{n_excluded} taxa not resolved to {sensitivity_rank} excluded"
            ax = method_parameter_sensitivity_heatmap_from_data_frame(
                plot_pivot, title=title, show=False,
            )
            ax.figure.savefig(
                join(method_plots_dir, f"method-parameter-sensitivity-{suffix}.pdf"),
                bbox_inches="tight",
            )
            plt.close(ax.figure)


def plot_results(cfg: dict) -> None:
    if not cfg.get("generate_plots", True):
        return
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.plotting_functions import (
        faceted_boxplot_from_data_frame,
        heatmap_from_data_frame,
        pointplot_from_data_frame,
        stacked_classification_barplot_from_data_frame,
    )
    from tax_credit.novel_evaluation import extract_per_level_classification_ratios
    from tax_credit.paths import list_assignment_result_dirs
    import matplotlib.pyplot as plt

    out = run_output_dir(cfg)
    plots_dir = join(out, cfg.get("plots_subdir", "plots"))
    os.makedirs(plots_dir, exist_ok=True)
    summaries_dir = join(out, "summaries")
    summary_names = cfg.get("summary_filenames") or {}
    plot_metrics = cfg.get(
        "plot_metrics", ["Precision", "Recall", "F-measure"]
    )
    plot_types = _plot_type_set(cfg)
    color_palette = cfg.get("plot_color_palette") or "tab10"
    heatmap_rows = cfg.get("plot_heatmap_rows") or ["Method", "Parameters"]
    heatmap_cols = cfg.get("plot_heatmap_cols") or ["Dataset", "level"]
    plot_ranks = _plot_ranks(cfg)
    best_run_rank = _best_run_rank(cfg)
    # stacked bars are narrow, so label them with rank initials
    stacked_rank_labels = {level: _rank_name(level)[0].upper() for level in range(1, 7)}
    stacked_rank_axis_label = "rank (" + ", ".join(
        f"{_rank_name(level)[0].upper()} = {_rank_name(level)}" for level in range(1, 7)
    ) + ")"
    default_summary = {
        "cross-validated": "evaluate_classification_summary_CV.csv",
        "cross-validated-taxa": "evaluate_classification_summary_CV.csv",
        "novel-taxa": "evaluate_classification_summary_novel.csv",
        "cross-validated-trad": "evaluate_classification_summary_CV_trad.csv",
        "self-validated": "evaluate_classification_summary_self_validated.csv",
    }

    for eval_method in cfg.get("evaluation_methods", []):
        if eval_method == "mock-community":
            continue
        plot_label = _eval_method_plot_label(eval_method)
        method_plots_dir = join(plots_dir, eval_method)
        os.makedirs(method_plots_dir, exist_ok=True)
        sub = analysis_data_subdir(eval_method)
        assignment_dirs = list_assignment_result_dirs(join(results_dir(cfg), sub))

        ratio_df = pd.DataFrame()
        if plot_types & {"stacked_bar", "best_run_stacked_bar"}:
            ratio_df = extract_per_level_classification_ratios(assignment_dirs)

        if "stacked_bar" in plot_types:
            stacked_df = ratio_df
            if not stacked_df.empty:
                novel_levels = sorted(stacked_df["novel_level"].dropna().unique())
                if novel_levels:
                    # one plot per novel level, only ranks above the novel rank
                    panels = [
                        (
                            f"{plot_label} L{novel_level}",
                            f"-L{novel_level}",
                            stacked_df[
                                (stacked_df["novel_level"] == novel_level)
                                & (stacked_df["level"] < novel_level)
                            ],
                        )
                        for novel_level in novel_levels
                    ]
                else:
                    panels = [(plot_label, "", stacked_df)]
                for label, suffix, panel_df in panels:
                    if panel_df.empty:
                        continue
                    ax = stacked_classification_barplot_from_data_frame(
                        panel_df,
                        title=f"{label}: classification ratios by rank",
                        level_labels=stacked_rank_labels,
                        level_axis_label=stacked_rank_axis_label,
                        show=False,
                    )
                    ax.figure.savefig(
                        join(
                            method_plots_dir,
                            f"classification-ratios-stacked-barplot{suffix}.pdf",
                        ),
                        bbox_inches="tight",
                    )
                    plt.close(ax.figure)

        candidates = [
            join(summaries_dir, summary_names.get(
                eval_method, default_summary.get(eval_method, f"evaluate_{sub}.csv")
            ))
        ]
        candidates.extend(glob(join(summaries_dir, f"*{sub}*")))
        if eval_method in ("cross-validated", "cross-validated-taxa"):
            candidates.extend(glob(join(summaries_dir, "*CV*")))

        seen = set()
        for fp in candidates:
            if fp in seen or not exists(fp) or not fp.endswith((".csv", ".tsv")):
                continue
            seen.add(fp)
            df = pd.read_csv(fp, index_col=0)
            if df.empty:
                continue
            base = Path(fp).stem
            if _summary_has_per_level_lists(df):
                # cross-validated / self-validated: metrics per fold and rank
                table = _per_level_plot_table(df, assignment_dirs)
                level_col = "rank"
                level_order = [
                    _rank_name(level) for level in sorted(table["level"].unique())
                ]
                box_panels = [rank_ for rank_ in plot_ranks if rank_ in level_order]
            else:
                # novel-taxa: metrics per fold and novel level
                table = _novel_plot_table(df)
                level_col = "novel_level"
                level_order = sorted(
                    table["novel_level"].unique(), key=lambda v: int(v[1:])
                )
                box_panels = level_order

            def _has_metric(metric):
                return metric in table.columns and table[metric].notna().any()

            if "boxplot" in plot_types and box_panels:
                box_df = table[table[level_col].isin(box_panels)]
                for metric in plot_metrics:
                    if not _has_metric(metric):
                        continue
                    grid = faceted_boxplot_from_data_frame(
                        box_df,
                        x="Dataset",
                        metric=metric,
                        hue="Method",
                        col=level_col,
                        col_order=box_panels,
                        color_palette=color_palette,
                        title=f"{plot_label}: {metric}",
                        show=False,
                    )
                    grid.savefig(
                        join(method_plots_dir, f"{base}-{metric}-boxplot.pdf"),
                        bbox_inches="tight",
                    )
                    plt.close(grid.figure)

            if "pointplot" in plot_types:
                for metric in plot_metrics:
                    if not _has_metric(metric):
                        continue
                    grid = pointplot_from_data_frame(
                        table,
                        level_col,
                        [metric],
                        group_by="Dataset",
                        color_by="Method",
                        color_palette=color_palette,
                        title_prefix=plot_label,
                        x_order=level_order,
                        show=False,
                    )
                    for y_var, facet in grid.items():
                        facet.savefig(
                            join(
                                method_plots_dir,
                                f"{base}-{y_var}-pointplot.pdf",
                            ),
                            bbox_inches="tight",
                        )
                        plt.close(facet.fig)

            if "heatmap" in plot_types:
                # config "level" means rank (or novel level); keep it in order
                heat_cols = [level_col if c == "level" else c for c in heatmap_cols]
                heat_df = table.copy()
                heat_df[level_col] = pd.Categorical(
                    heat_df[level_col], categories=level_order, ordered=True,
                )
                for metric in plot_metrics:
                    if not _has_metric(metric) or not all(
                        c in heat_df.columns for c in heatmap_rows + heat_cols
                    ):
                        continue
                    ax = heatmap_from_data_frame(
                        heat_df,
                        metric=metric,
                        rows=heatmap_rows,
                        cols=heat_cols,
                        title=f"{plot_label}: {metric}",
                        show=False,
                    )
                    ax.figure.savefig(
                        join(method_plots_dir, f"{base}-{metric}-heatmap.pdf"),
                        bbox_inches="tight",
                    )
                    plt.close(ax.figure)

            if "best_run_stacked_bar" in plot_types:
                _plot_best_runs(
                    table,
                    ratio_df,
                    level_col,
                    [metric for metric in plot_metrics if _has_metric(metric)],
                    best_run_rank,
                    plot_label,
                    base,
                    method_plots_dir,
                    join(summaries_dir, "best_runs", eval_method, f"{base}.csv"),
                    stacked_rank_labels,
                    stacked_rank_axis_label,
                )

    analyze_log_results(cfg)
    Path(join(plots_dir, ".done")).touch()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to config_04_tax_credit.yaml")
    parser.add_argument(
        "--phase",
        choices=["datasets", "manifest", "evaluate", "plot", "post-assign"],
        default="datasets",
    )
    args = parser.parse_args()
    cfg = load_config(args.config)

    if args.phase == "datasets":
        prepare_datasets(cfg)
    elif args.phase == "manifest":
        prepare_manifest(cfg)
    elif args.phase == "evaluate":
        evaluate(cfg)
    elif args.phase == "plot":
        plot_results(cfg)
    elif args.phase == "post-assign":
        evaluate(cfg)
        plot_results(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
