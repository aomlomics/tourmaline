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
    return cfg


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


def _param_value_list(cfg: dict, key: str, default: float) -> list[float]:
    val = cfg.get(key, default)
    if isinstance(val, (list, tuple)):
        return [float(v) for v in val]
    return [float(val)]


def consensus_param_combinations(cfg: dict) -> list[dict[str, float]]:
    """Cartesian product of consensus classifier parameters from config."""
    perc = _param_value_list(cfg, "perc_identity", 0.8)
    qc = _param_value_list(cfg, "query_cov", 0.8)
    mc = _param_value_list(cfg, "min_consensus", 0.51)
    return [
        {"perc_identity": pi, "query_cov": q, "min_consensus": m}
        for pi, q, m in itertools.product(perc, qc, mc)
    ]


def _format_sweep_param(value: float) -> str:
    return format(float(value), "g")


def consensus_param_id(combo: dict[str, float]) -> str:
    return (
        f"pi{_format_sweep_param(combo['perc_identity'])}-"
        f"qc{_format_sweep_param(combo['query_cov'])}-"
        f"mc{_format_sweep_param(combo['min_consensus'])}"
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
        "perc_identity": cfg.get("perc_identity", 0.8),
        "query_cov": cfg.get("query_cov", 0.8),
        "min_consensus": cfg.get("min_consensus", 0.51),
    }
    if method == "bt2-blca":
        settings["taxa_ranks"] = cfg.get(
            "taxa_ranks", "kingdom,phylum,class,order,family,genus,species"
        )
    return settings


def confidences_for_method(cfg: dict, method: str) -> list[float]:
    if method == "naive-bayes":
        return cfg.get("confidence_values") or [cfg.get("skl_confidence", 0.7)]
    if method == "bt2-blca":
        return cfg.get("confidence_values") or [cfg.get("confidence_thres", 0.8)]
    return [cfg.get("skl_confidence", 0.7)]


def param_id(
    method: str,
    method_cfg: dict,
    confidence: float,
    consensus_combo: dict[str, float] | None = None,
) -> str:
    if method in ("consensus-blast", "consensus-vsearch"):
        if consensus_combo is None:
            raise ValueError("consensus_combo is required for consensus methods")
        return consensus_param_id(consensus_combo)
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
) -> dict:
    if method in ("consensus-blast", "consensus-vsearch"):
        if consensus_combo is None:
            raise ValueError("consensus_combo is required for consensus methods")
        return {
            "classify_method": method,
            "fit_params": method_cfg.get("fit_params") or "",
            "classify_params": method_cfg.get("classify_params") or "",
            "perc_identity": float(consensus_combo["perc_identity"]),
            "query_cov": float(consensus_combo["query_cov"]),
            "min_consensus": float(consensus_combo["min_consensus"]),
            "taxa_ranks": "",
        }
    params = {
        "classify_method": method,
        "fit_params": method_cfg.get("fit_params") or "",
        "classify_params": method_cfg.get("classify_params") or "",
        "perc_identity": float(method_cfg.get("perc_identity", 0.8)),
        "query_cov": float(method_cfg.get("query_cov", 0.8)),
        "min_consensus": float(method_cfg.get("min_consensus", 0.51)),
        "taxa_ranks": "",
    }
    if method == "bt2-blca":
        params["taxa_ranks"] = method_cfg.get(
            "taxa_ranks", "kingdom,phylum,class,order,family,genus,species"
        )
    return params


def _empty_manifest_artifact_fields() -> dict:
    return {
        "classifier_qza": "",
        "bowtie_index_dir": "",
    }


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
            p = param_id(method, method_cfg, confidences[0], combo)
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
                "confidence": confidences[0],
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
            "confidence": confidences[0],
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
                "confidence": conf,
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
            "confidence": conf,
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
                            "confidence": confidences[0],
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
                                "confidence": conf,
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


def _prepare_plot_data(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame | None, pd.DataFrame]:
    """Return (per-rank summary, raw summary) for plotting."""
    from tax_credit.novel_evaluation import extract_per_level_accuracy

    if _summary_has_per_level_lists(df):
        return extract_per_level_accuracy(df), df
    return None, df


def _metric_plot_df(
    per_level: pd.DataFrame | None,
    raw: pd.DataFrame,
    metric: str,
) -> pd.DataFrame | None:
    """Pick the dataframe that contains *metric* for plotting."""
    if per_level is not None and metric in per_level.columns:
        return per_level
    if metric in raw.columns:
        return raw
    return None


def _plot_type_set(cfg: dict) -> set[str]:
    """Normalize configured plot types (boxplot, pointplot, heatmap, stacked_bar)."""
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

        for sensitivity_rank in ranks:
            sensitivity = summarize_method_parameter_sensitivity(
                log_df,
                rank=sensitivity_rank,
                min_obvs=min_obvs,
                metric="pct_mis",
            )
            if sensitivity.empty:
                continue
            sensitivity.to_csv(
                join(
                    method_summary_dir,
                    f"method_parameter_sensitivity-{sensitivity_rank}.csv",
                ),
            )
            plot_pivot = select_top_sensitivity_taxa(sensitivity, top_n=top_n)
            ax = method_parameter_sensitivity_heatmap_from_data_frame(
                plot_pivot,
                title=(
                    f"{plot_label}: misclassification by taxon and run "
                    f"({sensitivity_rank})"
                ),
                show=False,
            )
            ax.figure.savefig(
                join(
                    method_plots_dir,
                    f"method-parameter-sensitivity-{sensitivity_rank}.pdf",
                ),
                bbox_inches="tight",
            )
            plt.close(ax.figure)


def plot_results(cfg: dict) -> None:
    if not cfg.get("generate_plots", True):
        return
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.plotting_functions import (
        boxplot_from_data_frame,
        heatmap_from_data_frame,
        pointplot_from_data_frame,
        stacked_classification_barplot_from_data_frame,
    )
    from tax_credit.novel_evaluation import extract_per_level_classification_ratios
    from tax_credit.paths import list_assignment_result_dirs
    import matplotlib.pyplot as plt
    import seaborn as sns

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

        if "stacked_bar" in plot_types:
            computed = join(results_dir(cfg), sub)
            assignment_dirs = list_assignment_result_dirs(computed)
            stacked_df = extract_per_level_classification_ratios(assignment_dirs)
            if not stacked_df.empty:
                ax = stacked_classification_barplot_from_data_frame(
                    stacked_df,
                    title=f"{plot_label}: classification ratios by level",
                    show=False,
                )
                ax.figure.savefig(
                    join(
                        method_plots_dir,
                        "classification-ratios-stacked-barplot.pdf",
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
            per_level, raw = _prepare_plot_data(df)

            if "boxplot" in plot_types:
                for metric in plot_metrics:
                    metric_df = _metric_plot_df(per_level, raw, metric)
                    if metric_df is None:
                        continue
                    ax = boxplot_from_data_frame(
                        metric_df,
                        group_by="Dataset",
                        metric=metric,
                        hue="Method",
                        color_palette=color_palette,
                        plotf=sns.boxplot,
                        title=f"{plot_label}: {metric}",
                        show=False,
                    )
                    ax.figure.savefig(
                        join(method_plots_dir, f"{base}-{metric}-boxplot.pdf"),
                        bbox_inches="tight",
                    )
                    plt.close(ax.figure)

            if "pointplot" in plot_types:
                for metric in plot_metrics:
                    metric_df = _metric_plot_df(per_level, raw, metric)
                    if metric_df is None or not {
                        "Dataset", "Method", "level"
                    }.issubset(metric_df.columns):
                        continue
                    grid = pointplot_from_data_frame(
                        metric_df,
                        "level",
                        [metric],
                        group_by="Dataset",
                        color_by="Method",
                        color_palette=color_palette,
                        title_prefix=plot_label,
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
                for metric in plot_metrics:
                    metric_df = _metric_plot_df(per_level, raw, metric)
                    if metric_df is None or not all(
                        c in metric_df.columns
                        for c in heatmap_rows + heatmap_cols
                    ):
                        continue
                    ax = heatmap_from_data_frame(
                        metric_df,
                        metric=metric,
                        rows=heatmap_rows,
                        cols=heatmap_cols,
                        title=f"{plot_label}: {metric}",
                        show=False,
                    )
                    ax.figure.savefig(
                        join(method_plots_dir, f"{base}-{metric}-heatmap.pdf"),
                        bbox_inches="tight",
                    )
                    plt.close(ax.figure)

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
