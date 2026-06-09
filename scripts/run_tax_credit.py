#!/usr/bin/env python
"""Tourmaline tax-credit orchestrator: simulation, manifest, evaluation, plotting."""

from __future__ import annotations

import argparse
import hashlib
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
        elif m != "mock-community":
            raise ValueError(f"Unknown evaluation_method: {m}")
    return list(dict.fromkeys(methods))


def analysis_data_subdir(evaluation_method: str) -> str:
    if evaluation_method in ("cross-validated", "cross-validated-taxa"):
        return "cross-validated"
    if evaluation_method == "cross-validated-trad":
        return "cross-validated-trad"
    if evaluation_method == "novel-taxa":
        return "novel-taxa-simulations"
    raise ValueError(evaluation_method)


def param_id(cfg: dict, confidence: float) -> str:
    method = cfg.get("classify_method", "naive-bayes")
    fit_params = (cfg.get("fit_params") or "").strip()
    if method == "naive-bayes" and fit_params:
        digest = hashlib.md5(fit_params.encode()).hexdigest()[:8]
        return f"nb-{digest}-conf{confidence}"
    return f"{method}-conf{confidence}"


def fit_param_id(cfg: dict) -> str:
    """Shared classifier directory name (confidence-independent)."""
    method = cfg.get("classify_method", "naive-bayes")
    fit_params = (cfg.get("fit_params") or "").strip()
    if method == "naive-bayes" and fit_params:
        digest = hashlib.md5(fit_params.encode()).hexdigest()[:8]
        return f"nb-{digest}"
    return method.replace("-", "_")


def _append_fold_assignment_rows(
    rows: list,
    *,
    tmp_root: str,
    subdir: str,
    dataset_id: str,
    reference_id: str,
    query: str,
    ref_seqs: str,
    ref_taxa: str,
    eval_method: str,
    cfg: dict,
    confidences: list,
    job_id_prefix: str,
) -> None:
    """Add manifest rows for one simulated fold (or mock-community combo)."""
    method = cfg.get("classify_method", "naive-bayes")
    multi_conf_nb = method == "naive-bayes" and len(confidences) > 1

    if multi_conf_nb:
        classifier_dir = join(
            tmp_root, subdir, dataset_id, reference_id, method, fit_param_id(cfg)
        )
        fit_job_id = f"{job_id_prefix}-nb-fit"
        rows.append({
            "job_id": fit_job_id,
            "evaluation_method": eval_method,
            "dataset_id": dataset_id,
            "reference_id": reference_id,
            "query_qza": query,
            "ref_seqs": ref_seqs,
            "ref_taxa": ref_taxa,
            "output_dir": classifier_dir,
            "confidence": confidences[0],
            "skip_fit": False,
            "classifier_qza": "",
            "trad_fit": False,
            "fit_only": True,
            "fit_job_id": "",
        })
        classifier_qza = join(classifier_dir, "classifier.qza")
        for conf in confidences:
            p = param_id(cfg, conf)
            assign_dir = join(
                tmp_root, subdir, dataset_id, reference_id, method, p
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
                "classifier_qza": classifier_qza,
                "trad_fit": False,
                "fit_only": False,
                "fit_job_id": fit_job_id,
            })
        return

    for conf in confidences:
        p = param_id(cfg, conf)
        assign_dir = join(tmp_root, subdir, dataset_id, reference_id, method, p)
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
            "classifier_qza": "",
            "trad_fit": False,
            "fit_only": False,
            "fit_job_id": "",
        })


def prepare_datasets(cfg: dict) -> None:
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.framework_functions import generate_simulated_datasets

    out = run_output_dir(cfg)
    os.makedirs(out, exist_ok=True)
    ddir = data_dir(cfg)
    os.makedirs(ddir, exist_ok=True)

    df = reference_dataframe(cfg)
    eval_methods = cfg.get("evaluation_methods", [])
    sim_methods = simulation_methods_for_eval(eval_methods)

    if not sim_methods and "mock-community" not in eval_methods:
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
    from tax_credit.framework_functions import recall_simulated_taxa_dirs
    from tax_credit.paths import QUERY_TAX_ASSIGNMENTS_TXT

    out = run_output_dir(cfg)
    manifest_fp = join(out, "assignment_manifest.tsv")
    rows = []
    ddir = data_dir(cfg)
    tmp_root = join(out, cfg.get("results_tmp_subdir", "assignment-tmp"))
    db_ids = list(reference_dataframe(cfg).index)
    confidences = cfg.get("confidence_values") or [cfg.get("skl_confidence", 0.7)]
    method = cfg.get("classify_method", "naive-bayes")
    pid = param_id(cfg, confidences[0])

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
                        tmp_root=tmp_root,
                        subdir="mock-community",
                        dataset_id=mock_id,
                        reference_id=ref_id,
                        query=query,
                        ref_seqs=join(ddir, "ref_dbs", ref_id, "ref_seqs.qza"),
                        ref_taxa=join(ddir, "ref_dbs", ref_id, "ref_taxa.qza"),
                        eval_method=eval_method,
                        cfg=cfg,
                        confidences=confidences,
                        job_id_prefix=f"mock-{mock_id}-{ref_id}",
                    )
            continue

        subdir = analysis_data_subdir(eval_method)
        sim_dir = join(ddir, subdir)
        if eval_method in ("cross-validated", "cross-validated-taxa"):
            combos, ref_dbs = recall_simulated_taxa_dirs(
                sim_dir, db_ids, cfg["iterations"],
                ref_seqs="ref_seqs.qza", ref_taxa="ref_taxa.qza",
                max_level=6, min_level=cfg.get("cv_recall_min_level", 5),
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
            for db_id in db_ids:
                ref_seqs, ref_taxa = trad_cv_shared_reference_qzas(ddir, db_id)
                classifier_dir = join(
                    tmp_root, "trad-fit", db_id, method, pid
                )
                classifier_qza = join(classifier_dir, "classifier.qza")
                rows.append({
                    "job_id": f"trad-fit-{db_id}-{pid}",
                    "evaluation_method": eval_method,
                    "dataset_id": db_id,
                    "reference_id": db_id,
                    "query_qza": "",
                    "ref_seqs": ref_seqs,
                    "ref_taxa": ref_taxa,
                    "output_dir": classifier_dir,
                    "confidence": confidences[0],
                    "skip_fit": False,
                    "classifier_qza": "",
                    "trad_fit": True,
                    "fit_only": False,
                    "fit_job_id": "",
                })
            combos, ref_dbs = recall_simulated_taxa_dirs(
                sim_dir, db_ids, cfg["iterations"],
                ref_seqs="ref_seqs.qza", ref_taxa="ref_taxa.qza",
                max_level=6, min_level=cfg.get("cv_recall_min_level", 5),
                multilevel=False,
            )
            for dataset_id, reference_id in combos:
                fold_dir = join(sim_dir, dataset_id)
                query = join(fold_dir, "query.qza")
                classifier_qza = join(
                    tmp_root, "trad-fit", reference_id, method, pid, "classifier.qza"
                )
                for conf in confidences:
                    p = param_id(cfg, conf)
                    assign_dir = join(
                        tmp_root, subdir, dataset_id, reference_id, method, p
                    )
                    rows.append({
                        "job_id": f"trad-{dataset_id}-{p}",
                        "evaluation_method": eval_method,
                        "dataset_id": dataset_id,
                        "reference_id": reference_id,
                        "query_qza": query,
                        "ref_seqs": ref_dbs[dataset_id][0],
                        "ref_taxa": ref_dbs[dataset_id][1],
                        "output_dir": assign_dir,
                        "confidence": conf,
                        "skip_fit": True,
                        "classifier_qza": classifier_qza,
                        "trad_fit": False,
                        "fit_only": False,
                        "fit_job_id": "",
                    })
            continue
        else:
            raise ValueError(f"Unknown evaluation_method: {eval_method}")

        for dataset_id, reference_id in combos:
            ref_seqs, ref_taxa = ref_dbs[dataset_id]
            query = join(sim_dir, dataset_id, "query.qza")
            _append_fold_assignment_rows(
                rows,
                tmp_root=tmp_root,
                subdir=subdir,
                dataset_id=dataset_id,
                reference_id=reference_id,
                query=query,
                ref_seqs=ref_seqs,
                ref_taxa=ref_taxa,
                eval_method=eval_method,
                cfg=cfg,
                confidences=confidences,
                job_id_prefix=f"{subdir}-{dataset_id}",
            )

    manifest = pd.DataFrame(rows)
    manifest.to_csv(manifest_fp, sep="\t", index=False)
    return manifest_fp


def collect_results(cfg: dict) -> None:
    from tax_credit.framework_functions import move_results_to_repository
    from tax_credit.paths import list_assignment_result_dirs

    out = run_output_dir(cfg)
    tmp_root = join(out, cfg.get("results_tmp_subdir", "assignment-tmp"))
    repo = join(data_dir(cfg), cfg.get("results_repo_subdir", "self-results"))

    for eval_method in cfg.get("evaluation_methods", []):
        if eval_method == "mock-community":
            continue
        sub = analysis_data_subdir(eval_method)
        method_dirs = list_assignment_result_dirs(join(tmp_root, sub))
        dest = join(repo, sub)
        os.makedirs(dest, exist_ok=True)
        if method_dirs:
            move_results_to_repository(method_dirs, dest)


def evaluate(cfg: dict) -> None:
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.framework_functions import novel_taxa_classification_evaluation
    from tax_credit.mock_evaluation import evaluate_results
    from tax_credit.paths import list_assignment_result_dirs

    out = run_output_dir(cfg)
    summaries_dir = join(out, "summaries")
    os.makedirs(summaries_dir, exist_ok=True)
    ddir = data_dir(cfg)
    repo = join(ddir, cfg.get("results_repo_subdir", "self-results"))
    force = cfg.get("force_evaluation", False)
    summary_names = cfg.get("summary_filenames") or {}

    for eval_method in cfg.get("evaluation_methods", []):
        if eval_method == "mock-community":
            mock_root = cfg.get("mock_dir") or join(ddir, "mock-community")
            tmp_root = join(out, cfg.get("results_tmp_subdir", "assignment-tmp"))
            results_dirs = list_assignment_result_dirs(
                join(tmp_root, "mock-community")
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
        computed = join(repo, sub)
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
        }.get(eval_method, f"evaluate_{sub}.csv")
        summary_fp = join(summaries_dir, summary_names.get(eval_method, default_name))

        test_type = "novel-taxa"
        if eval_method in ("cross-validated", "cross-validated-taxa"):
            test_type = "cross-validated"
        elif eval_method == "cross-validated-trad":
            test_type = "cross-validated-trad"

        novel_taxa_classification_evaluation(
            results_dirs,
            expected,
            summary_fp,
            test_type=test_type,
        )


def plot_results(cfg: dict) -> None:
    if not cfg.get("generate_plots", True):
        return
    _ensure_tax_credit(cfg.get("tax_credit_package_dir", "../tax-credit"))
    from tax_credit.novel_evaluation import extract_per_level_accuracy
    from tax_credit.plotting_functions import boxplot_from_data_frame
    import matplotlib.pyplot as plt

    out = run_output_dir(cfg)
    plots_dir = join(out, cfg.get("plots_subdir", "plots"))
    os.makedirs(plots_dir, exist_ok=True)
    summaries_dir = join(out, "summaries")

    for eval_method in cfg.get("evaluation_methods", []):
        if eval_method == "mock-community":
            continue
        sub = analysis_data_subdir(eval_method)
        for fp in glob(join(summaries_dir, f"*{sub}*")) + glob(join(summaries_dir, "*CV*")):
            if not fp.endswith((".csv", ".tsv")):
                continue
            df = pd.read_csv(fp, index_col=0)
            if df.empty:
                continue
            per_level = extract_per_level_accuracy(df)
            for metric in cfg.get("plot_metrics", ["Precision", "Recall", "F-measure"]):
                if metric not in per_level.columns:
                    continue
                ax = boxplot_from_data_frame(
                    per_level, group_by="Dataset", metric=metric, hue="Method"
                )
                base = Path(fp).stem
                ax.figure.savefig(join(plots_dir, f"{base}-{metric}-boxplot.pdf"))
                plt.close(ax.figure)
    Path(join(plots_dir, ".done")).touch()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to config_04_tax_credit.yaml")
    parser.add_argument(
        "--phase",
        choices=["datasets", "manifest", "collect", "evaluate", "plot", "post-assign"],
        default="datasets",
    )
    args = parser.parse_args()
    cfg = load_config(args.config)

    if args.phase == "datasets":
        prepare_datasets(cfg)
    elif args.phase == "manifest":
        prepare_manifest(cfg)
    elif args.phase == "collect":
        collect_results(cfg)
    elif args.phase == "evaluate":
        evaluate(cfg)
    elif args.phase == "plot":
        plot_results(cfg)
    elif args.phase == "post-assign":
        collect_results(cfg)
        evaluate(cfg)
        plot_results(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
