#!/usr/bin/env python
"""Mock-community evaluation for the Tourmaline tax-credit step.

Used by ``run_tax_credit.py``: config validation, input staging, the
evaluation plan, per-job scoring, summaries and plots. The metrics themselves
live in ``tax_credit.mock_community``.

Staged layout under ``{run}-tax-credit/data/``::

    ref_dbs/{database}/ref_seqs.qza, ref_taxa.qza
    mock-community/
        dataset_log.txt                      what staging kept, dropped and warned about
        evaluation_plan.tsv                  dataset x expected set x database x samples
        datasets/{dataset}/rep_seqs.qza, feature_table.tsv, excluded_samples.tsv
        expected/{expected_set}/composition.tsv, asv_taxonomy.tsv,
                                backbone_check-{database}.tsv,
                                unscored_asvs-{dataset}.tsv
    results/mock-community/{dataset}/{database}/{method}/{parameters}/

Per-job scores go to ``summaries/mock-community/per-job/`` and are combined
into ``summaries/mock_community_metrics.tsv`` and
``summaries/mock_community_composition.tsv``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import textwrap
from datetime import datetime
from os.path import abspath, basename, exists, expandvars, getmtime, join

import pandas as pd

EVAL_METHOD = "mock-community"

LEGACY_KEYS = (
    "mock_communities",
    "mock_dir",
    "mock_taxonomy_level_range",
    "mock_per_seq_precision",
    "mock_force_eval",
)
TABLE_EXTENSIONS = (".tsv", ".txt", ".biom", ".qza")
SEQUENCE_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas", ".qza")
DEFAULT_RANKS = "kingdom,phylum,class,order,family,genus,species"
DEFAULT_PLOT_METRICS = [
    "Taxon Accuracy Rate",
    "Taxon Detection Rate",
    "Bray-Curtis",
    "Precision",
    "Recall",
    "F-measure",
]
BACKBONE_MODES = ("warn", "error")
AUTO_METRIC = "auto"
# tax_credit.mock_community.METRIC_COLUMNS; duplicated so the config can be
# validated before the tax-credit package is on sys.path
SELECTABLE_METRICS = (
    "Taxon Accuracy Rate",
    "Taxon Detection Rate",
    "Bray-Curtis",
    "Precision",
    "Recall",
    "F-measure",
    "ASV Precision",
    "ASV Recall",
    "ASV F-measure",
    "match_ratio",
    "underclassification_ratio",
    "overclassification_ratio",
    "misclassification_ratio",
)

METRICS_SUMMARY = "mock_community_metrics.tsv"
COMPOSITION_SUMMARY = "mock_community_composition.tsv"
ID_COLUMNS = ["MockDataset", "Reference", "ExpectedSet", "Method", "Parameters"]


# open dataset_log.txt while stage_mock_inputs runs; None otherwise
_dataset_log = None


def _wrap_title(title: str, width: int = 78) -> str:
    """Wrap a plot title so long lines stay clear of the figure legend."""
    return "\n".join(textwrap.fill(line, width) for line in title.split("\n"))


def _log(message: str) -> None:
    print(f"[mock-community] {message}", file=sys.stderr, flush=True)
    if _dataset_log is not None:
        _dataset_log.write(message + "\n")
        _dataset_log.flush()


def _path(value) -> str:
    return abspath(expandvars(str(value)))


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [part.strip() for part in value.split(",") if part.strip()]
    return [str(part).strip() for part in value]


def _run(cmd: str) -> None:
    print(cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True)


def mock_enabled(cfg: dict) -> bool:
    return EVAL_METHOD in (cfg.get("evaluation_methods") or [])


# --- Config ----------------------------------------------------------------

def mock_settings(cfg: dict) -> dict:
    """Validate the ``mock_community`` config block and fill in defaults.

    Raises ``ValueError`` listing every problem found.
    """
    errors = []
    legacy = [key for key in LEGACY_KEYS if cfg.get(key)]
    if legacy:
        errors.append(
            f"{', '.join(legacy)} are no longer used. Mock communities are now "
            "configured in the mock_community block; see "
            "config_04_tax_credit.yaml and docs/configuration.md.")

    block = cfg.get("mock_community")
    if not isinstance(block, dict):
        raise ValueError(
            "mock-community is in evaluation_methods but the mock_community "
            "config block is missing or empty.\n" + "\n".join(errors))

    db_ids = [db["id"] for db in cfg.get("reference_databases") or []]
    excluded = set(cfg.get("exclude_databases") or [])
    default_ranks = _as_list(cfg.get("taxa_ranks") or DEFAULT_RANKS)

    datasets, seen = [], set()
    for entry in block.get("datasets") or []:
        dataset_id = str(entry.get("id") or "").strip()
        label = f"mock_community.datasets[{dataset_id or '?'}]"
        errors.extend(_check_id(dataset_id, label, seen))
        dataset = {"id": dataset_id, "samples": _as_list(entry.get("samples"))}
        for key, extensions in (("feature_table", TABLE_EXTENSIONS),
                                ("rep_seqs", SEQUENCE_EXTENSIONS)):
            errors.extend(_check_file(entry.get(key), f"{label}.{key}",
                                      extensions, required=True))
            dataset[key] = _path(entry[key]) if entry.get(key) else None
        datasets.append(dataset)
    if not datasets:
        errors.append("mock_community.datasets needs at least one dataset.")

    expected_sets, seen, db_to_set = [], set(), {}
    for entry in block.get("expected_sets") or []:
        set_id = str(entry.get("id") or "").strip()
        label = f"mock_community.expected_sets[{set_id or '?'}]"
        errors.extend(_check_id(set_id, label, seen))
        if not entry.get("composition") and not entry.get("asv_taxonomy"):
            errors.append(f"{label}: set composition, asv_taxonomy, or both.")
        for key in ("composition", "asv_taxonomy"):
            errors.extend(_check_file(entry.get(key), f"{label}.{key}",
                                      (".tsv", ".txt"), required=False))
        databases = []
        for db_id in _as_list(entry.get("databases")):
            if db_id not in db_ids:
                errors.append(
                    f"{label}: database {db_id!r} is not in reference_databases "
                    f"({db_ids}).")
            elif db_id in db_to_set:
                errors.append(
                    f"{label}: database {db_id!r} is already used by expected "
                    f"set {db_to_set[db_id]!r}. Each database follows one "
                    "taxonomic backbone, so it can belong to only one set.")
            elif db_id in excluded:
                _log(f"{label}: skipping {db_id} (in exclude_databases)")
            else:
                db_to_set[db_id] = set_id
                databases.append(db_id)
        if not _as_list(entry.get("databases")):
            errors.append(f"{label}: databases must list at least one database.")
        ranks = _as_list(entry.get("ranks")) or default_ranks
        expected_sets.append({
            "id": set_id,
            "databases": databases,
            "composition": _path(entry["composition"]) if entry.get("composition") else None,
            "asv_taxonomy": _path(entry["asv_taxonomy"]) if entry.get("asv_taxonomy") else None,
            "ranks": ranks,
        })
    if not expected_sets:
        errors.append("mock_community.expected_sets needs at least one set.")
    elif not db_to_set and not errors:
        errors.append("No reference database is assigned to an expected set.")

    eval_ranks = _as_list(block.get("eval_ranks"))
    if not eval_ranks:
        eval_ranks = [r for r in ("family", "genus", "species") if r in default_ranks]
        eval_ranks = eval_ranks or default_ranks[-3:]
    for expected_set in expected_sets:
        missing = [r for r in eval_ranks if r not in expected_set["ranks"]]
        if missing:
            errors.append(
                f"mock_community.expected_sets[{expected_set['id']}]: eval_ranks "
                f"{missing} are not in its ranks {expected_set['ranks']}.")

    min_abundance = block.get("min_relative_abundance")
    min_abundance = 0.0 if min_abundance is None else min_abundance
    try:
        min_abundance = float(min_abundance)
        if not 0 <= min_abundance < 1:
            raise ValueError
    except (TypeError, ValueError):
        errors.append("mock_community.min_relative_abundance must be in [0, 1).")

    backbone_check = str(block.get("backbone_check") or "warn").strip().lower()
    if backbone_check not in BACKBONE_MODES:
        errors.append(f"mock_community.backbone_check must be one of {BACKBONE_MODES}.")

    best_run_metric = str(block.get("best_run_metric") or AUTO_METRIC).strip()
    if best_run_metric != AUTO_METRIC and best_run_metric not in SELECTABLE_METRICS:
        errors.append(
            f"mock_community.best_run_metric must be {AUTO_METRIC!r} or one of: "
            f"{', '.join(SELECTABLE_METRICS)}.")

    for db in cfg.get("reference_databases") or []:
        if db.get("pretrained_classifier"):
            errors.extend(_check_file(
                db["pretrained_classifier"],
                f"reference_databases[{db['id']}].pretrained_classifier",
                (".qza",), required=True))

    if errors:
        raise ValueError("Invalid mock-community config:\n  - " + "\n  - ".join(errors))

    return {
        "datasets": datasets,
        "expected_sets": expected_sets,
        "database_to_set": db_to_set,
        "eval_ranks": eval_ranks,
        "min_relative_abundance": min_abundance,
        "legacy_unresolved_taxa": bool(block.get("legacy_unresolved_taxa")),
        "backbone_check": backbone_check,
        "best_run_metric": best_run_metric,
        "plot_metrics": _as_list(block.get("plot_metrics")) or DEFAULT_PLOT_METRICS,
        "composition_top_n": int(block.get("composition_top_n") or 12),
    }


def _check_id(value: str, label: str, seen: set) -> list[str]:
    if not value:
        return [f"{label}: id is required."]
    if any(ch in value for ch in "/\\ \t"):
        return [f"{label}: id {value!r} cannot contain slashes or whitespace."]
    if value in seen:
        return [f"{label}: duplicate id {value!r}."]
    seen.add(value)
    return []


def _check_file(value, label: str, extensions, required: bool) -> list[str]:
    if not value:
        return [f"{label} is required."] if required else []
    fp = _path(value)
    if not exists(fp):
        return [f"{label}: file not found: {fp}"]
    if not fp.lower().endswith(tuple(extensions)):
        return [f"{label}: expected one of {', '.join(extensions)}, got {basename(fp)}"]
    return []


# --- Paths -----------------------------------------------------------------

def mock_data_dir(data_root: str) -> str:
    return join(data_root, EVAL_METHOD)


def dataset_dir(data_root: str, dataset_id: str) -> str:
    return join(mock_data_dir(data_root), "datasets", dataset_id)


def expected_dir(data_root: str, set_id: str) -> str:
    return join(mock_data_dir(data_root), "expected", set_id)


def dataset_log_path(data_root: str) -> str:
    return join(mock_data_dir(data_root), "dataset_log.txt")


def plan_path(data_root: str) -> str:
    return join(mock_data_dir(data_root), "evaluation_plan.tsv")


def reference_qzas(data_root: str, db_id: str) -> tuple[str, str]:
    db_dir = join(data_root, "ref_dbs", db_id)
    return join(db_dir, "ref_seqs.qza"), join(db_dir, "ref_taxa.qza")


def per_job_paths(run_output: str, job_id: str) -> tuple[str, str]:
    job_dir = join(run_output, "summaries", EVAL_METHOD, "per-job")
    return (join(job_dir, f"{job_id}.metrics.tsv"),
            join(job_dir, f"{job_id}.composition.tsv"))


# --- Staging ---------------------------------------------------------------

def _link(src: str, dest: str) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.lexists(dest):
        if os.path.islink(dest) and os.readlink(dest) == src:
            return
        os.remove(dest)
    os.symlink(src, dest)


def _up_to_date(dest: str, src: str) -> bool:
    return exists(dest) and not os.path.islink(dest) and getmtime(dest) >= getmtime(src)


def _stage_sequences_qza(src: str, dest: str) -> None:
    if src.lower().endswith(".qza"):
        _link(src, dest)
    elif not _up_to_date(dest, src):
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if os.path.lexists(dest):
            os.remove(dest)
        _run(f"qiime tools import --type 'FeatureData[Sequence]' "
             f"--input-path '{src}' --output-path '{dest}'")


def _stage_taxonomy_qza(src: str, dest: str) -> None:
    if src.lower().endswith(".qza"):
        _link(src, dest)
        return
    if _up_to_date(dest, src):
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.lexists(dest):
        os.remove(dest)
    with open(src, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().startswith("Feature ID")
    input_format = "TSVTaxonomyFormat" if header else "HeaderlessTSVTaxonomyFormat"
    _run(f"qiime tools import --type 'FeatureData[Taxonomy]' "
         f"--input-format {input_format} --input-path '{src}' --output-path '{dest}'")


def _load_feature_table(src: str, work_dir: str) -> pd.DataFrame:
    from tax_credit import mock_community as mc

    if not src.lower().endswith(".qza"):
        return mc.read_feature_table(src)
    export_dir = join(work_dir, "feature-table-export")
    shutil.rmtree(export_dir, ignore_errors=True)
    _run(f"qiime tools export --input-path '{src}' --output-path '{export_dir}'")
    table = mc.read_feature_table(join(export_dir, "feature-table.biom"))
    shutil.rmtree(export_dir, ignore_errors=True)
    return table


def _fasta_ids(fp: str) -> set[str]:
    with open(fp, encoding="utf-8", errors="replace") as fh:
        return {line[1:].split()[0] for line in fh if line.startswith(">") and len(line) > 2}


def stage_mock_inputs(cfg: dict, data_root: str, reference_taxonomy_text: dict) -> None:
    """Stage mock inputs, run the backbone check and write the evaluation plan.

    *reference_taxonomy_text* maps database id to a plain-text taxonomy file.
    Every message is also written to ``dataset_log.txt``, including the error
    that stops staging, if any.
    """
    global _dataset_log
    os.makedirs(mock_data_dir(data_root), exist_ok=True)
    with open(dataset_log_path(data_root), "w", encoding="utf-8") as log:
        log.write(f"# mock-community staging, {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        _dataset_log = log
        try:
            _stage_mock_inputs(cfg, data_root, reference_taxonomy_text)
        except Exception as exc:
            log.write(f"ERROR: {exc}\n")
            raise
        finally:
            _dataset_log = None
    print(f"[mock-community] staging log: {dataset_log_path(data_root)}", file=sys.stderr)


def _stage_mock_inputs(cfg: dict, data_root: str, reference_taxonomy_text: dict) -> None:
    from tax_credit import mock_community as mc

    settings = mock_settings(cfg)
    db_entries = {db["id"]: db for db in cfg["reference_databases"]}

    for db_id in settings["database_to_set"]:
        seqs_qza, taxa_qza = reference_qzas(data_root, db_id)
        _stage_sequences_qza(_path(db_entries[db_id]["refseqs_file"]), seqs_qza)
        _stage_taxonomy_qza(_path(db_entries[db_id]["taxa_file"]), taxa_qza)

    counts_by_dataset = {}
    for dataset in settings["datasets"]:
        ddir = dataset_dir(data_root, dataset["id"])
        os.makedirs(ddir, exist_ok=True)
        _stage_sequences_qza(dataset["rep_seqs"], join(ddir, "rep_seqs.qza"))
        counts = _load_feature_table(dataset["feature_table"], ddir)
        counts.to_csv(join(ddir, "feature_table.tsv"), sep="\t")
        counts_by_dataset[dataset["id"]] = counts
        if not dataset["rep_seqs"].lower().endswith(".qza"):
            missing = counts.index.difference(list(_fasta_ids(dataset["rep_seqs"])))
            if len(missing):
                _log(f"{dataset['id']}: {len(missing)} feature(s) in the feature "
                     f"table have no sequence in rep_seqs and will count as "
                     f"Unassigned: {', '.join(missing)}")
        _log(f"{dataset['id']}: {counts.shape[0]} features, {counts.shape[1]} samples")

    excluded_by_dataset = {dataset["id"]: [] for dataset in settings["datasets"]}
    plan_rows, backbone_problems = [], []
    for expected_set in settings["expected_sets"]:
        if not expected_set["databases"]:
            continue
        edir = expected_dir(data_root, expected_set["id"])
        os.makedirs(edir, exist_ok=True)
        composition = asv_taxonomy = None
        expected_lineages = set()
        if expected_set["composition"]:
            composition = mc.read_composition(expected_set["composition"])
            composition.to_csv(join(edir, "composition.tsv"), sep="\t")
            expected_lineages.update(composition.index)
        if expected_set["asv_taxonomy"]:
            asv_taxonomy = mc.read_taxonomy(expected_set["asv_taxonomy"])
            asv_taxonomy.rename_axis("Feature ID").to_frame().to_csv(
                join(edir, "asv_taxonomy.tsv"), sep="\t")
            expected_lineages.update(asv_taxonomy[asv_taxonomy != ""])

        for db_id in expected_set["databases"]:
            backbone_problems.extend(_backbone_check(
                expected_set, db_id, expected_lineages,
                reference_taxonomy_text[db_id], settings["eval_ranks"], edir))

        for dataset in settings["datasets"]:
            counts = counts_by_dataset[dataset["id"]]
            label = f"{dataset['id']} / {expected_set['id']}"
            try:
                samples, excluded = mc.select_mock_samples(
                    counts, composition, asv_taxonomy, dataset["samples"])
            except ValueError as exc:
                raise ValueError(
                    f"dataset {dataset['id']}, expected set {expected_set['id']}: {exc}"
                ) from exc
            excluded_by_dataset[dataset["id"]].extend(
                {**row, "expected_set": expected_set["id"]} for row in excluded)
            for reason, rows in pd.DataFrame(excluded, columns=["sample_id", "reason"]).groupby("reason"):
                _log(f"{label}: {len(rows)} sample(s) excluded ({reason}): "
                     f"{', '.join(rows['sample_id'])}")

            unscored_fp = join(edir, f"unscored_asvs-{dataset['id']}.tsv")
            if not samples:
                _log(f"{label}: no mock samples in common; this pair is not evaluated")
                if exists(unscored_fp):
                    os.remove(unscored_fp)
                continue
            _log(f"{label}: {len(samples)} mock sample(s) evaluated: {', '.join(samples)}")
            if asv_taxonomy is not None:
                unscored = mc.unscored_asvs(counts, asv_taxonomy, samples)
                unscored.to_csv(unscored_fp, sep="\t", index=False)
                if len(unscored):
                    in_mock = int((unscored["reads_in_mock_samples"] > 0).sum())
                    _log(f"{label}: {len(unscored)} ASV(s) have no known taxonomy "
                         f"({in_mock} with reads in the mock samples) and are left out "
                         f"of precision / recall / F-measure; see {unscored_fp}")
            elif exists(unscored_fp):
                os.remove(unscored_fp)
            for db_id in expected_set["databases"]:
                plan_rows.append({
                    "dataset_id": dataset["id"],
                    "expected_set": expected_set["id"],
                    "database": db_id,
                    "samples": ",".join(samples),
                })

    for dataset_id, rows in excluded_by_dataset.items():
        excluded_fp = join(dataset_dir(data_root, dataset_id), "excluded_samples.tsv")
        pd.DataFrame(rows, columns=["sample_id", "expected_set", "reason"]).to_csv(
            excluded_fp, sep="\t", index=False)
        _log(f"{dataset_id}: excluded samples written to {excluded_fp}")

    if not plan_rows:
        raise ValueError(
            "No mock samples to evaluate: no dataset shares samples with any "
            "expected set. Check sample ids in the feature table and composition.")
    pd.DataFrame(plan_rows).to_csv(plan_path(data_root), sep="\t", index=False)
    _log(f"evaluation plan: {plan_path(data_root)}")

    if backbone_problems and settings["backbone_check"] == "error":
        raise ValueError(
            "Backbone check failed (backbone_check: error):\n  - "
            + "\n  - ".join(backbone_problems))


def _backbone_check(expected_set, db_id, expected_lineages, reference_taxonomy_fp,
                    eval_ranks, out_dir) -> list[str]:
    from tax_credit import mock_community as mc

    reference = mc.read_reference_taxonomy(reference_taxonomy_fp)
    report, messages = mc.check_backbone(
        expected_lineages, reference, expected_set["ranks"], eval_ranks)
    report_fp = join(out_dir, f"backbone_check-{db_id}.tsv")
    report.to_csv(report_fp, sep="\t", index=False)

    label = f"backbone check {expected_set['id']} vs {db_id}"
    problems = [f"{label}: {message}" for message in messages]
    for message in messages:
        _log(f"{label}: WARNING {message}")
    counts = report.groupby(["rank", "status"]).size().unstack(fill_value=0)
    for rank in eval_ranks:
        if rank not in counts.index:
            continue
        row = counts.loc[rank]
        found = int(row.get("found", 0))
        different = int(row.get("different_lineage", 0))
        absent = int(row.get("not_in_database", 0))
        _log(f"{label} at {rank}: {found} found, {different} under a different "
             f"lineage, {absent} not in database")
        if different or absent:
            problems.append(
                f"{label} at {rank}: {different} different_lineage, "
                f"{absent} not_in_database (see {report_fp})")
    return problems


# --- Manifest and scoring --------------------------------------------------

def read_plan(data_root: str) -> pd.DataFrame:
    fp = plan_path(data_root)
    if not exists(fp):
        raise FileNotFoundError(
            f"{fp} not found; run the datasets phase first.")
    return pd.read_csv(fp, sep="\t", dtype=str)


def assignment_jobs(cfg: dict, data_root: str) -> list[dict]:
    """One assignment input per dataset and database in the evaluation plan."""
    plan = read_plan(data_root)
    db_entries = {db["id"]: db for db in cfg["reference_databases"]}
    jobs = []
    for (dataset_id, db_id), _ in plan.groupby(["dataset_id", "database"], sort=False):
        ref_seqs, ref_taxa = reference_qzas(data_root, db_id)
        pretrained = db_entries[db_id].get("pretrained_classifier")
        jobs.append({
            "dataset_id": dataset_id,
            "reference_id": db_id,
            "query": join(dataset_dir(data_root, dataset_id), "rep_seqs.qza"),
            "ref_seqs": ref_seqs,
            "ref_taxa": ref_taxa,
            "pretrained_classifier": _path(pretrained) if pretrained else None,
        })
    return jobs


def _is_true(value) -> bool:
    return str(value).strip().lower() == "true"


def scoring_rows(manifest: pd.DataFrame) -> pd.DataFrame:
    """Manifest rows that produce taxonomy assignments for mock communities."""
    if manifest.empty or "evaluation_method" not in manifest:
        return manifest.iloc[0:0]
    mask = manifest["evaluation_method"] == EVAL_METHOD
    mask &= ~manifest["fit_only"].map(_is_true) & ~manifest["trad_fit"].map(_is_true)
    return manifest[mask]


def evaluate_job(cfg: dict, row, data_root: str, run_output: str) -> tuple[str, str]:
    """Score one assignment job and write its per-job metric tables."""
    from tax_credit import mock_community as mc

    settings = mock_settings(cfg)
    plan = read_plan(data_root)
    entry = plan[(plan["dataset_id"] == row["dataset_id"])
                 & (plan["database"] == row["reference_id"])]
    if entry.empty:
        raise ValueError(
            f"{row['job_id']}: {row['dataset_id']} / {row['reference_id']} is not "
            "in the evaluation plan")
    entry = entry.iloc[0]
    expected_set = next(s for s in settings["expected_sets"]
                        if s["id"] == entry["expected_set"])

    counts = mc.read_feature_table(
        join(dataset_dir(data_root, row["dataset_id"]), "feature_table.tsv"))
    assignments = mc.read_taxonomy(join(row["output_dir"], "query_tax_assignments.txt"))
    unassigned = counts.index.difference(assignments.index)
    if len(unassigned):
        _log(f"{row['job_id']}: {len(unassigned)} feature(s) have no assignment "
             "and count as Unassigned")
    edir = expected_dir(data_root, expected_set["id"])
    composition_fp = join(edir, "composition.tsv")
    asv_taxonomy_fp = join(edir, "asv_taxonomy.tsv")

    metrics, composition = mc.evaluate_mock_samples(
        counts,
        assignments,
        ranks=expected_set["ranks"],
        eval_ranks=settings["eval_ranks"],
        samples=entry["samples"].split(","),
        composition=mc.read_composition(composition_fp) if exists(composition_fp) else None,
        asv_taxonomy=mc.read_taxonomy(asv_taxonomy_fp) if exists(asv_taxonomy_fp) else None,
        min_relative_abundance=settings["min_relative_abundance"],
        legacy_unresolved_taxa=settings["legacy_unresolved_taxa"],
    )
    ids = {
        "MockDataset": row["dataset_id"],
        "Reference": row["reference_id"],
        "ExpectedSet": expected_set["id"],
        "Method": row["classify_method"],
        "Parameters": basename(str(row["output_dir"]).rstrip("/")),
    }
    metrics_fp, composition_fp = per_job_paths(run_output, row["job_id"])
    os.makedirs(os.path.dirname(metrics_fp), exist_ok=True)
    for table, fp in ((metrics, metrics_fp), (composition, composition_fp)):
        for position, (column, value) in enumerate(ids.items()):
            table.insert(position, column, value)
        table.to_csv(fp, sep="\t", index=False)
    return metrics_fp, composition_fp


def write_summaries(cfg: dict, manifest_fp: str, data_root: str, run_output: str) -> None:
    """Combine per-job scores (computing any that are missing) into summaries."""
    manifest = pd.read_csv(manifest_fp, sep="\t", dtype=str, keep_default_na=False)
    rows = scoring_rows(manifest)
    metric_tables, composition_tables = [], []
    for _, row in rows.iterrows():
        metrics_fp, composition_fp = per_job_paths(run_output, row["job_id"])
        if not (exists(metrics_fp) and exists(composition_fp)):
            evaluate_job(cfg, row, data_root, run_output)
        metric_tables.append(pd.read_csv(metrics_fp, sep="\t"))
        composition_tables.append(pd.read_csv(composition_fp, sep="\t"))
    summaries_dir = join(run_output, "summaries")
    os.makedirs(summaries_dir, exist_ok=True)
    if not metric_tables:
        _log("no mock-community assignment jobs in the manifest")
    metrics = pd.concat(metric_tables, ignore_index=True) if metric_tables else pd.DataFrame(columns=ID_COLUMNS)
    composition = (pd.concat(composition_tables, ignore_index=True)
                   if composition_tables else pd.DataFrame(columns=ID_COLUMNS))
    metrics.to_csv(join(summaries_dir, METRICS_SUMMARY), sep="\t", index=False)
    composition.to_csv(join(summaries_dir, COMPOSITION_SUMMARY), sep="\t", index=False)


# --- Plots -----------------------------------------------------------------

def _best_runs_per_method(table, settings, best_rank):
    """Best parameter set for each dataset, database and method.

    The metric is ``mock_community.best_run_metric``; with ``auto`` it is
    F-measure where the expected set has an ASV taxonomy, and Bray-Curtis
    otherwise, so it is chosen per database. Returns the rows of
    ``select_best_runs`` for every group, with a ``metric`` column.
    """
    from tax_credit.novel_evaluation import select_best_runs

    configured = settings["best_run_metric"]
    at_rank = table[table["rank"] == best_rank]
    picks = []
    for (dataset_id, reference), scores in at_rank.groupby(["MockDataset", "Reference"]):
        if configured == AUTO_METRIC:
            metric = "F-measure" if scores["F-measure"].notna().any() else "Bray-Curtis"
        else:
            metric = configured
            if metric not in scores.columns or scores[metric].isna().all():
                _log(f"{dataset_id} / {reference}: no {metric} scores at {best_rank} "
                     "rank; no best run selected (set mock_community.best_run_metric "
                     "to a metric this expected set produces)")
                continue
        group = select_best_runs(scores, [metric],
                                 group_cols=["MockDataset", "Reference", "Method"],
                                 run_cols=["Parameters"])
        if not group.empty:
            picks.append(group)
    if not picks:
        return pd.DataFrame()
    return pd.concat(picks, ignore_index=True)


def plot_mock(cfg: dict, run_output: str, plots_dir: str, plot_types: set,
              palette_override=None) -> None:
    """Metric plots and best-run composition bars for mock communities."""
    import matplotlib.pyplot as plt
    from tax_credit.novel_evaluation import select_best_runs
    from tax_credit.plot_theme import method_palette, metric_label
    from tax_credit.plotting_functions import (
        composition_barplot_from_data_frame,
        faceted_boxplot_from_data_frame,
        heatmap_from_data_frame,
        pointplot_from_data_frame,
        stacked_classification_barplot_from_data_frame,
    )

    settings = mock_settings(cfg)
    summaries_dir = join(run_output, "summaries")
    metrics_fp = join(summaries_dir, METRICS_SUMMARY)
    if not exists(metrics_fp):
        _log(f"no summary at {metrics_fp}; skipping plots")
        return
    table = pd.read_csv(metrics_fp, sep="\t")
    if table.empty:
        return
    out_dir = join(plots_dir, EVAL_METHOD)
    os.makedirs(out_dir, exist_ok=True)

    eval_ranks = settings["eval_ranks"]
    box_ranks = [r for r in _as_list(cfg.get("plot_ranks")) if r in eval_ranks] or eval_ranks
    best_rank = str(cfg.get("best_run_rank") or "").strip().lower()
    best_rank = best_rank if best_rank in eval_ranks else eval_ranks[-1]
    metrics = [m for m in settings["plot_metrics"]
               if m in table.columns and table[m].notna().any()]
    palette = method_palette(table["Method"].unique(), palette_override)
    ratio_cols = ["match_ratio", "underclassification_ratio",
                  "overclassification_ratio", "misclassification_ratio"]

    # best parameter set per dataset, database and method: plotted on their own
    # and used for the composition bars
    picks = _best_runs_per_method(table, settings, best_rank)
    best_dir = join(summaries_dir, "best_runs", EVAL_METHOD)
    best_table = pd.DataFrame()
    if not picks.empty:
        os.makedirs(best_dir, exist_ok=True)
        picks.insert(3, "rank", best_rank)
        picks.to_csv(join(best_dir, "mock_community_best_per_method.csv"), index=False)
        best_table = table.merge(
            picks[["MockDataset", "Reference", "Method", "Parameters"]],
            on=["MockDataset", "Reference", "Method", "Parameters"],
        )

    def _save(fig, filename):
        fig.savefig(join(out_dir, filename))
        plt.close(fig)

    def _selection_note(dataset_picks, label_metric=metric_label):
        """How the best runs were chosen, e.g. 'selected at species by highest
        mean F-measure'.

        With `auto`, databases whose expected set has no ASV taxonomy fall back
        to Bray-Curtis, so when the metric differs each one names its
        databases. Beyond two, the CSV is cited instead of listing them.
        """
        chosen = dataset_picks[["metric", "direction", "Reference"]].drop_duplicates()
        groups = list(chosen.groupby(["metric", "direction"]))
        if len(groups) > 2:
            return (f"selected at {best_rank} by the metric each database supports; "
                    "see mock_community_best_per_method.csv")
        parts = []
        for (metric_name, direction), rows in groups:
            part = f"{direction} mean {label_metric(metric_name)}"
            if len(groups) > 1:
                part += f" ({', '.join(sorted(rows['Reference']))})"
            parts.append(part)
        return f"selected at {best_rank} by {' and '.join(parts)}"

    for dataset_id, df in table.groupby("MockDataset", sort=True):
        label = f"Mock community {dataset_id}"
        prefix = f"{EVAL_METHOD}-{dataset_id}"
        best_df = (best_table[best_table["MockDataset"] == dataset_id]
                   if not best_table.empty else pd.DataFrame())
        note = (_selection_note(picks[picks["MockDataset"] == dataset_id])
                if not picks.empty else "")
        for metric in metrics:
            if "boxplot" in plot_types:
                _save(faceted_boxplot_from_data_frame(
                    df[df["rank"].isin(box_ranks)], x="Reference", metric=metric,
                    hue="Method", col="rank", col_order=box_ranks, palette=palette,
                    title=_wrap_title(f"{label}: {metric_label(metric)}\n(one point per sample and parameter set)"),
                ), f"{prefix}-{metric}-boxplot.pdf")
                if not best_df.empty:
                    _save(faceted_boxplot_from_data_frame(
                        best_df[best_df["rank"].isin(box_ranks)], x="Reference",
                        metric=metric, hue="Method", col="rank", col_order=box_ranks,
                        palette=palette,
                        title=_wrap_title(f"{label}: {metric_label(metric)}, one point per sample\n"
                                          f"best parameters per method ({note})"),
                    ), f"{prefix}-{metric}-best-boxplot.pdf")
            if "pointplot" in plot_types:
                _save(pointplot_from_data_frame(
                    df, x="rank", metric=metric, hue="Method", col="Reference",
                    x_order=eval_ranks, palette=palette, x_label="rank",
                    title=_wrap_title(f"{label}: {metric_label(metric)}"),
                ), f"{prefix}-{metric}-pointplot.pdf")
                if not best_df.empty:
                    _save(pointplot_from_data_frame(
                        best_df, x="rank", metric=metric, hue="Method", col="Reference",
                        x_order=eval_ranks, palette=palette, x_label="rank",
                        title=_wrap_title(f"{label}: {metric_label(metric)}\n"
                                          f"best parameters per method ({note})"),
                    ), f"{prefix}-{metric}-best-pointplot.pdf")
            if "heatmap" in plot_types:
                heat_df = df.copy()
                heat_df["rank"] = pd.Categorical(heat_df["rank"], categories=eval_ranks,
                                                 ordered=True)
                _save(heatmap_from_data_frame(
                    heat_df, metric=metric, rows=["Method", "Parameters"],
                    cols=["Reference", "rank"],
                    title=_wrap_title(f"{label}: mean {metric_label(metric)} across samples"),
                ), f"{prefix}-{metric}-heatmap.pdf")

        if "stacked_bar" in plot_types and df[ratio_cols].notna().any().any():
            ratios = (df.groupby(["Reference", "Method", "Parameters", "level"],
                                 as_index=False)[ratio_cols].mean())
            rank_labels = dict(zip(df["level"], df["rank"]))
            _save(stacked_classification_barplot_from_data_frame(
                ratios, col="Reference", level_col="level",
                level_labels={lvl: name[0].upper() for lvl, name in rank_labels.items()},
                level_axis_label="rank (" + ", ".join(
                    f"{name[0].upper()} = {name}" for name in eval_ranks) + ")",
                title=_wrap_title(f"{label}: read-weighted classification ratios by rank"),
            ), f"{prefix}-classification-ratios-stacked-barplot.pdf")

    if "best_run_stacked_bar" not in plot_types or picks.empty:
        return
    # the summary tables rank every run per metric, alongside the per-method picks
    at_rank = table[table["rank"] == best_rank]
    ranked = select_best_runs(at_rank, metrics, group_cols=["MockDataset", "Reference"])
    if not ranked.empty:
        ranked.insert(2, "rank", best_rank)
        ranked.to_csv(join(best_dir, "mock_community_best_runs.csv"), index=False)
    _plot_best_compositions(settings, picks, summaries_dir, out_dir, best_rank,
                            composition_barplot_from_data_frame, _selection_note)


def _plot_best_compositions(settings, picks, summaries_dir, out_dir, best_rank,
                            barplot, note) -> None:
    """Expected vs the best parameter set of each method, per sample.

    One figure per mock dataset: one row of panels per reference database, one
    panel per mock sample, and one bar per run. Every panel shares the taxon
    colours, so databases can be compared directly.
    """
    import matplotlib.pyplot as plt

    composition_fp = join(summaries_dir, COMPOSITION_SUMMARY)
    if not exists(composition_fp):
        return
    composition = pd.read_csv(composition_fp, sep="\t")
    composition = composition[composition["rank"] == best_rank]

    for dataset_id, dataset_picks in picks.groupby("MockDataset"):
        frames, runs_by_reference = [], {}
        for reference, group in dataset_picks.groupby("Reference"):
            subset = composition[(composition["MockDataset"] == dataset_id)
                                 & (composition["Reference"] == reference)]
            if subset.empty:
                continue
            # expected abundance is the same for every run of a database
            first = group.iloc[0]
            expected = subset[(subset["Method"] == first["Method"])
                              & (subset["Parameters"] == first["Parameters"])]
            frames.append(expected.assign(run="Expected", abundance=expected["expected"]))
            runs = ["Expected"]
            for pick in group.sort_values("Method").itertuples(index=False):
                run = f"{pick.Method}\n{pick.Parameters}"
                observed = subset[(subset["Method"] == pick.Method)
                                  & (subset["Parameters"] == pick.Parameters)]
                frames.append(observed.assign(run=run, abundance=observed["observed"]))
                runs.append(run)
            runs_by_reference[reference] = runs
        if not frames:
            continue
        plot_df = pd.concat(frames, ignore_index=True)
        plot_df = plot_df[plot_df["abundance"] > 0]
        if plot_df.empty:
            continue
        fig = barplot(
            plot_df, runs_by_reference, row_col="Reference",
            top_n=settings["composition_top_n"],
            title=_wrap_title(f"Mock community {dataset_id}: expected vs best parameters "
                              f"per method\n({note(dataset_picks)})"),
        )
        fig.savefig(join(out_dir, f"{EVAL_METHOD}-{dataset_id}-composition-{best_rank}.pdf"))
        plt.close(fig)
