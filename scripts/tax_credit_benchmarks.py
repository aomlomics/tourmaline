#!/usr/bin/env python
"""Summarise the tax-credit step's Snakemake benchmark files.

Every tax-credit rule declares a ``benchmark:`` file under
``<run>-tax-credit/benchmarks/``: one TSV per job, holding wall time, CPU time,
peak memory and mean CPU load for the job's whole process tree. On its own a
directory of hundreds of per-job files says little, and none of them records
what the job *was*. This module joins the per-assignment files to
``assignment_manifest.tsv`` on ``job_id`` and writes two tables to
``summaries/``:

``job_benchmarks.tsv``
    One row per assignment job: manifest identity columns, a ``job_kind``
    (what the job actually ran), and the benchmark measurements.
``benchmark_summary.tsv``
    Totals per evaluation method x classify method x job kind, plus one row per
    pipeline phase (dataset generation, manifest, mock scoring, evaluation,
    plotting), sorted by total wall time, so the dominant cost is the first
    row.

Snakemake writes a benchmark file only for a job that succeeds, and records
only the attempt that succeeded. Jobs that finished before benchmarking was
added, or that ran outside Snakemake, have no file; they are kept with blank
measurements and counted in ``n_measured``.

Pure pandas: no QIIME 2 import, so it can be tested anywhere.
"""

from __future__ import annotations

import os
from os.path import exists, join

import pandas as pd

BENCHMARKS_SUBDIR = "benchmarks"
ASSIGN_SUBDIR = "assign"
MOCK_SCORE_SUBDIR = "mock_score"

# Single-job phases: benchmark file name (without .tsv) -> label in the summary.
PHASE_FILES = {
    "prepare_datasets": "dataset generation",
    "prepare_manifest": "manifest",
    "evaluate": "evaluation",
    "plot": "plotting",
}

# Snakemake benchmark columns kept, in output order. Units: seconds, MB.
# mean_load is CPU percent averaged over the job (100 = one full core).
MEASURE_COLUMNS = [
    "s", "cpu_time", "max_rss", "max_pss", "max_uss", "mean_load",
    "io_in", "io_out",
]
MANIFEST_COLUMNS = [
    "job_id", "evaluation_method", "reference_id", "dataset_id",
    "classify_method",
]
GROUP_COLUMNS = ["evaluation_method", "classify_method", "job_kind"]


def benchmarks_dir(run_output: str) -> str:
    return join(run_output, BENCHMARKS_SUBDIR)


def read_benchmark(fp: str) -> dict:
    """Measurements from one Snakemake benchmark TSV, as floats.

    A file normally holds one row; with ``repeat()`` it holds one per repeat,
    and the last is used. ``NA`` / ``-`` (no data collected, e.g. a job too
    short to sample) become NaN.
    """
    df = pd.read_csv(fp, sep="\t", na_values=["NA", "-"])
    if df.empty:
        return {col: float("nan") for col in MEASURE_COLUMNS}
    row = df.iloc[-1]
    record = {
        col: pd.to_numeric(row.get(col), errors="coerce") for col in MEASURE_COLUMNS
    }
    if resources_unreadable(record):
        # Only the wall time is real; the rest are placeholders, not zeros.
        record.update({col: float("nan") for col in MEASURE_COLUMNS if col != "s"})
    return record


def resources_unreadable(record: dict) -> bool:
    """True when Snakemake could time a job but not read its processes.

    Snakemake reads each process with psutil's ``memory_full_info()`` and skips
    any process that call is refused for -- dropping its CPU time, memory and
    I/O along with it. macOS refuses that call for every process (it needs
    privileges for the USS figure), so on a Mac every job records 0 for all of
    them, plus a meaningless partial ``mean_load``. No real job runs for
    seconds with zero resident memory and zero CPU, so that combination means
    "not measured". Linux reads its own processes, and reports all columns.
    """
    return (
        record.get("s", 0) > 0
        and record.get("max_rss") == 0
        and record.get("cpu_time") == 0
    )


def _is_true(value) -> bool:
    return str(value).strip().lower() == "true"


def _blank(value) -> bool:
    return str(value).strip().lower() in ("", "nan", "na", "none")


def job_kind(row) -> str:
    """What one manifest row actually ran, from its dispatch columns.

    Mirrors the dispatch in tax_credit_step.Snakefile and
    run_taxonomy_assignment_fold.py: ``blca_stage`` takes precedence, then the
    fit-only flags, then whether the job builds its own classifier, index or
    database before classifying.
    """
    stage = row.get("blca_stage", "")
    if not _blank(stage):
        # align = bowtie2 index + alignment; blca = muscle + bootstraps;
        # reformat = confidence threshold only
        return f"blca-{str(stage).strip()}"
    if _is_true(row.get("fit_only")) or _is_true(row.get("trad_fit")):
        return "fit"
    method = str(row.get("classify_method", "")).strip()
    if method == "consensus-vsearch":
        # vsearch has no separate fit; every job just searches the reference
        return "classify"
    if method == "revamp" and not _blank(row.get("revamp_blast_results", "")):
        return "classify"
    if _is_true(row.get("skip_fit")):
        return "classify"
    if method == "bt2-blca":
        # a single-combination fold runs align, BLCA and reformat in one job
        return "blca-all-stages"
    return "fit+classify"


def job_benchmarks(manifest_fp: str, assign_dir: str) -> pd.DataFrame:
    """One row per manifest job, with its benchmark measurements (or NaN)."""
    manifest = pd.read_csv(manifest_fp, sep="\t", dtype=str, keep_default_na=False)
    rows = []
    for _, m in manifest.iterrows():
        record = {col: m.get(col, "") for col in MANIFEST_COLUMNS}
        record["job_kind"] = job_kind(m)
        fp = join(assign_dir, f"{m['job_id']}.tsv")
        record["measured"] = exists(fp)
        record.update(
            read_benchmark(fp) if exists(fp)
            else {col: float("nan") for col in MEASURE_COLUMNS}
        )
        rows.append(record)
    columns = MANIFEST_COLUMNS + ["job_kind", "measured"] + MEASURE_COLUMNS
    return pd.DataFrame(rows, columns=columns)


def _stat(series: pd.Series, how: str) -> float:
    """*how* ('sum', 'median', 'max', 'mean') over the measured values, else NaN.

    Checked up front: numpy warns "Mean of empty slice" on an all-NaN median.
    """
    values = series.dropna()
    return float(getattr(values, how)()) if len(values) else float("nan")


def _aggregate(df: pd.DataFrame) -> dict:
    """Totals for a group of jobs; unmeasured jobs are skipped."""
    return {
        "n_jobs": len(df),
        "n_measured": int(df["measured"].sum()),
        "total_wall_s": _stat(df["s"], "sum"),
        "total_cpu_s": _stat(df["cpu_time"], "sum"),
        "median_wall_s": _stat(df["s"], "median"),
        "max_wall_s": _stat(df["s"], "max"),
        "peak_rss_mb": _stat(df["max_rss"], "max"),
        "peak_pss_mb": _stat(df["max_pss"], "max"),
        "mean_load_pct": _stat(df["mean_load"], "mean"),
    }


def phase_rows(bench_dir: str) -> list[dict]:
    """Summary rows for the non-assignment rules."""
    rows = []
    for name, label in PHASE_FILES.items():
        fp = join(bench_dir, f"{name}.tsv")
        measured = exists(fp)
        m = read_benchmark(fp) if measured else {c: float("nan") for c in MEASURE_COLUMNS}
        df = pd.DataFrame([{**m, "measured": measured}])
        rows.append({"evaluation_method": "-", "classify_method": "-",
                     "job_kind": label, **_aggregate(df)})
    mock_dir = join(bench_dir, MOCK_SCORE_SUBDIR)
    mock_files = sorted(
        join(mock_dir, f) for f in os.listdir(mock_dir) if f.endswith(".tsv")
    ) if os.path.isdir(mock_dir) else []
    if mock_files:
        df = pd.DataFrame([{**read_benchmark(f), "measured": True} for f in mock_files])
        rows.append({"evaluation_method": "mock-community", "classify_method": "-",
                     "job_kind": "mock scoring", **_aggregate(df)})
    return rows


def benchmark_summary(jobs: pd.DataFrame, bench_dir: str) -> pd.DataFrame:
    """Per-group totals for assignment jobs plus one row per pipeline phase."""
    rows = [
        {**dict(zip(GROUP_COLUMNS, key)), **_aggregate(group)}
        for key, group in jobs.groupby(GROUP_COLUMNS, sort=False)
    ]
    rows += phase_rows(bench_dir)
    summary = pd.DataFrame(rows)
    total = summary["total_wall_s"].sum(min_count=1)
    summary["pct_of_total_wall"] = (
        100 * summary["total_wall_s"] / total if pd.notna(total) and total > 0
        else float("nan")
    )
    return summary.sort_values(
        "total_wall_s", ascending=False, na_position="last"
    ).reset_index(drop=True)


def write_benchmark_summaries(run_output: str, summaries_dir: str) -> tuple[str, str]:
    """Write job_benchmarks.tsv and benchmark_summary.tsv; return their paths."""
    bench_dir = benchmarks_dir(run_output)
    jobs = job_benchmarks(
        join(run_output, "assignment_manifest.tsv"), join(bench_dir, ASSIGN_SUBDIR)
    )
    summary = benchmark_summary(jobs, bench_dir)
    os.makedirs(summaries_dir, exist_ok=True)
    jobs_fp = join(summaries_dir, "job_benchmarks.tsv")
    summary_fp = join(summaries_dir, "benchmark_summary.tsv")
    jobs.to_csv(jobs_fp, sep="\t", index=False, float_format="%.2f")
    summary.to_csv(summary_fp, sep="\t", index=False, float_format="%.2f")

    n_missing = int((~jobs["measured"]).sum())
    print(f"Benchmarks: {len(jobs) - n_missing} of {len(jobs)} assignment jobs measured; "
          f"wrote {jobs_fp} and {summary_fp}")
    if n_missing:
        print(f"  {n_missing} job(s) have no benchmark file (finished before benchmarking "
              "was added, failed, or ran outside Snakemake); their measurements are blank.")
    n_no_resources = int((jobs["measured"] & jobs["s"].notna()
                          & jobs["max_rss"].isna()).sum())
    if n_no_resources:
        print(f"  {n_no_resources} job(s) have wall time only: Snakemake could not read "
              "their processes' CPU and memory (expected on macOS, which refuses it; "
              "Linux reports them). Those columns are blank, not zero.")
    top = summary.dropna(subset=["total_wall_s"]).head(5)
    for _, r in top.iterrows():
        print(f"  {r['pct_of_total_wall']:5.1f}%  {r['total_wall_s'] / 3600:8.2f} h  "
              f"{r['evaluation_method']} / {r['classify_method']} / {r['job_kind']} "
              f"({int(r['n_jobs'])} jobs)")
    return jobs_fp, summary_fp
