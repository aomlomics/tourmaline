## Tourmaline tax-credit Snakemake workflow.
## Invoked via `tourmaline.sh --step tax-credit`.
import os
import shutil

import pandas as pd

output_dir = config["output_dir"] + "/"
run_output = output_dir + config["run_name"] + "-tax-credit/"
manifest_fp = run_output + "assignment_manifest.tsv"
datasets_done = run_output + ".datasets.done"
assign_done_dir = run_output + "assignment-done/"
summaries_dir = run_output + "summaries/"
evaluate_done = summaries_dir + ".evaluate.done"
mock_jobs_dir = summaries_dir + "mock-community/per-job/"
plots_done = run_output + config.get("plots_subdir", "plots") + "/.done"
# Per-job wall time / CPU / peak memory, written by Snakemake's benchmark:
# directive, and summarised by tax_credit_benchmark_summary. See
# scripts/tax_credit_benchmarks.py for the layout and columns.
benchmarks_dir = run_output + "benchmarks/"
benchmark_summary_tsv = summaries_dir + "benchmark_summary.tsv"
job_benchmarks_tsv = summaries_dir + "job_benchmarks.tsv"

# bt2-blca has its own cutoffs; refuse older configs that shared perc_identity /
# query_cov with the consensus methods (checked before anything is written).
_classify_methods = config.get("classify_methods") or config.get(
    "classify_method", "naive-bayes"
)
if isinstance(_classify_methods, str):
    _classify_methods = [_classify_methods]
if "bt2-blca" in _classify_methods:
    _missing_blca = [
        key for key in ("blca_perc_identity", "blca_query_cov") if config.get(key) is None
    ]
    if _missing_blca:
        raise ValueError(
            f"bt2-blca is in classify_methods but {', '.join(_missing_blca)} is not set. "
            "bt2-blca no longer reads perc_identity / query_cov (those are for "
            "consensus-blast / consensus-vsearch only). Add blca_perc_identity and "
            "blca_query_cov to the BT2-BLCA OPTIONS section of the config; see "
            "config_04_tax_credit.yaml."
        )

# bt2-blca cost multipliers. Single run-wide values, not sweeps: bowtie2 -k is part
# of the align stage that every BLCA parameter set in a fold shares. Defaults are
# the values that used to be hardcoded. Checked before anything is written.
_BLCA_MULTIPLIER_DEFAULTS = {
    "blca_bowtie2_max_hits": 100,  # bowtie2 -k   (align stage)
    "blca_bootstraps": 100,        # BLCA -n      (blca stage)
    "blca_muscle_maxiters": 16,    # muscle -maxiters (blca stage)
}


def _blca_multiplier(key):
    val = config.get(key)
    if val is None or (isinstance(val, str) and not val.strip()):
        return _BLCA_MULTIPLIER_DEFAULTS[key]
    if isinstance(val, bool) or not isinstance(val, int) or val < 1:
        raise ValueError(
            f"{key} must be a single positive whole number (default "
            f"{_BLCA_MULTIPLIER_DEFAULTS[key]}), got {val!r}. It is a run-wide "
            "setting, not a list to sweep."
        )
    return val


_blca_multipliers = {key: _blca_multiplier(key) for key in _BLCA_MULTIPLIER_DEFAULTS}

os.makedirs(run_output, exist_ok=True)
config_output_path = run_output + config["run_name"] + "-tax-credit_config.yaml"
shutil.copy(workflow.configfiles[0], config_output_path)

CONFIGFILE = workflow.configfiles[0]


def _row_str(row, field, default=""):
    val = row.get(field, default)
    if pd.isna(val):
        return default
    s = str(val).strip()
    if s.lower() in ("", "nan", "na"):
        return default
    return s


def _row_optional_float(row, field):
    val = row.get(field)
    if pd.isna(val):
        return None
    s = str(val).strip()
    if s.lower() in ("", "nan", "na"):
        return None
    return float(s)


def _config_max_accepts(row):
    """maxaccepts for a consensus job; blank manifest cell falls back to config.

    Kept as a string: consensus-vsearch also accepts ``all``.
    """
    val = _row_str(row, "max_accepts", "")
    if not val:
        cfg_val = config.get("max_accepts", 10)
        if isinstance(cfg_val, (list, tuple)):
            cfg_val = cfg_val[0]
        val = str(cfg_val if cfg_val is not None else 10).strip()
    return val if val.lower() == "all" else str(int(float(val)))


def _config_float(key, default):
    val = config.get(key, default)
    if isinstance(val, (list, tuple)):
        return float(val[0])
    return float(val)


def _assign_params(wildcards):
    row = _manifest_row(wildcards.job_id)
    default_method = config.get("classify_methods") or config.get(
        "classify_method", "naive-bayes"
    )
    if isinstance(default_method, list):
        default_method = default_method[0]
    method = _row_str(row, "classify_method", default_method)
    params = {
        "classify_method": method,
        "fit_params": _row_str(row, "fit_params", ""),
        "classify_params": _row_str(row, "classify_params", ""),
        "bowtie_index_dir": _row_str(row, "bowtie_index_dir", ""),
        # consensus-blast shared database; blank for every other method and for
        # single-combination folds, which build it in-job
        "blastdb_qza": _row_str(row, "blastdb_qza", ""),
        # bt2-blca stage split; blank for every other method and for single-job folds
        "blca_stage": _row_str(row, "blca_stage", ""),
        "blca_shared_dir": _row_str(row, "blca_shared_dir", ""),
        "blca_raw_taxonomy": _row_str(row, "blca_raw_taxonomy", ""),
        "confidence": _row_str(row, "confidence", ""),
        # consensus methods only; overwritten below (blank for every other method)
        "max_accepts": "",
        # revamp (mock-community only) and revamp-custom; blank for every other method
        "revamp_dir": _row_str(row, "revamp_dir", ""),
        "revamp_blastdb": _row_str(row, "revamp_blastdb", ""),
        "revamp_blast_results": _row_str(row, "revamp_blast_results", ""),
        "revamp_blast_mode": _row_str(row, "revamp_blast_mode", "mostEnvOUT"),
        "revamp_query_cov": _row_str(row, "revamp_query_cov", "90"),
        "revamp_taxonomy_cutoffs": _row_str(
            row, "revamp_taxonomy_cutoffs", "97,95,90,80,70,60"
        ),
        # bt2-blca cost multipliers; overwritten below (blank for every other method)
        "blca_opts": "",
    }
    if method in ("consensus-blast", "consensus-vsearch"):
        params["perc_identity"] = _row_optional_float(row, "perc_identity")
        params["query_cov"] = _row_optional_float(row, "query_cov")
        params["min_consensus"] = _row_optional_float(row, "min_consensus")
        params["max_accepts"] = _config_max_accepts(row)
        params["taxa_ranks"] = ""
    elif method in ("revamp", "revamp-custom"):
        # The shared cutoff columns are unused; the revamp_* params above carry what
        # the assignment needs. taxa_ranks is the database's own rank list, which
        # sets how many ranks (and so cutoffs) revamp-custom works with.
        params["perc_identity"] = 0.8
        params["query_cov"] = 0.8
        params["min_consensus"] = 0.51
        params["taxa_ranks"] = _row_str(
            row,
            "taxa_ranks",
            config.get("taxa_ranks", "kingdom,phylum,class,order,family,genus,species"),
        )
    elif method == "bt2-blca":
        # manifest columns are shared; config fallbacks are the required blca_* keys
        params["perc_identity"] = _row_optional_float(row, "perc_identity")
        if params["perc_identity"] is None:
            params["perc_identity"] = _config_float("blca_perc_identity", None)
        params["query_cov"] = _row_optional_float(row, "query_cov")
        if params["query_cov"] is None:
            params["query_cov"] = _config_float("blca_query_cov", None)
        params["min_consensus"] = 0.51
        params["taxa_ranks"] = _row_str(
            row,
            "taxa_ranks",
            config.get("taxa_ranks", "kingdom,phylum,class,order,family,genus,species"),
        )
        params["blca_opts"] = _blca_stage_opts(row)
    else:
        params["perc_identity"] = 0.8
        params["query_cov"] = 0.8
        params["min_consensus"] = 0.51
        params["taxa_ranks"] = ""
    for field in ("perc_identity", "query_cov", "min_consensus"):
        if params.get(field) is None:
            params[field] = _config_float(field, 0.8 if field != "min_consensus" else 0.51)
    return params


def _row_is_true(row, field):
    return _row_str(row, field, "").lower() == "true"


def _blca_stage_opts(row):
    """bt2-blca multiplier flags for the stage(s) this job runs, and only those.

    They are params, so changing one reruns exactly the jobs it affects:
    bowtie2 -k the align stage, bootstraps and muscle iterations the BLCA stage.
    The reformat stage and index-only fit jobs use none of them.
    """
    stage = _row_str(row, "blca_stage", "")
    k = f"--bowtie2-max-hits {_blca_multipliers['blca_bowtie2_max_hits']}"
    blca = (
        f"--blca-bootstraps {_blca_multipliers['blca_bootstraps']} "
        f"--muscle-maxiters {_blca_multipliers['blca_muscle_maxiters']}"
    )
    if stage == "align":
        return k
    if stage == "blca":
        return blca
    if stage == "reformat" or _row_is_true(row, "fit_only") or _row_is_true(row, "trad_fit"):
        return ""
    # single-combination fold: all three stages in one job
    return f"{k} {blca}"


# One rule dispatches every kind of assignment work, but only some of the tools
# it reaches are multithreaded (see scripts/run_taxonomy_assignment_fold.py), so
# the thread count is decided per job rather than per rule. Jobs running a
# single-threaded tool ask for one core, letting Snakemake pack more of them
# into the same --cores.
# Concurrency is bounded by mem_mb (below), not by threads, so a job running a
# single-threaded tool can safely ask for one core.
def _assign_threads(wildcards):
    row = _manifest_row(wildcards.job_id)
    method = _row_str(row, "classify_method", "")
    blca_stage = _row_str(row, "blca_stage", "")
    fit_only = _row_is_true(row, "fit_only")
    trad_fit = _row_is_true(row, "trad_fit")

    # BLCA spreads its queries (muscle + bootstraps each) over classify_threads
    # worker processes, like the align stage before it; the confidence
    # reformat takes no thread option.
    if blca_stage == "reformat":
        return 1
    # fit-classifier-naive-bayes has no thread option. The bt2-blca fit jobs
    # that share this flag build a bowtie2 index, which does.
    if method == "naive-bayes" and (fit_only or trad_fit):
        return 1
    # makeblastdb, the shared consensus-blast fit
    if method == "consensus-blast" and fit_only:
        return 1
    # revamp assignment is the post-BLAST LCA; BLAST is its own fit job, whose
    # btab this row names. Without one it would BLAST here and need the threads.
    if (
        method in ("revamp", "revamp-custom")
        and not fit_only
        and _row_str(row, "revamp_blast_results", "")
    ):
        return 1
    return config["classify_threads"]


# Peak RSS per job type, in MB. These jobs differ by more than an order of
# magnitude, and threads alone do not bound them: the BLCA stage loads the whole
# fold reference into a dict (blca_from_bowtie.py), and a naive-bayes fit over a
# full database with a wide fit_params ngram-range is larger again. Running a
# single-threaded job on every core therefore exhausts memory long before CPU.
#
# Snakemake only ENFORCES this when the run passes --resources mem_mb=<total>;
# without that the numbers are recorded and ignored. See docs/steps/tax_credit.md.
#
# The defaults are starting points, not measurements -- memory scales with the
# reference, so a large database needs larger values. Measure one job with
#   /usr/bin/time -v <the command Snakemake printed>
# and read "Maximum resident set size", then set the matching config key.
# Only one job type's memory scales with its thread count. classify-sklearn
# parallelises through joblib's loky backend (q2_feature_classifier/_skl.py:89),
# which forks n_jobs separate interpreters, each loading its own copy of the
# fitted classifier -- so its reservation is per worker and multiplied below.
# Everything else is thread-parallel inside one process, sharing the one large
# structure it holds: blastn (-num_threads), vsearch (--threads) and bowtie2
# (-p / --threads) all keep a single database/index and give each thread only
# scratch space, so their reservations stay flat. The rest run single-threaded.
# Each extra BLCA worker process: a Python interpreter with Biopython plus one
# muscle subprocess aligning a query with its hits. Workers are sent each
# query's hits rather than a copy of the reference, so this does not scale with
# the database; the reference stays in the one parent counted by mem_mb_blca.
_BLCA_WORKER_MB = 500

_MEM_MB_DEFAULTS = {
    "naive_bayes_fit": 16000,
    "blca": 4000,
    "align": 4000,
    "classify_per_worker": 4000,
    "classify": 4000,
    "makeblastdb": 2000,
    "default": 2000,
}


def _mem_mb(key):
    return int(config.get(f"mem_mb_{key}", _MEM_MB_DEFAULTS[key]))


def _assign_mem_mb(wildcards, threads, attempt):
    """Memory reservation for one assignment job.

    Scaled by *attempt* so a job killed for memory asks for more on each retry
    (snakemake --retries N), which matters because the right number depends on
    a reference size this cannot know. *threads* is used only by
    classify-sklearn, the one job whose memory scales with its own parallelism.
    """
    row = _manifest_row(wildcards.job_id)
    method = _row_str(row, "classify_method", "")
    blca_stage = _row_str(row, "blca_stage", "")
    fit_only = _row_is_true(row, "fit_only")
    trad_fit = _row_is_true(row, "trad_fit")
    # skip_fit means some fit job already built this job's classifier, index or
    # database, so this job only classifies against it. Without it the job
    # builds its own first, which dominates the reservation.
    skip_fit = _row_is_true(row, "skip_fit")

    if blca_stage == "blca":
        base = _mem_mb("blca") + _BLCA_WORKER_MB * (threads - 1)
    elif blca_stage == "align":
        # bowtie2-build over the fold reference
        base = _mem_mb("align")
    elif blca_stage == "reformat":
        # reads one raw-taxonomy TSV; nothing reference-sized
        base = _mem_mb("default")
    elif method == "naive-bayes":
        if skip_fit:
            # classify-sklearn: n_jobs loky workers each hold a classifier
            # copy, and so does the parent that dispatched to them. joblib
            # memmaps large arrays, so some of that is shared and RSS
            # over-reports it -- but the reservation has to assume it is not.
            # n_jobs=1 runs in-process, spawning no workers at all.
            workers = threads + 1 if threads > 1 else 1
            base = _mem_mb("classify_per_worker") * workers
        else:
            base = _mem_mb("naive_bayes_fit")
    elif method == "bt2-blca":
        if fit_only or trad_fit:
            base = _mem_mb("align")
        else:
            # No blca_stage and not a fit job: a fold with one parameter set
            # and one confidence runs all three stages here, so it needs the
            # largest of them.
            base = max(
                _mem_mb("align"), _mem_mb("blca") + _BLCA_WORKER_MB * (threads - 1)
            )
    elif method == "consensus-blast":
        if fit_only:
            base = _mem_mb("makeblastdb")
        elif skip_fit:
            base = _mem_mb("classify")
        else:
            # single-combination fold: makeblastdb then classify, in one job
            base = max(_mem_mb("makeblastdb"), _mem_mb("classify"))
    elif method == "revamp":
        # BLAST against nt, or the post-BLAST LCA over its btab
        base = _mem_mb("default")
    elif method == "revamp-custom":
        # fit: makeblastdb over the fold reference, then BLAST; otherwise the
        # post-BLAST reformat + LCA over its btab
        base = _mem_mb("makeblastdb") if fit_only else _mem_mb("default")
    else:
        # consensus-vsearch: loads the reference and classifies against it
        base = _mem_mb("classify")
    return base * int(attempt)


def _summary_targets():
    names = []
    summary_cfg = config.get("summary_filenames") or {}
    for method in config.get("evaluation_methods", []):
        if method == "mock-community":
            names.extend(
                ["mock_community_metrics.tsv", "mock_community_composition.tsv"]
            )
        elif method in ("cross-validated", "cross-validated-taxa"):
            names.append(
                summary_cfg.get("cross-validated", "evaluate_classification_summary_CV.csv")
            )
        elif method == "novel-taxa":
            names.append(
                summary_cfg.get("novel-taxa", "evaluate_classification_summary_novel.csv")
            )
        elif method == "cross-validated-trad":
            names.append(
                summary_cfg.get(
                    "cross-validated-trad",
                    "evaluate_classification_summary_CV_trad.csv",
                )
            )
        elif method == "self-validated":
            names.append(
                summary_cfg.get(
                    "self-validated",
                    "evaluate_classification_summary_self_validated.csv",
                )
            )
    return names


def _plots_input():
    if config.get("generate_plots", True):
        return plots_done
    return []


_manifest_cache = {"key": None, "rows": None}


def _manifest_row(job_id):
    # Read the manifest once (re-read only if it changes); with thousands of jobs,
    # per-call reads made DAG building / checkpoint updates extremely slow.
    st = os.stat(manifest_fp)
    key = (st.st_mtime_ns, st.st_size)
    if _manifest_cache["key"] != key:
        df = pd.read_csv(manifest_fp, sep="\t")
        _manifest_cache["rows"] = {
            str(jid): row for jid, (_, row) in zip(df["job_id"].astype(str), df.iterrows())
        }
        _manifest_cache["key"] = key
    row = _manifest_cache["rows"].get(str(job_id))
    if row is None:
        raise ValueError(f"job_id not in manifest: {job_id}")
    return row


def _assign_fit_done(wildcards):
    row = _manifest_row(wildcards.job_id)
    fit_job_id = str(row.get("fit_job_id", "") or "").strip()
    if fit_job_id and fit_job_id.lower() != "nan":
        return assign_done_dir + fit_job_id + ".done"
    return []


def _assignment_done_inputs(wildcards):
    ck = checkpoints.tax_credit_prepare_manifest.get().output[0]
    df = pd.read_csv(ck, sep="\t")
    return expand(assign_done_dir + "{job_id}.done", job_id=df["job_id"].astype(str))


def _mock_score_inputs(wildcards):
    """Per-job mock-community scores (jobs that classify, not fit-only jobs)."""
    ck = checkpoints.tax_credit_prepare_manifest.get().output[0]
    df = pd.read_csv(ck, sep="\t", dtype=str, keep_default_na=False)
    if df.empty:
        return []
    is_true = lambda col: df[col].str.lower() == "true"
    jobs = df[
        (df["evaluation_method"] == "mock-community")
        & ~is_true("fit_only")
        & ~is_true("trad_fit")
    ]["job_id"]
    return expand(mock_jobs_dir + "{job_id}.metrics.tsv", job_id=jobs)


rule run_tax_credit:
    input:
        datasets_done,
        manifest_fp,
        _assignment_done_inputs,
        expand(
            summaries_dir + "{summary}",
            summary=_summary_targets(),
        ),
        _plots_input(),
        benchmark_summary_tsv,


rule tax_credit_prepare_datasets:
    output:
        datasets_done,
    benchmark:
        benchmarks_dir + "prepare_datasets.tsv"
    params:
        cfg=CONFIGFILE,
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "python scripts/run_tax_credit.py --config {params.cfg} --phase datasets && "
        "touch {output}"


checkpoint tax_credit_prepare_manifest:
    input:
        datasets_done,
    output:
        manifest_fp,
    benchmark:
        benchmarks_dir + "prepare_manifest.tsv"
    params:
        cfg=CONFIGFILE,
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "python scripts/run_tax_credit.py --config {params.cfg} --phase manifest"


rule tax_credit_assign_fold:
    input:
        # plain path, not checkpoints...get(): these jobs are only requested via
        # the checkpoint-aware aggregators, so re-evaluating each one on checkpoint
        # update is unnecessary (and very slow with thousands of jobs)
        manifest=manifest_fp,
        fit_done=_assign_fit_done,
    output:
        touch(assign_done_dir + "{job_id}.done"),
    benchmark:
        benchmarks_dir + "assign/{job_id}.tsv"
    params:
        row=lambda wildcards: _manifest_row(wildcards.job_id),
        assign=lambda wildcards: _assign_params(wildcards),
    conda:
        "qiime2-amplicon-2024.10"
    threads: _assign_threads
    resources:
        mem_mb=_assign_mem_mb,
    shell:
        """
        set -euo pipefail
        ROW_QUERY=$(echo "{params.row.query_qza}" | tr -d '"')
        ROW_REFS=$(echo "{params.row.ref_seqs}" | tr -d '"')
        ROW_REFT=$(echo "{params.row.ref_taxa}" | tr -d '"')
        ROW_OUT=$(echo "{params.row.output_dir}" | tr -d '"')
        ROW_CONF=$(echo "{params.row.confidence}" | tr -d '"')
        ROW_SKIP=$(echo "{params.row.skip_fit}" | tr -d '"')
        ROW_CLS=$(echo "{params.row.classifier_qza}" | tr -d '"')
        ROW_BT2=$(echo "{params.assign[bowtie_index_dir]}" | tr -d '"')
        ROW_BLASTDB=$(echo "{params.assign[blastdb_qza]}" | tr -d '"')
        ROW_BLCA_STAGE=$(echo "{params.assign[blca_stage]}" | tr -d '"')
        ROW_BLCA_SHARED=$(echo "{params.assign[blca_shared_dir]}" | tr -d '"')
        ROW_BLCA_RAW=$(echo "{params.assign[blca_raw_taxonomy]}" | tr -d '"')
        ROW_RVBLAST=$(echo "{params.assign[revamp_blast_results]}" | tr -d '"')
        TRAD_FIT=$(echo "{params.row.trad_fit}" | tr -d '"')
        ROW_FIT_ONLY=$(echo "{params.row.fit_only}" | tr -d '"')

        mkdir -p "$ROW_OUT"
        EXTRA_CLS=""
        if [ -n "$ROW_CLS" ] && [ "$ROW_CLS" != "nan" ] && [ "$ROW_CLS" != "" ]; then
            EXTRA_CLS="--classifier-qza $ROW_CLS"
        fi
        EXTRA_BT2=""
        if [ -n "$ROW_BT2" ] && [ "$ROW_BT2" != "nan" ] && [ "$ROW_BT2" != "" ]; then
            EXTRA_BT2="--bowtie-index-dir $ROW_BT2"
        fi
        EXTRA_CONSENSUS=""
        if [ "{params.assign[classify_method]}" = "consensus-blast" ] || [ "{params.assign[classify_method]}" = "consensus-vsearch" ]; then
            EXTRA_CONSENSUS="--max-accepts {params.assign[max_accepts]}"
        fi
        EXTRA_BLASTDB=""
        if [ -n "$ROW_BLASTDB" ] && [ "$ROW_BLASTDB" != "nan" ] && [ "$ROW_BLASTDB" != "" ]; then
            EXTRA_BLASTDB="--blastdb-qza $ROW_BLASTDB"
        fi
        EXTRA_REVAMP=""
        if [ "{params.assign[classify_method]}" = "revamp" ]; then
            EXTRA_REVAMP="--revamp-dir {params.assign[revamp_dir]} \
                --revamp-blastdb {params.assign[revamp_blastdb]} \
                --revamp-blast-mode {params.assign[revamp_blast_mode]} \
                --revamp-query-cov {params.assign[revamp_query_cov]} \
                --revamp-taxonomy-cutoffs {params.assign[revamp_taxonomy_cutoffs]}"
            if [ -n "$ROW_RVBLAST" ] && [ "$ROW_RVBLAST" != "nan" ]; then
                EXTRA_REVAMP="$EXTRA_REVAMP --revamp-blast-results $ROW_RVBLAST"
            fi
        elif [ "{params.assign[classify_method]}" = "revamp-custom" ]; then
            EXTRA_REVAMP="--revamp-dir {params.assign[revamp_dir]} \
                --revamp-query-cov {params.assign[revamp_query_cov]} \
                --revamp-taxonomy-cutoffs {params.assign[revamp_taxonomy_cutoffs]}"
            if [ -n "$ROW_RVBLAST" ] && [ "$ROW_RVBLAST" != "nan" ]; then
                EXTRA_REVAMP="$EXTRA_REVAMP --revamp-blast-results $ROW_RVBLAST"
            fi
        fi
        CONF_FLAG=""
        if [ -n "$ROW_CONF" ] && [ "$ROW_CONF" != "nan" ] && [ "$ROW_CONF" != "NA" ]; then
            CONF_FLAG="--confidence $ROW_CONF"
        fi
        SKIP_FLAG=""
        if [ "$ROW_SKIP" = "True" ] || [ "$ROW_SKIP" = "true" ]; then
            SKIP_FLAG="--skip-fit"
        fi
        EXTRA_BLCA=""
        if [ -n "$ROW_BLCA_STAGE" ] && [ "$ROW_BLCA_STAGE" != "nan" ]; then
            EXTRA_BLCA="--blca-stage $ROW_BLCA_STAGE"
            if [ -n "$ROW_BLCA_SHARED" ] && [ "$ROW_BLCA_SHARED" != "nan" ]; then
                EXTRA_BLCA="$EXTRA_BLCA --blca-shared-dir $ROW_BLCA_SHARED"
            fi
            if [ -n "$ROW_BLCA_RAW" ] && [ "$ROW_BLCA_RAW" != "nan" ]; then
                EXTRA_BLCA="$EXTRA_BLCA --blca-raw-taxonomy $ROW_BLCA_RAW"
            fi
        fi
        FIT_ONLY=""
        if [ "$TRAD_FIT" = "True" ] || [ "$TRAD_FIT" = "true" ] || [ "$ROW_FIT_ONLY" = "True" ] || [ "$ROW_FIT_ONLY" = "true" ]; then
            # --blca-stage takes over dispatch, so don't add --fit-only for it.
            # The align stage needs the real query, so don't substitute the
            # reference either (a naive-bayes fit ignores the query). Nor for the
            # REVAMP methods, whose fit job is the BLAST search of the query.
            if [ -z "$EXTRA_BLCA" ]; then
                FIT_ONLY="--fit-only"
                case "{params.assign[classify_method]}" in
                    revamp|revamp-custom) ;;
                    *) ROW_QUERY="$ROW_REFS" ;;
                esac
            fi
        fi

        python scripts/run_taxonomy_assignment_fold.py \
            --query-qza "$ROW_QUERY" \
            --ref-seqs "$ROW_REFS" \
            --ref-taxa "$ROW_REFT" \
            --output-dir "$ROW_OUT" \
            --classify-method {params.assign[classify_method]} \
            $CONF_FLAG \
            --fit-params '{params.assign[fit_params]}' \
            --classify-params '{params.assign[classify_params]}' \
            --classify-threads {threads} \
            --perc-identity {params.assign[perc_identity]} \
            --query-cov {params.assign[query_cov]} \
            --min-consensus {params.assign[min_consensus]} \
            --taxa-ranks '{params.assign[taxa_ranks]}' \
            {params.assign[blca_opts]} \
            $SKIP_FLAG $EXTRA_CLS $EXTRA_BT2 $EXTRA_CONSENSUS $EXTRA_BLASTDB $EXTRA_REVAMP $EXTRA_BLCA $FIT_ONLY

        touch {output}
        """


rule tax_credit_mock_score_job:
    input:
        assign_done_dir + "{job_id}.done",
    output:
        mock_jobs_dir + "{job_id}.metrics.tsv",
        mock_jobs_dir + "{job_id}.composition.tsv",
    benchmark:
        benchmarks_dir + "mock_score/{job_id}.tsv"
    params:
        cfg=CONFIGFILE,
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "python scripts/run_tax_credit.py --config {params.cfg} "
        "--phase mock-evaluate-job --job-id {wildcards.job_id}"


rule tax_credit_evaluate:
    input:
        _assignment_done_inputs,
        _mock_score_inputs,
    output:
        expand(
            summaries_dir + "{summary}",
            summary=_summary_targets(),
        ),
        evaluate_done=touch(evaluate_done),
    benchmark:
        benchmarks_dir + "evaluate.tsv"
    params:
        cfg=CONFIGFILE,
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        """
        python scripts/run_tax_credit.py --config {params.cfg} --phase evaluate
        touch {output.evaluate_done}
        """


rule tax_credit_plot:
    input:
        evaluate_done,
    output:
        plots_done,
    benchmark:
        benchmarks_dir + "plot.tsv"
    params:
        cfg=CONFIGFILE,
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        """
        if [ "{config[generate_plots]}" = "True" ] || [ "{config[generate_plots]}" = "true" ]; then
            python scripts/run_tax_credit.py --config {params.cfg} --phase plot
            touch {output}
        else
            mkdir -p $(dirname {output})
            touch {output}
        fi
        """


# Runs last: every other rule's benchmark file is written when that job
# finishes, so by the time evaluation (and plotting, if enabled) is done, all of
# them exist. Reads them by directory rather than as inputs, because the
# assignment jobs are only known after the manifest checkpoint.
rule tax_credit_benchmark_summary:
    input:
        evaluate_done,
        _plots_input(),
    output:
        benchmark_summary_tsv,
        job_benchmarks_tsv,
    params:
        cfg=CONFIGFILE,
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "python scripts/run_tax_credit.py --config {params.cfg} --phase benchmark-summary"
