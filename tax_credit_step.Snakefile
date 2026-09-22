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
        "confidence": _row_str(row, "confidence", ""),
        # revamp (mock-community only); blank for every other method
        "revamp_dir": _row_str(row, "revamp_dir", ""),
        "revamp_blastdb": _row_str(row, "revamp_blastdb", ""),
        "revamp_blast_results": _row_str(row, "revamp_blast_results", ""),
        "revamp_blast_mode": _row_str(row, "revamp_blast_mode", "mostEnvOUT"),
        "revamp_query_cov": _row_str(row, "revamp_query_cov", "90"),
        "revamp_taxonomy_cutoffs": _row_str(
            row, "revamp_taxonomy_cutoffs", "97,95,90,80,70,60"
        ),
    }
    if method in ("consensus-blast", "consensus-vsearch"):
        params["perc_identity"] = _row_optional_float(row, "perc_identity")
        params["query_cov"] = _row_optional_float(row, "query_cov")
        params["min_consensus"] = _row_optional_float(row, "min_consensus")
        params["taxa_ranks"] = ""
    elif method == "revamp":
        # nt is the reference database: the shared cutoff columns are unused, and the
        # revamp_* params above carry everything the assignment needs.
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
    else:
        params["perc_identity"] = 0.8
        params["query_cov"] = 0.8
        params["min_consensus"] = 0.51
        params["taxa_ranks"] = ""
    for field in ("perc_identity", "query_cov", "min_consensus"):
        if params.get(field) is None:
            params[field] = _config_float(field, 0.8 if field != "min_consensus" else 0.51)
    return params


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


rule tax_credit_prepare_datasets:
    output:
        datasets_done,
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
    params:
        row=lambda wildcards: _manifest_row(wildcards.job_id),
        assign=lambda wildcards: _assign_params(wildcards),
    conda:
        "qiime2-amplicon-2024.10"
    threads: config["classify_threads"]
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
        fi
        CONF_FLAG=""
        if [ -n "$ROW_CONF" ] && [ "$ROW_CONF" != "nan" ] && [ "$ROW_CONF" != "NA" ]; then
            CONF_FLAG="--confidence $ROW_CONF"
        fi
        SKIP_FLAG=""
        if [ "$ROW_SKIP" = "True" ] || [ "$ROW_SKIP" = "true" ]; then
            SKIP_FLAG="--skip-fit"
        fi
        FIT_ONLY=""
        if [ "$TRAD_FIT" = "True" ] || [ "$TRAD_FIT" = "true" ] || [ "$ROW_FIT_ONLY" = "True" ] || [ "$ROW_FIT_ONLY" = "true" ]; then
            FIT_ONLY="--fit-only"
            ROW_QUERY="$ROW_REFS"
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
            $SKIP_FLAG $EXTRA_CLS $EXTRA_BT2 $EXTRA_REVAMP $FIT_ONLY

        touch {output}
        """


rule tax_credit_mock_score_job:
    input:
        assign_done_dir + "{job_id}.done",
    output:
        mock_jobs_dir + "{job_id}.metrics.tsv",
        mock_jobs_dir + "{job_id}.composition.tsv",
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
