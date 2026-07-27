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
plots_done = run_output + config.get("plots_subdir", "plots") + "/.done"

os.makedirs(run_output, exist_ok=True)
config_output_path = run_output + config["run_name"] + "-tax-credit_config.yaml"
shutil.copy(workflow.configfiles[0], config_output_path)

CONFIGFILE = workflow.configfiles[0]


def _row_str(row, field, default=""):
    val = row.get(field, default)
    if pd.isna(val):
        return default
    return str(val)


def _config_float(key, default):
    val = config.get(key, default)
    if isinstance(val, (list, tuple)):
        return float(val[0])
    return float(val)


def _row_float(row, field, default):
    val = row.get(field, default)
    if pd.isna(val):
        return _config_float(field, default)
    return float(val)


def _assign_params(wildcards):
    row = _manifest_row(wildcards.job_id)
    default_method = config.get("classify_methods") or config.get(
        "classify_method", "naive-bayes"
    )
    if isinstance(default_method, list):
        default_method = default_method[0]
    return {
        "classify_method": _row_str(row, "classify_method", default_method),
        "fit_params": _row_str(row, "fit_params", config.get("fit_params", "") or ""),
        "classify_params": _row_str(
            row, "classify_params", config.get("classify_params", "") or ""
        ),
        "perc_identity": _row_float(row, "perc_identity", 0.8),
        "query_cov": _row_float(row, "query_cov", 0.8),
        "min_consensus": _row_float(row, "min_consensus", 0.51),
        "taxa_ranks": _row_str(
            row,
            "taxa_ranks",
            config.get("taxa_ranks", "kingdom,phylum,class,order,family,genus,species"),
        ),
        "bowtie_index_dir": _row_str(row, "bowtie_index_dir", ""),
    }


def _summary_targets():
    names = []
    summary_cfg = config.get("summary_filenames") or {}
    for method in config.get("evaluation_methods", []):
        if method == "mock-community":
            names.append(
                summary_cfg.get("mock-community", "mock_evaluation_summary.tsv")
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


def _manifest_row(job_id):
    df = pd.read_csv(manifest_fp, sep="\t")
    rows = df[df["job_id"].astype(str) == str(job_id)]
    if rows.empty:
        raise ValueError(f"job_id not in manifest: {job_id}")
    return rows.iloc[0]


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
        manifest=lambda wildcards: checkpoints.tax_credit_prepare_manifest.get().output[0],
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
            --confidence "$ROW_CONF" \
            --fit-params '{params.assign[fit_params]}' \
            --classify-params '{params.assign[classify_params]}' \
            --classify-threads {threads} \
            --perc-identity {params.assign[perc_identity]} \
            --query-cov {params.assign[query_cov]} \
            --min-consensus {params.assign[min_consensus]} \
            --taxa-ranks '{params.assign[taxa_ranks]}' \
            $SKIP_FLAG $EXTRA_CLS $EXTRA_BT2 $FIT_ONLY

        touch {output}
        """


rule tax_credit_evaluate:
    input:
        _assignment_done_inputs,
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
