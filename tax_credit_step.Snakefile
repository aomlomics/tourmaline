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
fit_params = config.get("fit_params", "") or ""
classify_params = config.get("classify_params", "") or ""


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
        run_output + ".collected.done",
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
        TRAD_FIT=$(echo "{params.row.trad_fit}" | tr -d '"')
        ROW_FIT_ONLY=$(echo "{params.row.fit_only}" | tr -d '"')

        mkdir -p "$ROW_OUT"
        EXTRA_CLS=""
        if [ -n "$ROW_CLS" ] && [ "$ROW_CLS" != "nan" ] && [ "$ROW_CLS" != "" ]; then
            EXTRA_CLS="--classifier-qza $ROW_CLS"
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
            --classify-method {config[classify_method]} \
            --confidence "$ROW_CONF" \
            --fit-params '{fit_params}' \
            --classify-params '{classify_params}' \
            --classify-threads {threads} \
            --perc-identity {config[perc_identity]} \
            --query-cov {config[query_cov]} \
            --min-consensus {config[min_consensus]} \
            $SKIP_FLAG $EXTRA_CLS $FIT_ONLY

        touch {output}
        """


rule tax_credit_collect_results:
    input:
        _assignment_done_inputs,
    output:
        touch(run_output + ".collected.done"),
    params:
        cfg=CONFIGFILE,
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "python scripts/run_tax_credit.py --config {params.cfg} --phase collect && "
        "touch {output}"


rule tax_credit_evaluate:
    input:
        run_output + ".collected.done",
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
