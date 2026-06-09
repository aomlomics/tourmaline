## Shared taxonomy assignment rules for Tourmaline taxonomy and tax-credit steps.
## Parent Snakefile must define: output_dir, config, input_repseqs, output_seq, output_tax,
## use_classifier, classify_method, fasta_repseqs (bt2-blca), taxonomy_dir, taxonomy_qza,
## taxonomy_tsv, classifier_qza, fit_params (optional string).

fit_params = config.get("fit_params", "") or ""

if config["classify_method"] == "naive-bayes":
    if use_classifier != "yes":
        rule fit_classifier:
            input:
                refseq=output_seq,
                reftax=output_tax
            output:
                classifier_qza
            conda:
                "qiime2-amplicon-2024.10"
            threads: config["classify_threads"]
            params:
                fitparams=fit_params
            shell:
                "qiime feature-classifier fit-classifier-naive-bayes "
                "--i-reference-reads {input.refseq} "
                "--i-reference-taxonomy {input.reftax} "
                "--o-classifier {output} "
                "{params.fitparams};"
    else:
        rule import_classifier:
            input:
                config["pretrained_classifier"]
            output:
                classifier_qza
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                "ln -sf $(readlink -f {input}) {output}"

    rule feature_classifier_nb:
        input:
            repseqs=input_repseqs,
            classifier=classifier_qza
        output:
            taxonomy_qza,
        params:
            classifyparams=config["classify_params"] if config["classify_params"] else "",
            conf=config["skl_confidence"]
        conda:
            "qiime2-amplicon-2024.10"
        threads: config["classify_threads"]
        shell:
            """
            qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.repseqs} \
            --p-confidence {params.conf} \
            --o-classification {output} \
            --p-n-jobs {threads} \
            {params.classifyparams};
            """

elif config["classify_method"] == "bt2-blca":
    ruleorder: bt2 > export_taxonomy_to_tsv
    bt2_index_path = taxonomy_dir + "bowtie2_index/bowtie2_index.1.bt2"
    if not os.path.exists(bt2_index_path):
        rule bt2_index:
            input:
                refseq=output_seq,
            output:
                directory(taxonomy_dir + "bowtie2_index/")
            params:
                prefix=taxonomy_dir + "bowtie2_index/bowtie2_index"
            conda:
                "bt2-blca"
            threads: config["classify_threads"]
            shell:
                "mkdir -p {output}; "
                "bowtie2-build "
                "--threads {threads} "
                "-f {input.refseq} "
                "{params.prefix}; "

    rule bt2:
        input:
            repseqs=fasta_repseqs,
            index=taxonomy_dir + "bowtie2_index/",
            refseq=output_seq,
            reftax=output_tax
        output:
            sam=taxonomy_dir + config["run_name"] + "_bowtie2_all.sam",
            raw_taxonomy=taxonomy_dir + config["run_name"] + "-raw-taxonomy.tsv",
            taxonomy=taxonomy_tsv,
        params:
            prefix=taxonomy_dir + "bowtie2_index/bowtie2_index",
            classifyparams=config["classify_params"] if config["classify_params"] else "",
            temp_dir=taxonomy_dir + "temp",
            percID=config["perc_identity"],
            querycov=config["query_cov"],
            taxaranks=config["taxa_ranks"],
            conf=config["confidence_thres"]
        conda:
            "bt2-blca"
        threads: config["classify_threads"]
        shell:
            """
            mkdir -p {params.temp_dir};
            bowtie2 -x {params.prefix} -f -U {input.repseqs} -S {params.temp_dir}/end_to_end.sam --no-hd --no-sq --very-sensitive --end-to-end --no-unal -p 120 -k 100 --un {params.temp_dir}/end_to_end_reject.fasta
            bowtie2 -x {params.prefix}  -f -U {params.temp_dir}/end_to_end_reject.fasta -S {params.temp_dir}/local.sam --no-hd --no-sq --very-sensitive --local --no-unal -p 120 -k 100 --un {params.temp_dir}/end_to_end_and_local_reject.fasta
            cat {params.temp_dir}/*.sam > {output.sam}
            python scripts/blca_from_bowtie.py -i {output.sam} -r {input.reftax} -q {input.refseq} -b {params.percID} -l {params.querycov} -p muscle -n 100 -m 1.0 -f 2.5 -g -2 -tr {params.taxaranks} -o {output.raw_taxonomy}
            python scripts/reformat_summary_for_r.py {output.raw_taxonomy} {output.taxonomy} {params.conf} {params.taxaranks}
            if grep -q '^>' {params.temp_dir}/end_to_end_and_local_reject.fasta > /dev/null 2>&1; then
                grep '^>' {params.temp_dir}/end_to_end_and_local_reject.fasta | sed 's/^>//; s/$/\tUnassigned\t0/' >> {output.taxonomy}
            fi
            """

    rule import_taxonomy_to_qza:
        input:
            taxonomy_tsv
        output:
            taxonomy_qza
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "qiime tools import "
            "--type 'FeatureData[Taxonomy]' "
            "--input-format TSVTaxonomyFormat "
            "--input-path {input} "
            "--output-path {output}"

elif config["classify_method"] == "consensus-blast":
    rule feature_classifier_cb:
        input:
            repseqs=input_repseqs,
            refseq=output_seq,
            reftax=output_tax
        output:
            taxonomy_qza,
        params:
            classifyparams=config["classify_params"] if config["classify_params"] else "",
            searchout=taxonomy_dir + "search_results.qza",
            percID=config["perc_identity"],
            querycov=config["query_cov"],
            consensus=config["min_consensus"]
        conda:
            "qiime2-amplicon-2024.10"
        threads: config["classify_threads"]
        shell:
            """
            qiime feature-classifier classify-consensus-blast \
            --i-reference-reads {input.refseq} \
            --i-reference-taxonomy {input.reftax} \
            --i-query {input.repseqs} \
            --p-perc-identity {params.percID} \
            --p-query-cov {params.querycov} \
            --p-min-consensus {params.consensus} \
            --o-classification {output} \
            --o-search-results {params.searchout} \
            {params.classifyparams};
            """

elif config["classify_method"] == "consensus-vsearch":
    rule feature_classifier_cv:
        input:
            repseqs=input_repseqs,
            refseq=output_seq,
            reftax=output_tax,
        output:
            taxonomy_qza,
        params:
            classifyparams=config["classify_params"] if config["classify_params"] else "",
            searchout=taxonomy_dir + "search_results.qza",
            percID=config["perc_identity"],
            querycov=config["query_cov"],
            consensus=config["min_consensus"]
        conda:
            "qiime2-amplicon-2024.10"
        threads: config["classify_threads"]
        shell:
            """
            qiime feature-classifier classify-consensus-vsearch \
            --i-reference-reads {input.refseq} \
            --i-reference-taxonomy {input.reftax} \
            --i-query {input.repseqs} \
            --p-perc-identity {params.percID} \
            --p-query-cov {params.querycov} \
            --p-min-consensus {params.consensus} \
            --o-classification {output} \
            --o-search-results {params.searchout} \
            --p-threads {threads} \
            {params.classifyparams};
            """

rule export_taxonomy_to_tsv:
    input:
        taxonomy_qza
    output:
        taxonomy_tsv
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "--output-format TSVTaxonomyFormat"
