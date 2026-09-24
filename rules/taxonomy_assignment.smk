## Shared taxonomy assignment rules, included by taxonomy_step.Snakefile.
## Parent Snakefile must define: output_dir, config, input_repseqs, output_seq, output_tax,
## use_classifier, classify_method, fasta_repseqs (bt2-blca), input_table (revamp),
## taxonomy_dir, taxonomy_qza, taxonomy_tsv, classifier_qza, fit_params (optional string).

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

elif config["classify_method"] == "revamp":
    ## REVAMP (https://github.com/McAllister-NOAA/REVAMP) BLASTs ASVs against a local
    ## NCBI nt database and merges the best hits to their lowest common ancestor.
    ## Tourmaline calls REVAMP's taxonomy scripts directly (see
    ## scripts/run_revamp_taxonomy.sh) rather than revamp.sh, which also runs cutadapt,
    ## DADA2, figures and interactive prompts. Nothing in the REVAMP clone is modified.
    ruleorder: revamp_to_taxonomy_tsv > export_taxonomy_to_tsv

    # REVAMP's scripts use relative paths between these directories, so its output
    # layout is mirrored under the taxonomy run directory.
    revamp_dir = taxonomy_dir + "revamp/"
    revamp_fasta = revamp_dir + "dada2/ASVs.fa"
    revamp_counts = revamp_dir + "dada2/ASVs_counts.tsv"
    revamp_check = revamp_dir + "input_check.txt"
    revamp_btab = revamp_dir + "blast_results/ASV_blastn_nt.btab"
    revamp_formatted_blast = revamp_dir + "blast_results/ASV_blastn_nt_formatted.txt"
    revamp_asv_taxonomy = revamp_dir + "ASV2Taxonomy/" + config["run_name"] + "_asvTaxonomyTable.txt"

    # BLASTing nt usually happens on the machine that holds the database, so BLAST
    # results may be supplied instead of being produced here.
    revamp_blast_results = config.get("revamp_blast_results") or None

    rule revamp_prep:
        input:
            repseqs=input_repseqs,
            table=input_table,
            btab=revamp_blast_results if revamp_blast_results else []
        output:
            fasta=revamp_fasta,
            counts=revamp_counts,
            check=revamp_check
        params:
            btabarg=("--blast-results " + revamp_blast_results) if revamp_blast_results else ""
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            """
            qiime tools export \
            --input-path {input.repseqs} \
            --output-path {output.fasta} \
            --output-format DNAFASTAFormat
            python scripts/revamp_check_inputs.py \
            --repseqs-fasta {output.fasta} \
            --table {input.table} \
            {params.btabarg} \
            --output {output.check}
            printf 'x\tplaceholder_sample\n' > {output.counts}
            grep '^>' {output.fasta} | sed 's/^>//' | awk '{{print $1"\t1"}}' >> {output.counts}
            """

    if revamp_blast_results:
        rule revamp_stage_blast:
            input:
                btab=revamp_blast_results,
                check=revamp_check
            output:
                revamp_btab
            shell:
                "cp {input.btab} {output}"
    else:
        rule revamp_blast:
            input:
                fasta=revamp_fasta,
                check=revamp_check
            output:
                revamp_btab
            params:
                revampdir=config["revamp_dir"],
                workdir=revamp_dir,
                blastdb=config["revamp_blastdb"],
                blastmode=config.get("revamp_blast_mode") or "mostEnvOUT"
            conda:
                "revamp"
            threads: config["classify_threads"]
            shell:
                "bash scripts/run_revamp_taxonomy.sh "
                "--mode blast "
                "--revamp-dir {params.revampdir} "
                "--workdir {params.workdir} "
                "--blastdb {params.blastdb} "
                "--blast-mode {params.blastmode} "
                "--threads {threads}"

    rule revamp_assign:
        input:
            btab=revamp_btab,
            fasta=revamp_fasta,
            counts=revamp_counts
        output:
            asvtaxonomy=revamp_asv_taxonomy,
            formattedblast=revamp_formatted_blast
        params:
            revampdir=config["revamp_dir"],
            workdir=revamp_dir,
            blastdb=config["revamp_blastdb"],
            runname=config["run_name"],
            querycov=config.get("revamp_query_cov") or 90,
            cutoffs=config.get("revamp_taxonomy_cutoffs") or "97,95,90,80,70,60"
        conda:
            "revamp"
        shell:
            "bash scripts/run_revamp_taxonomy.sh "
            "--mode assign "
            "--revamp-dir {params.revampdir} "
            "--workdir {params.workdir} "
            "--blastdb {params.blastdb} "
            "--run-name {params.runname} "
            "--query-cov {params.querycov} "
            "--cutoffs {params.cutoffs}"

    rule revamp_to_taxonomy_tsv:
        input:
            asvtaxonomy=revamp_asv_taxonomy,
            formattedblast=revamp_formatted_blast,
            fasta=revamp_fasta
        output:
            taxonomy_tsv
        params:
            taxaranks=config["taxa_ranks"]
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "python scripts/revamp_to_qiime_taxonomy.py "
            "--asv-taxonomy-table {input.asvtaxonomy} "
            "--formatted-blast {input.formattedblast} "
            "--repseqs-fasta {input.fasta} "
            "--output {output} "
            "--taxaranks {params.taxaranks}"

    rule revamp_import_taxonomy:
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
