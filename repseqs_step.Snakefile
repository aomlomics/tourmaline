## STILL NEED TO ADD some RULES FOR FILTERING SEQUENCES

## Snakefile for repseqs step of Tourmaline V2 pipeline
output_dir = config["output_dir"]+"/"

if config["sample_metadata_file"] != None:
    use_metadata="yes"
else:
    use_metadata="no"

# set run name
if config["sample_run_name"] != None:
    sample_run_name=config["sample_run_name"]
    input_fastq=output_dir+sample_run_name+"-samples/"+sample_run_name+"_fastq_pe.qza"
elif config["fastq_qza_file"] != None:
    input_fastq=config["fastq_qza_file"]
else:
    input_fastq=output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza"

# set Filtering
if config["to_filter"] == True:
    temp_table = output_dir+config["run_name"]+"-repseqs/temp-table.qza"
    temp_repseqs = output_dir+config["run_name"]+"-repseqs/temp-repseqs.qza"
else:
    temp_table = output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza"
    temp_repseqs = output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.qza"


## Master RULES
rule run_dada2_pe_denoise:
    """Run paired end dada2"""
    input:
       #output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
        output_dir+config["run_name"]+"-repseqs/stats/table_summary.qzv",
        output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.tsv",
        output_dir+config["run_name"]+"-repseqs/stats/table_summary_samples.txt",
        output_dir+config["run_name"]+"-repseqs/stats/table_summary_features.txt",
        output_dir+config["run_name"]+"-repseqs/stats/repseqs.qzv",
        output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths_describe.md"
        #config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.fasta",
        #config["run_name"]+"-repseqs/"+config["run_name"]+"-table.biom",
        #config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs-stats.tsv",
        # add figures here

#rule run_dada2_se_denoise:

#rule run_deblur_se_denoise:

rule denoise_dada2_pe:
    input:
        input_fastq
    params:
        trunclenf=config["dada2pe_trunc_len_f"],
        trunclenr=config["dada2pe_trunc_len_r"],
        trimleftf=config["dada2pe_trim_left_f"],
        trimleftr=config["dada2pe_trim_left_r"],
        maxeef=config["dada2pe_max_ee_f"],
        maxeer=config["dada2pe_max_ee_r"],
        truncq=config["dada2pe_trunc_q"],
        poolingmethod=config["dada2pe_pooling_method"],        
        chimeramethod=config["dada2pe_chimera_method"],
        minfoldparentoverabundance=config["dada2pe_min_fold_parent_over_abundance"],
        nreadslearn=config["dada2pe_n_reads_learn"],
        hashedfeatureids=config["dada2pe_hashed_feature_ids"],
        minlength=config["repseq_min_length"],
        maxlength=config["repseq_max_length"],
        minabund=config["repseq_min_abundance"],
        minprev=config["repseq_min_prevalence"],
        to_filter=config["to_filter"]
    output:
        table=temp_table,
        repseqs=temp_repseqs,
        stats=output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qza",
    conda:
        "qiime2-2023.5"
    threads: config["asv_threads"]
    shell:
        """
        qiime dada2 denoise-paired \
        --i-demultiplexed-seqs {input[0]} \
        --p-trunc-len-f {params.trunclenf} \
        --p-trunc-len-r {params.trunclenr} \
        --p-trim-left-f {params.trimleftf} \
        --p-trim-left-r {params.trimleftr} \
        --p-max-ee-f {params.maxeef} \
        --p-max-ee-r {params.maxeer} \
        --p-trunc-q {params.truncq} \
        --p-pooling-method {params.poolingmethod} \
        --p-chimera-method {params.chimeramethod} \
        --p-min-fold-parent-over-abundance {params.minfoldparentoverabundance} \
        --p-n-reads-learn {params.nreadslearn} \
        --p-n-threads {threads} \
        {params.hashedfeatureids} \
        --o-table {output.table} \
        --o-representative-sequences {output.repseqs} \
        --o-denoising-stats {output.stats} \
        --verbose  
        """

# FILTER
if config["to_filter"] == True:
    rule filter_sequences:
        input:
            table=temp_table,
            repseqs=temp_repseqs,
             #repseqstofilter="00-data/repseqs_to_filter_{method}.tsv",
             #samplestofilter="00-data/samples_to_filter_{method}.tsv",
             #metadata="00-data/metadata.tsv"
        params:
            minlength=config["repseq_min_length"],
            maxlength=config["repseq_max_length"],
            minabund=config["repseq_min_abundance"],
            minprev=config["repseq_min_prevalence"],
        output:
            output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
            output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.qza",
        conda:
            "qiime2-2023.5"
        shell:
            # FILTER SEQUENCES BY LENGTH
            "qiime feature-table filter-seqs "
            "--i-data {input.repseqs} "
            "--m-metadata-file {input.repseqs} "
            "--p-where 'length(sequence) >= {params.minlength} AND length(sequence) <= {params.maxlength}' "
            "--o-filtered-data temp_repseqs1.qza; "
            "/bin/rm {input.repseqs}; "
            # FILTER TABLE BY FEATURE IDS
            "qiime feature-table filter-features "
            "--i-table {input.table} "
            "--m-metadata-file temp_repseqs1.qza "
            "--o-filtered-table temp_table.qza; "
            "/bin/rm {input.table}; "
            # FILTER TABLE BY ABUNDANCE & PREVALENCE
            "qiime feature-table filter-features-conditionally "
            "--i-table temp_table.qza "
            "--p-abundance {params.minabund} "
            "--p-prevalence {params.minprev} "
            "--o-filtered-table {output[0]}; "
            # FILTER SEQUENCES USING TABLE
            "qiime feature-table filter-seqs "
            "--i-data temp_repseqs1.qza "
            "--i-table {output[0]} "
            "--p-no-exclude-ids "
            "--o-filtered-data {output[1]}; "
            # REMOVE TEMP FILES
            "/bin/rm temp_repseqs1.qza; "
            "/bin/rm temp_table.qza; "

# RULES: SUMMARIZE FEATURE TABLE -----------------------------------------------

rule summarize_feature_table:
    input:
        table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
    output:
        output_dir+config["run_name"]+"-repseqs/stats/table_summary.qzv"
    params:
        metadata=config["sample_metadata_file"]
    conda:
        "qiime2-2023.5"
    threads: config["asv_threads"]
    shell:
        """
        if [ {use_metadata} = yes ]; then
            qiime feature-table summarize \
            --i-table {input.table} \
            --m-sample-metadata-file {params.metadata} \
            --o-visualization {output}
        else
            qiime feature-table summarize \
            --i-table {input.table} \
            --o-visualization {output}
        fi
        """

rule summarize_repseqs:
    input:
        stats=output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qza"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        "qiime metadata tabulate "
        "--m-input-file {input.stats} "
        "--o-visualization {output}"

rule export_repseqs_summary_to_tsv:
    input:
        output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qzv"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.tsv"
    conda:
        "qiime2-2023.5"
    shell:
        "unzip -qq -o {input} -d temp0; "
        "mv temp0/*/data/metadata.tsv {output}; "
        "/bin/rm -r temp0"


rule export_table_to_biom:
    input:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza"
    output:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.biom"
    conda:
        "qiime2-2023.5"
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "--output-format BIOMV210Format"

rule summarize_biom_samples:
    input:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.biom"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/table_summary_samples.txt"
    conda:
        "qiime2-2023.5"
    shell:
        "biom summarize-table "
        "--input-fp {input} "
        "--output-fp {output}; "
        "cat {output} | sed 's/observation/feature/g' | sed 's/.000$//' > temp; "
        "mv temp {output}"

rule summarize_biom_features:
    input:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.biom"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/table_summary_features.txt"
    conda:
        "qiime2-2023.5"
    threads: config["asv_threads"]
    shell:
        "biom summarize-table "
        "--observations "
        "--input-fp {input} "
        "--output-fp {output}; "
        "cat {output} | sed 's/observation/feature/g' | sed 's|Counts/sample|Counts/feature|g' | sed 's/.000$//' > temp; "
        "mv temp {output}"

rule visualize_repseqs:
    input:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.qza"
    output:
       output_dir+config["run_name"]+"-repseqs/stats/repseqs.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["asv_threads"]
    shell:
        "qiime feature-table tabulate-seqs "
        "--i-data {input} "
        "--o-visualization {output}"

rule export_repseqs_to_fasta:
    input:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.qza"
    output:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.fasta"
    conda:
        "qiime2-2023.5"
    threads: config["asv_threads"]
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "--output-format DNAFASTAFormat"

rule repseqs_lengths:
    input:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.fasta"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["asv_threads"]
    shell:
        "perl scripts/fastaLengths.pl {input} > {output}"

rule repseqs_lengths_describe:
    input:
        output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths.tsv"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths_describe.md"
    conda:
        "qiime2-2023.5"
    threads: config["asv_threads"]
    shell:
        "python scripts/repseqs_lengths_describe.py {input} {output}"
