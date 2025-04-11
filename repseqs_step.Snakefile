## STILL NEED TO ADD some RULES FOR FILTERING SEQUENCES

# NEED TO ADD TABLE TSV OUTPUT

## Snakefile for repseqs step of Tourmaline V2 pipeline
output_dir = config["output_dir"]+"/"

if config["sample_metadata_file"] != None:
    use_metadata="yes"
else:
    use_metadata="no"

# set run name
if config["sample_run_name"] != None:
    sample_run_name=config["sample_run_name"]
    input_fastq=output_dir+sample_run_name+"-qaqc/"+sample_run_name+"_fastq.qza"
elif config["fastq_qza_file"] != None:
    input_fastq=config["fastq_qza_file"]
else:
    input_fastq=output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq.qza"

# set Filtering
if config["to_filter"] == True:
    temp_table = output_dir+config["run_name"]+"-repseqs/temp-table.qza"
    temp_repseqs = output_dir+config["run_name"]+"-repseqs/temp-repseqs.qza"
else:
    temp_table = output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza"
    temp_repseqs = output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.qza"

## Master RULES

# Helper function to determine required input files
def get_required_inputs(config):
    required = []
    
    # Base outputs always required
    required.append(output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.tsv")
    required.append(output_dir + config["run_name"] + "-repseqs/stats/table_summary.qzv")
    required.append(output_dir+config["run_name"]+"-repseqs/stats/table_summary_samples.txt")
    required.append(output_dir+config["run_name"]+"-repseqs/stats/table_summary_features.txt")
    required.append(output_dir+config["run_name"]+"-repseqs/stats/repseqs.qzv")
    required.append(output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths_describe.md")
    required.append(output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.tsv")
    
    # Conditional outputs based on config
    if config.get("plot_diversity", True):
        required.append(output_dir+config["run_name"]+"-repseqs/stats/alpha_rarefaction.qzv")
        required.append(output_dir+config["run_name"]+"-repseqs/stats/rarefied_table.qza")
        
    return required


rule run_denoise:
    """Run paired end dada2"""
    input:
        get_required_inputs(config)
        #output_dir+config["run_name"]+"-repseqs/stats/table_summary.qzv",
        #output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.tsv",
        #output_dir+config["run_name"]+"-repseqs/stats/table_summary_samples.txt",
        #output_dir+config["run_name"]+"-repseqs/stats/table_summary_features.txt",
        #output_dir+config["run_name"]+"-repseqs/stats/repseqs.qzv",
        #output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths_describe.md",
        #output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.tsv"
        # add figures here

if config["asv_method"] == "dada2pe":
    rule denoise_dada2_pe:
        input:
            input_fastq
        params:
            trunclenf=config["dada2_trunc_len_f"],
            trunclenr=config["dada2pe_trunc_len_r"],
            trimleftf=config["dada2_trim_left_f"],
            trimleftr=config["dada2pe_trim_left_r"],
            maxeef=config["dada2_max_ee_f"],
            maxeer=config["dada2pe_max_ee_r"],
            truncq=config["dada2_trunc_q"],
            poolingmethod=config["dada2_pooling_method"],        
            chimeramethod=config["dada2_chimera_method"],
            minfoldparentoverabundance=config["dada2_min_fold_parent_over_abundance"],
            nreadslearn=config["dada2_n_reads_learn"],
            hashedfeatureids=config["dada2_hashed_feature_ids"]
        output:
            table=temp_table,
            repseqs=temp_repseqs,
            stats=output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qza",
        conda:
            "qiime2-amplicon-2024.10"
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
elif config["asv_method"] == "dada2se":
    rule denoise_dada2_se:
        input:
            input_fastq
        params:
            trunclenf=config["dada2_trunc_len_f"],
            trimleftf=config["dada2_trim_left_f"],
            maxeef=config["dada2_max_ee_f"],
            truncq=config["dada2_trunc_q"],
            poolingmethod=config["dada2_pooling_method"],        
            chimeramethod=config["dada2_chimera_method"],
            minfoldparentoverabundance=config["dada2_min_fold_parent_over_abundance"],
            nreadslearn=config["dada2_n_reads_learn"],
            hashedfeatureids=config["dada2_hashed_feature_ids"]
        output:
            table=temp_table,
            repseqs=temp_repseqs,
            stats=output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qza",
        conda:
            "qiime2-amplicon-2024.10"
        threads: config["asv_threads"]
        shell:
            """
            qiime dada2 denoise-single \
            --i-demultiplexed-seqs {input[0]} \
            --p-trunc-len {params.trunclenf} \
            --p-trim-left {params.trimleftf} \
            --p-max-ee {params.maxeef} \
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
else:
    raise ValueError("Invalid ASV method specified")

#repseq_min_samples: 0 #qiime feature-table filter-features --p-min-samples

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
            "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
    threads: config["asv_threads"]
    shell:
        """
        if [ {use_metadata} == "yes" ]; then
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
        "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
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
        "qiime2-amplicon-2024.10"
    threads: config["asv_threads"]
    shell:
        "perl scripts/fastaLengths.pl {input} > {output}"

rule repseqs_lengths_describe:
    input:
        output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths.tsv"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/repseqs_lengths_describe.md"
    conda:
        "qiime2-amplicon-2024.10"
    threads: config["asv_threads"]
    shell:
        "python scripts/repseqs_lengths_describe.py {input} {output}"

rule export_biom_tsv:
    input:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.biom"
    output:
        output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.tsv"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "biom convert "
        "-i {input} "
        "-o {output} "
        "--to-tsv"

rule diversity_alpha_rarefaction:
    input:
        table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
    params:
        maxdepth=config["alpha_max_depth"],
        metadata=config["sample_metadata_file"]
    output:
        output_dir+config["run_name"]+"-repseqs/stats/alpha_rarefaction.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        """
        if [ {use_metadata} == "yes" ]; then
            qiime diversity alpha-rarefaction \
            --i-table {input.table} \
            --p-max-depth {params.maxdepth} \
            --p-metrics observed_features \
            --p-metrics shannon \
            --p-metrics pielou_e \
            --m-metadata-file {params.metadata} \
            --o-visualization {output}
        else
            qiime diversity alpha-rarefaction \
            --i-table {input.table} \
            --p-max-depth {params.maxdepth} \
            --p-metrics observed_features \
            --p-metrics shannon \
            --p-metrics pielou_e \
            --o-visualization {output}
        fi
        """
rule diversity_core_metrics:
    input:
        table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
    params:
        samplingdepth=config["core_sampling_depth"],
        metadata=config["sample_metadata_file"]
    output:
        rarefiedtable=output_dir+config["run_name"]+"-repseqs/stats/rarefied_table.qza",
        observedfeaturesvector=output_dir+config["run_name"]+"-repseqs/stats/observed_features_vector.qza",
        shannonvector=output_dir+config["run_name"]+"-repseqs/stats/shannon_vector.qza",
        evennessvector=output_dir+config["run_name"]+"-repseqs/stats/evenness_vector.qza",
        jaccarddistancematrix=output_dir+config["run_name"]+"-repseqs/stats/jaccard_distance_matrix.qza",
        braycurtisdistancematrix=output_dir+config["run_name"]+"-repseqs/stats/bray_curtis_distance_matrix.qza",
        jaccardpcoaresults=output_dir+config["run_name"]+"-repseqs/stats/jaccard_pcoa_results.qza",
        braycurtispcoaresults=output_dir+config["run_name"]+"-repseqs/stats/bray_curtis_pcoa_results.qza",
        jaccardemperor=output_dir+config["run_name"]+"-repseqs/stats/jaccard_emperor.qzv",
        braycurtisemperor=output_dir+config["run_name"]+"-repseqs/stats/bray_curtis_emperor.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    threads: config["asv_threads"]
    shell:
        "qiime diversity core-metrics "
        "--i-table {input.table} "
        "--p-sampling-depth {params.samplingdepth} "
        "--m-metadata-file {params.metadata} "
        "--o-rarefied-table {output.rarefiedtable} "
        "--o-observed-features-vector {output.observedfeaturesvector} "
        "--o-shannon-vector {output.shannonvector} "
        "--o-evenness-vector {output.evennessvector} "
        "--o-jaccard-distance-matrix {output.jaccarddistancematrix} "
        "--o-bray-curtis-distance-matrix {output.braycurtisdistancematrix} "
        "--o-jaccard-pcoa-results {output.jaccardpcoaresults} "
        "--o-bray-curtis-pcoa-results {output.braycurtispcoaresults} "
        "--o-jaccard-emperor {output.jaccardemperor} "
        "--o-bray-curtis-emperor {output.braycurtisemperor} "
        "--p-n-jobs {threads}"
