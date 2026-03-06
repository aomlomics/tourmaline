## Tourmaline repseqs Snakemake workflow.
## Invoked via `tourmaline.sh --step repseqs`; denoises reads, optionally filters
## features, optionally runs diversity statistics, and exports summaries/visualizations for representative sequences.

import shutil

## STILL NEED TO ADD some RULES FOR FILTERING SEQUENCES
output_dir = config["output_dir"]+"/"

# Copy config file to output directory
config_output_path = output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs_config.yaml"
os.makedirs(os.path.dirname(config_output_path), exist_ok=True)
shutil.copy(workflow.configfiles[0], config_output_path)

# set run name
if config["qaqc_run_name"] != None:
    qaqc_run_name=config["qaqc_run_name"]
    input_fastq=output_dir+qaqc_run_name+"-qaqc/"+qaqc_run_name+"_fastq.qza"
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
    
    # Add stats file based on asv_method
    if config.get("asv_method") in ["dada2pe", "dada2se"]:
        required.append(output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.tsv")
    elif config.get("asv_method") == "deblur":
        required.append(output_dir+config["run_name"]+"-repseqs/stats/deblur_stats.tsv")
    
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


# Aggregate rule: ensures denoising outputs plus requested visualizations exist.
rule run_denoise:
    """Run denoising step"""
    input:
        get_required_inputs(config)


# set metadata file option
if config["sample_metadata_file"] != None:
    use_metadata="yes"

    rule cp_metadata:
        input:
            config["sample_metadata_file"]
        output: 
            output_dir+config["run_name"]+"-repseqs/stats/metadata_used.tsv"
        conda: "qiime2-amplicon-2024.10"
        threads: config["asv_threads"]
        shell: "cp {input} {output}" 
else:
    use_metadata="no"
    # Generate minimal metadata when none is supplied by the user.
    rule autogenerate_metadata:
        input:
            input_fastq
        output:
            output_dir+config["run_name"]+"-repseqs/stats/metadata_used.tsv"
        conda: "qiime2-amplicon-2024.10"
        threads: config["asv_threads"]
        shell:
            """
            qiime tools export --input-path {input} --output-path exported_data
            awk -F',' 'NR>1 && !seen[$1]++ {{print $1}}' exported_data/MANIFEST | awk 'BEGIN{{print "sample_name"}} {{print}}' > {output}
            /bin/rm -r exported_data
            """


if config["asv_method"] == "dada2pe":
    print(f"Running DADA2 paired-end.\n\n")
    # Run paired-end DADA2 denoising with validation against read lengths.
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
            export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=; 
            # Only check if trunclen > 0
            if ([ {params.trunclenf} -gt 0 ] || [ {params.trunclenr} -gt 0 ]); then
                echo "Checking that truncation lengths are less than maximum read length."
                qiime demux summarize \
                --i-data {input[0]} \
                --o-visualization temp0-fastq;
                unzip -qq -o temp0-fastq.qzv -d temp0
                fwdresult=$(python ./scripts/get_last_column_int.py temp0/*/data/forward-seven-number-summaries.tsv)
                revresult=$(python ./scripts/get_last_column_int.py temp0/*/data/reverse-seven-number-summaries.tsv)

                # Check if the result is less than trunclen
                if [ "$fwdresult" -le {params.trunclenf} ]; then
                    echo "ERROR: Forward read length ($fwdresult) is less than dada2_trunc_len_f ({params.trunclenf}). Fix your config file so that dada2_trunc_len_f is less than the read length."
                    /bin/rm -r temp0
                    /bin/rm -r temp0-fastq.qzv
                    exit 1
                else
                    echo "SUCCESS: Forward read length ($fwdresult) is more than dada2_trunc_len_f ({params.trunclenf})."
                fi;
                if [ "$revresult" -le {params.trunclenr} ]; then
                    echo "ERROR: Reverse read length ($revresult) is less than dada2pe_trunc_len_r ({params.trunclenr}). Fix your config file so that dada2pe_trunc_len_r is less than the read length."
                    /bin/rm -r temp0
                    /bin/rm -r temp0-fastq.qzv
                    exit 1
                else
                    echo "SUCCESS: Reverse read length ($revresult) is less than dada2_trunc_len_r ({params.trunclenr})."
                fi;
                /bin/rm -r temp0
                /bin/rm -r temp0-fastq.qzv
            fi;

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
    print(f"Running DADA2 single-end.\n\n")
    # Run single-end DADA2 denoising using configured trimming and error params.
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
            export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=;
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
elif config["asv_method"] == "deblur":
    print(f"Running Deblur.\n\n")
    # Run Deblur workflow leveraging reference sequences for positive filtering.
    rule denoise_deblur:
        input:
            input_fastq,
            reference_seqs=config["reference_seqs"]
        params:
            trim_length=config["deblur_trim_length"],
            trim_left=config["deblur_trim_left"],
            mean_error=config["deblur_mean_error"],
            min_reads=config["deblur_min_reads"],
            min_size=config["deblur_min_size"],
            indel_max=config["deblur_indel_max"]
        output:
            table=temp_table,
            repseqs=temp_repseqs,
            stats=output_dir+config["run_name"]+"-repseqs/stats/deblur_stats.qza",
        conda:
            "qiime2-amplicon-2024.10"
        threads: config["asv_threads"]
        shell:
            """
            export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=;
            qiime deblur denoise-other \
            --i-demultiplexed-seqs {input[0]} \
            --i-reference-seqs {input.reference_seqs} \
            --p-trim-length {params.trim_length} \
            --p-left-trim-len {params.trim_left} \
            --p-mean-error {params.mean_error} \
            --p-min-reads {params.min_reads} \
            --p-min-size {params.min_size} \
            --p-indel-max {params.indel_max} \
            --p-hashed-feature-ids \
            --p-sample-stats \
            --p-jobs-to-start {threads} \
            --o-table {output.table} \
            --o-representative-sequences {output.repseqs} \
            --o-stats {output.stats} \
            --verbose  \
            """
else:
    raise ValueError("Invalid ASV method specified")

#repseq_min_samples: 0 #qiime feature-table filter-features --p-min-samples

# FILTER
if config["to_filter"] == True:
    print(f"Filtering table and/or sequences.\n\n")
    # Apply optional post-denoising length, abundance, and prevalence filters.
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
            minfreq=config["repseq_min_frequency"],
            minsamps=config["repseq_min_samples"]
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
            "--p-min-frequency {params.minfreq} "
            "--p-min-samples {params.minsamps} "
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

# Summarize feature table statistics and optional sample metadata.
rule summarize_feature_table:
    input:
        table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
        metadata=output_dir+config["run_name"]+"-repseqs/stats/metadata_used.tsv"
    output:
        output_dir+config["run_name"]+"-repseqs/stats/table_summary.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    threads: config["asv_threads"]
    shell:
        """
        if [ {use_metadata} == "yes" ]; then
            qiime feature-table summarize \
            --i-table {input.table} \
            --m-sample-metadata-file {input.metadata} \
            --o-visualization {output}
        else
            qiime feature-table summarize \
            --i-table {input.table} \
            --o-visualization {output}
        fi
        """

# Summarize denoising stats (dada2 or deblur)
if config["asv_method"] in ["dada2pe", "dada2se"]:
    stats_tsv = output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.tsv"
    # Visualize DADA2 denoising stats as a QIIME2 tabulation.
    rule summarize_dada2_repseqs:
        input:
            stats=output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qza"
        output:
            output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qzv"
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "qiime metadata tabulate --m-input-file {input.stats} --o-visualization {output}"
    
    # Export DADA2 stats visualization to a TSV for downstream inspection.
    rule export_dada2_summary_to_tsv:
        input:
            output_dir+config["run_name"]+"-repseqs/stats/dada2_stats.qzv"
        output:
            stats_tsv
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "unzip -qq -o {input} -d temp0; "
            "mv temp0/*/data/metadata.tsv {output}; "
            "/bin/rm -r temp0"
    
elif config["asv_method"] == "deblur":
    stats_tsv = output_dir+config["run_name"]+"-repseqs/stats/deblur_stats.tsv"
    # Visualize Deblur stats in QIIME2.
    rule summarize_deblur_repseqs:
        input:
            stats=output_dir+config["run_name"]+"-repseqs/stats/deblur_stats.qza"
        output:
            output_dir+config["run_name"]+"-repseqs/stats/deblur_stats.qzv"
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "qiime deblur visualize-stats --i-deblur-stats {input.stats} --o-visualization {output}"
    
    # Export Deblur statistics to TSV format.
    rule export_deblur_summary_to_tsv:
        input:
            output_dir+config["run_name"]+"-repseqs/stats/deblur_stats.qza"
        output:
            stats_tsv
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "unzip -qq -o {input} -d temp0; "
            "sed 's/,/\t/g' temp0/*/data/stats.csv > {output}; "
            "/bin/rm -r temp0"
    
else:
    raise ValueError("Invalid ASV method specified for stats export")


# Convert feature table to BIOM format for downstream utilities.
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

# Summarize sample depths from BIOM table.
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

# Summarize feature counts from BIOM table.
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

# Render representative sequences to an interactive QIIME2 visualization.
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

# Export representative sequences to FASTA for external analysis.
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

# Calculate sequence length distribution for representative sequences.
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

# Summarize representative sequence length statistics.
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

# Export feature table to TSV with Tourmaline-friendly header.
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
        "-o TEMP.tsv "
        "--to-tsv "
        "&& cat TEMP.tsv | tail -n +2 | sed 's/^#OTU ID/featureid/' > {output} "
        "&& /bin/rm TEMP.tsv"

# Generate alpha-rarefaction curves (optional).
rule diversity_alpha_rarefaction:
    input:
        table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
        metadata=output_dir+config["run_name"]+"-repseqs/stats/metadata_used.tsv"
    params:
        maxdepth=config["alpha_max_depth"],
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
            --m-metadata-file {input.metadata} \
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
# Run QIIME2 core-metrics workflow to produce rarefied tables and ordinations (optional).
rule diversity_core_metrics:
    input:
        table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza",
        metadata=output_dir+config["run_name"]+"-repseqs/stats/metadata_used.tsv"
    params:
        samplingdepth=config["core_sampling_depth"],
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
        "--m-metadata-file {input.metadata} "
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
