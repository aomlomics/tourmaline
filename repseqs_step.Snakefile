rule run_dada2_pe_denoise:
    """Run paired end dada2"""
    input:
        onfig["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza"
        #config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.fasta",
        #config["run_name"]+"-repseqs/"+config["run_name"]+"-table.biom",
        #config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs-stats.tsv",
        # add figures here

rule run_dada2_se_denoise:
    """Trim all reads with all supplied trimming parameters"""
    input:
        expand(config["run_name"]+"-samples/trimmed/fastqc_R1/{sample}.1_fastqc.html",sample=SAMPLES),
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest",
        config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt"

rule run_deblur_se_denoise:
    """Trim all reads with all supplied trimming parameters"""
    input:
        expand(config["run_name"]+"-samples/trimmed/fastqc_R1/{sample}.1_fastqc.html",sample=SAMPLES),
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest",
        config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt"


rule denoise_dada2_pe:
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza"
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
        hashedfeatureids=config["dada2pe_hashed_feature_ids"]
    output:
        table=config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza"
        repseqs=config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.qza",
        stats=config["run_name"]+"-repseqs/stats/dada2_stats.qza"
    conda:
        "qiime2-2023.5"
    threads: config["dada2pe_threads"]
    shell:
        "qiime dada2 denoise-paired "
        "--i-demultiplexed-seqs {input[0]} "
        "--p-trunc-len-f {params.trunclenf} "
        "--p-trunc-len-r {params.trunclenr} "
        "--p-trim-left-f {params.trimleftf} "
        "--p-trim-left-r {params.trimleftr} "
        "--p-max-ee-f {params.maxeef} "
        "--p-max-ee-r {params.maxeer} "
        "--p-trunc-q {params.truncq} "
        "--p-pooling-method {params.poolingmethod} "
        "--p-chimera-method {params.chimeramethod} "
        "--p-min-fold-parent-over-abundance {params.minfoldparentoverabundance} "
        "--p-n-reads-learn {params.nreadslearn} "
        "--p-n-threads {threads} "
        "{params.hashedfeatureids} "
        "--o-table {output.table} "
        "--o-representative-sequences {output.repseqs} "
        "--o-denoising-stats {output.stats} "
        "--verbose"