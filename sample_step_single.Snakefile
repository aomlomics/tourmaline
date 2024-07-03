from Bio.Seq import Seq
import os
import shutil

if config["sample_manifest_file"] != None:
    if config["to_trim"] == True:
        print("Manifest provided, trimming reads\n")
        ruleorder: cutadapt_se > import_fastq_demux_se
        SAMPLES=['nothing']
    else:
        print("Manifest provided, not trimming reads\n")
        ruleorder: import_fastq_demux_se > cutadapt_se
        SAMPLES=['nothing']
else:
    with open("create_manifest.txt", 'w') as file:
        pass
    if config["raw_fastq_path"] != None:
        print("no manifest, trimming reads\n")
        SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz")
        ruleorder: cutadapt_se > import_fastq_demux_se
    elif config["trimmed_fastq_path"] != None:
        print("no manifest, not trimming reads\n")
        SAMPLES, =glob_wildcards(config["trimmed_fastq_path"]+"/{sample}.1.fastq.gz")
        ruleorder: import_fastq_demux_se > cutadapt_se


ruleorder: make_raw_manifest_se_path > make_raw_manifest_se_file
ruleorder: make_manifest_se > make_manifest_se_file


primerF=config["fwd_primer"]
primerR=config["rev_primer"]
revcomp_primerF=Seq(config["fwd_primer"]).reverse_complement()
revcomp_primerR=Seq(config["rev_primer"]).reverse_complement()

rule trim_se_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_se.manifest",
        config["run_name"]+"-samples/raw_se_fastq.qza",
        config["run_name"]+"-samples/stats/fastq_summary.qzv",
        #config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule no_trim_se_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest",
        config["run_name"]+"-samples/stats/fastq_summary.qzv",
        #config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule make_raw_manifest_se_path:
    input:
        fwdreads=expand(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz",sample=SAMPLES),
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_se.manifest"
    shell:
        """
        echo -e "sample-id\tabsolute-filepath" > {output}
        for f in {input.fwdreads}
        do
            echo "$(basename $f | sed 's/_R1.fastq.gz//')\t$f" >> {output}
        done
       
        """

rule make_raw_manifest_se_file:
    input:
        "create_manifest.txt"
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_se.manifest"
    params:
        file=config["sample_manifest_file"]
    run:
        if config["sample_manifest_file"] != None:
            os.makedirs(config["run_name"]+"-samples", exist_ok=True)
            # copy the file
            shutil.copy(params.file, output[0])
        os.remove("create_manifest.txt")

rule import_raw_fastq_demux_se:
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_se.manifest",
    output:
        config["run_name"]+"-samples/raw_se_fastq.qza"
    conda:
        "qiime2-2023.5"
    shell:
        "qiime tools import "
        "--type 'SampleData[SequencesWithQuality]' "
        "--input-path {input[0]} "
        "--output-path {output} "
        "--input-format SingleEndFastqManifestPhred33V2"


rule cutadapt_se:
    input:
        #rules.import_raw_fastq_demux_se.output
        config["run_name"]+"-samples/raw_se_fastq.qza"
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_fastq_se.qza",
        config["run_name"]+"-samples/stats/cutadapt_summary.txt",
    conda:
        "qiime2-2023.5"
    threads: workflow.cores
    shell:
        """
        qiime cutadapt trim-single \
        --p-cores {threads} \
        --i-demultiplexed-sequences {input} \
        --p-adapter {revcomp_primerR} \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-minimum-length 50 \
        --verbose \
        --o-trimmed-sequences trimmed_1.qza 1> {output[1]}

        qiime cutadapt trim-single \
        --p-cores {threads} \
        --i-demultiplexed-sequences trimmed_1.qza \
        --p-front {primerF} \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-discard-untrimmed \
        --p-minimum-length 50 \
        --verbose \
        --o-trimmed-sequences {output[0]} 1>> {output[1]}

        rm trimmed_1.qza
        """

# make updated manifest file

rule make_manifest_se:
    input:
        fwdreads=expand(config["run_name"]+"-samples/trimmed/{sample}.1.fastq",sample=SAMPLES),
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest"
    shell:
        """
        echo -e "sample-id\tabsolute-filepath" > {output}
        for f in {input.fwdreads}
        do
            echo "$(basename $f | sed 's/.1.fastq//')\t$f" >> {output}
        done
       
        """

rule make_manifest_se_file:
    input:
        "create_manifest.txt"
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest"
    params:
        file=config["sample_manifest_file"]
    run:
        if config["sample_manifest_file"] != None:
            os.makedirs(config["run_name"]+"-samples", exist_ok=True)
            # copy the file
            shutil.copy(params.file, output[0])
        os.remove("create_manifest.txt")

rule import_fastq_demux_se:
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest"
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_fastq_se.qza"
    conda:
        "qiime2-2023.5"
    shell:
        "qiime tools import "
        "--type 'SampleData[SequencesWithQuality]' "
        "--input-path {input[0]} "
        "--output-path {output} "
        "--input-format SingleEndFastqManifestPhred33V2"

rule summarize_fastq_demux_se:
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_fastq_se.qza"
    output:
        config["run_name"]+"-samples/stats/fastq_summary.qzv"
    conda:
        "qiime2-2023.5"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"


# make stats file with script

rule check_seq_qual_dropoff:
    input:
        R1=config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
    output:
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",
    params:
        seq_cutoff = config["seq_quality_cutoff"],
        seq_qual_R1 = config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt",
        seq_qual_R2 = config["run_name"]+"-samples/trimmed/fastqc_R2/multiqc_report_R2_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt",
        stats_dir = config["run_name"]+"-samples/stats/",
    conda:
        "qiime2-2023.5"
    threads: 1
    shell:
        """
        mkdir -p {params.stats_dir}
        echo "R1 seq quality" > {output};
        python scripts/fastqc_per_base_sequence_quality_dropoff.py \
        --input {params.seq_qual_R1} --cutoff {params.seq_cutoff} >> {output};
        echo "R2 seq quality" >> {output}
        python scripts/fastqc_per_base_sequence_quality_dropoff.py \
        --input {params.seq_qual_R2} --cutoff {params.seq_cutoff} >> {output};
        """









