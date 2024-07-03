from Bio.Seq import Seq
import os
import shutil

if config["sample_manifest_file"] != None:
    if config["to_trim"] == True:
        print("Manifest provided, trimming reads\n")
        ruleorder: cutadapt_pe > import_fastq_demux_pe
        SAMPLES=['nothing']
    else:
        print("Manifest provided, not trimming reads\n")
        ruleorder: import_fastq_demux_pe > cutadapt_pe
        SAMPLES=['nothing']
else:
    with open("create_manifest.txt", 'w') as file:
        pass
    if config["raw_fastq_path"] != None:
        print("no manifest, trimming reads\n")
        SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz")
        cutadapt_pe > import_fastq_demux_pe
    elif config["trimmed_fastq_path"] != None:
        print("no manifest, not trimming reads\n")
        SAMPLES, =glob_wildcards(config["trimmed_fastq_path"]+"/{sample}.1.fastq.gz")
        import_fastq_demux_pe > cutadapt_pe


ruleorder: make_raw_manifest_pe_path > make_raw_manifest_pe_file
ruleorder: make_manifest_pe > make_manifest_pe_file

#         make manifest_pe.tsv
#         import to config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza"
#         demux summarize to stats/fastq_quality_summary.qzv
#         extract forward-seven-number-summaries.tsv from qzv, run code on it


primerF=config["fwd_primer"]
primerR=config["rev_primer"]
revcomp_primerF=Seq(config["fwd_primer"]).reverse_complement()
revcomp_primerR=Seq(config["rev_primer"]).reverse_complement()

rule trim_pe_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_pe.manifest",
        config["run_name"]+"-samples/raw_pe_fastq.qza",
        config["run_name"]+"-samples/stats/fastq_summary.qzv",
        #config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule no_trim_pe_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest",
        config["run_name"]+"-samples/stats/fastq_summary.qzv",
        #config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule make_raw_manifest_pe_path:
    input:
        fwdreads=expand(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz",sample=SAMPLES),
        revreads=expand(config["raw_fastq_path"]+"/{sample}_R2.fastq.gz",sample=SAMPLES),
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_pe.manifest"
    shell:
        """
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > {output}
        for f in {input.fwdreads}
        do
            echo -e "$(basename $f | sed 's/_R1.fastq.gz//')\t$f\t$(echo $f | sed 's/_R1.fastq.gz/_R2.fastq.gz/g')" >> {output}
        done
       
        """

rule make_raw_manifest_pe_file:
    input:
        config["sample_manifest_file"]
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_pe.manifest"
    params:
        outdir=config["run_name"]+"-samples"
    run:
        if config["sample_manifest_file"] != None:
            os.makedirs(config["run_name"]+"-samples", exist_ok=True)
            # copy the file
            shutil.copy(input[0], output[0])

rule import_raw_fastq_demux_pe:
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_raw_pe.manifest",
    output:
        config["run_name"]+"-samples/raw_pe_fastq.qza"
    conda:
        "qiime2-2023.5"
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input[0]} "
        "--output-path {output} "
        "--input-format PairedEndFastqManifestPhred33V2"


rule cutadapt_pe:
    input:
        #rules.import_raw_fastq_demux_pe.output
        config["run_name"]+"-samples/raw_pe_fastq.qza"
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza",
        config["run_name"]+"-samples/stats/cutadapt_summary.txt",
    conda:
        "qiime2-2023.5"
    threads: workflow.cores
    shell:
        """
        qiime cutadapt trim-paired \
        --p-cores {threads} \
        --i-demultiplexed-sequences {input} \
        --p-adapter-f {revcomp_primerR} \
        --p-adapter-r {revcomp_primerF} \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-minimum-length 50 \
        --verbose \
        --o-trimmed-sequences trimmed_1.qza 1> {output[1]}

        qiime cutadapt trim-paired \
        --p-cores {threads} \
        --i-demultiplexed-sequences trimmed_1.qza \
        --p-front-f {primerF} \
        --p-front-r {primerR} \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-discard-untrimmed \
        --p-minimum-length 50 \
        --verbose \
        --o-trimmed-sequences {output[0]} 1>> {output[1]}

        rm trimmed_1.qza
        """

# make updated manifest file

rule make_manifest_pe:
    input:
        fwdreads=expand(config["run_name"]+"-samples/trimmed/{sample}.1.fastq",sample=SAMPLES),
        revreads=expand(config["run_name"]+"-samples/trimmed/{sample}.2.fastq",sample=SAMPLES),
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest"
    shell:
        """
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > {output}
        for f in {input.fwdreads}
        do
            echo -e "$(basename $f | sed 's/.1.fastq//')\t$(pwd)/$f\t$(pwd)/$(echo $f | sed 's/.1.fastq/.2.fastq/g')" >> {output}
        done
       
        """

rule make_manifest_pe_file:
    input:
        config["sample_manifest_file"]
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest"
    run:
        if config["sample_manifest_file"] != None:
            os.makedirs(config["run_name"]+"-samples", exist_ok=True)
            # copy the file
            shutil.copy(input[0], output[0])

rule import_fastq_demux_pe:
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest"
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza"
    conda:
        "qiime2-2023.5"
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input} "
        "--output-path {output} "
        "--input-format PairedEndFastqManifestPhred33V2"

rule summarize_fastq_demux_pe:
    input:
        config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza"
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
        R2=config["run_name"]+"-samples/trimmed/fastqc_R2/multiqc_report_R2.html",
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


