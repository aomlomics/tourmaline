from Bio.Seq import Seq
import os

if config["raw_fastq_path"] != "":
    SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz")
elif config["trimmed_fastq_path"] != "":
    SAMPLES, =glob_wildcards(config["trimmed_fastq_path"]+"/{sample}.1.fastq.gz")
else:
    print("ERROR: Need to provide either a path to raw fastq files or trimmed fastq files")

primerF=config["fwd_primer"]
primerR=config["rev_primer"]
revcomp_primerF=Seq(config["fwd_primer"]).reverse_complement()
revcomp_primerR=Seq(config["rev_primer"]).reverse_complement()

a_primer = f"{primerF};required...{revcomp_primerR};optional"
A_primer = f"{primerR};required...{revcomp_primerF};optional"

rule cutadapt_se_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        expand(config["run_name"]+"-samples/trimmed/fastqc_R1/{sample}.1_fastqc.html",sample=SAMPLES),
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest",
        config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule cutadapt_se:
    input:
        fastqF=config["raw_fastq_path"]+"/{sample}_R1.fastq.gz"
    output:
        fastq1=config["run_name"]+"-samples/trimmed/{sample}.1.fastq",
        qc=config["run_name"]+"-samples/trimmed/{sample}_single.qc.txt",
    params:
        run_name=config["run_name"],
        extra=config["trim_params"],
    conda:
        "qiime2-2023.5"
    threads: config["trimming_threads"]
    shell:
        """
        mkdir -p {params.run_name}-samples;
        mkdir -p {params.run_name}-samples/trimmed/;
        cutadapt \
        -a '{a_primer}' \
        {params.extra} \
        -o {output.fastq1} \
        {input.fastqF} > {output.qc}
        """

# make updated manifest file

rule make_manifest_se:
    input:
        fwdreads=expand(config["run_name"]+"-samples/trimmed/{sample}.1.fastq",sample=SAMPLES),
    output:
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest"
    shell:
        """
        echo "sample-id,forward-absolute-filepath" > {output}
        for f in {input.fwdreads}
        do
            echo "$(basename $f | sed 's/.1.fastq//'),$(pwd)/$f" >> {output}
        done
        """

# make stats file with multiqc

# make stats file with multiqc

rule fastqc:
    input:
        config["run_name"]+"-samples/trimmed/{sample}.{read}.fastq"
    output:
        html=config["run_name"]+"-samples/trimmed/fastqc_R{read}/{sample}.{read}_fastqc.html",
        #zip=config["run_name"]+"-trimmed/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet",
        outdir = config["run_name"]+"-samples/trimmed/fastqc_R{read}"
    conda:
        "multiqc"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir};
        fastqc --threads {threads} {input} --outdir {params.outdir}
        """

rule multiqc:
    input:
        expand(config["run_name"]+"-samples/trimmed/fastqc_R{read}/{sample}.{read}_fastqc.html",sample=SAMPLES,read=[1,2]),
    output:
        R1=config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        #zip=config["run_name"]+"-trimmed/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        R1dir = config["run_name"]+"-samples/trimmed/fastqc_R1",
    conda:
        "multiqc"
    threads: 1
    shell:
        """
        multiqc {params.R1dir} --export -n multiqc_report_R1 -o {params.R1dir};
        """

rule check_seq_qual_dropoff:
    input:
        R1=config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
    output:
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",
    params:
        seq_cutoff = config["seq_quality_cutoff"],
        seq_qual_R1 = config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt",
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
        """

