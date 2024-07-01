from Bio.Seq import Seq
import os

SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz")

primerF=config["fwd_primer"]
primerR=config["rev_primer"]
revcomp_primerF=Seq(config["fwd_primer"]).reverse_complement()
revcomp_primerR=Seq(config["rev_primer"]).reverse_complement()

a_primer = f"{primerF};required...{revcomp_primerR};optional"
A_primer = f"{primerR};required...{revcomp_primerF};optional"

rule cutadapt_se_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        #expand(config["run_name"]+"-trimmed/{sample}_single.qc.txt",sample=SAMPLES),
        config["run_name"]+"_se.manifest"

rule cutadapt_se:
    input:
        fastqF=config["raw_fastq_path"]+"/{sample}_R1.fastq.gz"
    output:
        fastq1=config["run_name"]+"-trimmed/{sample}.1.fastq",
        qc=config["run_name"]+"-trimmed/{sample}_single.qc.txt",
    params:
        run_name=config["run_name"],
        extra=config["trim_params"],
        paired=config["paired_end"],
    conda:
        "qiime2-2023.5"
    threads: config["trimming_threads"]
    shell:
        """
        mkdir -p {params.run_name}-trimmed/;
        cutadapt \
        -a '{a_primer}' \
        {params.extra} \
        -o {output.fastq1} \
        {input.fastqF} > {output.qc}
        """

# make updated manifest file

rule make_manifest_se:
    input:
        fwdreads=expand(config["run_name"]+"-trimmed/{sample}.1.fastq",sample=SAMPLES),
    output:
        config["run_name"]+"_se.manifest"
    shell:
        """
        echo "sample-id,forward-absolute-filepath" > {output}
        for f in {input.fwdreads}
        do
            echo "$(basename $f | sed 's/.1.fastq//'),$(pwd)/$f" >> {output}
        done
        """

# make stats file with multiqc

