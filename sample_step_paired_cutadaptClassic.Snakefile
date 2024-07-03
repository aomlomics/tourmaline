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

rule cutadapt_pe_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        expand(config["run_name"]+"-samples/trimmed/fastqc_R1/{sample}.1_fastqc.html",sample=SAMPLES),
        expand(config["run_name"]+"-samples/trimmed/fastqc_R2/{sample}.2_fastqc.html",sample=SAMPLES),
        config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest",
        config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        config["run_name"]+"-samples/trimmed/fastqc_R2/multiqc_report_R2.html",
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",
        config["run_name"]+"-samples/stats/fastq_summary.qzv"



rule cutadapt_pe:
    input:
        fastqF=config["raw_fastq_path"]+"/{sample}_R1.fastq.gz",
        fastqR=config["raw_fastq_path"]+"/{sample}_R2.fastq.gz"
    output:
        fastq1=config["run_name"]+"-samples/trimmed/{sample}.1.fastq",
        fastq2=config["run_name"]+"-samples/trimmed/{sample}.2.fastq",
        qc=config["run_name"]+"-samples/trimmed/{sample}_paired.qc.txt",
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
        -A '{A_primer}' \
        {params.extra} \
        -o {output.fastq1} -p {output.fastq2} \
        {input.fastqF} {input.fastqR} > {output.qc}
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
        R2=config["run_name"]+"-samples/trimmed/fastqc_R2/multiqc_report_R2.html",
        #zip=config["run_name"]+"-trimmed/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        R1dir = config["run_name"]+"-samples/trimmed/fastqc_R1",
        R2dir = config["run_name"]+"-samples/trimmed/fastqc_R2",
    conda:
        "multiqc"
    threads: 1
    shell:
        """
        multiqc {params.R1dir} --export -n multiqc_report_R1 -o {params.R1dir};
        multiqc {params.R2dir} --export -n multiqc_report_R2 -o {params.R2dir};
        """

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


