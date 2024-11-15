## changed import toPairedEndFastqManifestPhred33 instead of V2, for legacy manifests
## NEED TO CHANGE BACK

## Snakefile for sample step (paired-end reads) of Tourmalne V2 pipeline
from Bio.Seq import Seq
import glob
import os
#import shutil
import subprocess

output_dir = config["output_dir"]+"/"

if config["sample_manifest_file"] != None:
    if config["to_trim"] == True:
        print("Manifest provided, trimming reads\n")
        rule raw_fastq_demux_pe:
            input:
                config["sample_manifest_file"],
            output:
                output_dir+config["run_name"]+"-samples/raw_pe_fastq.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                qiime tools import \
                --type 'SampleData[PairedEndSequencesWithQuality]' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format PairedEndFastqManifestPhred33      
                """
        ruleorder: cutadapt_pe > import_fastq_demux_pe
        SAMPLES=['nothing']
    else:
        print("Manifest provided, not trimming reads\n")
        rule notrim_fastq_demux_pe:
            input:
                config["sample_manifest_file"],
            output:
                output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                qiime tools import \
                --type 'SampleData[PairedEndSequencesWithQuality]' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format PairedEndFastqManifestPhred33      
                """
        #ruleorder: import_fastq_demux_pe > cutadapt_pe
        SAMPLES=['nothing']

else:
    if config["raw_fastq_path"] != None:
        print("no manifest, trimming reads\n")
        SAMPLES, SFX = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1{suffix}.fastq.gz")
        SFX = list(set(SFX))
        rule make_raw_manifest_pe_path:
            input:
                fwdreads=expand(config["raw_fastq_path"]+"/{sample}_R1{suffix}.fastq.gz",sample=SAMPLES,suffix=SFX),
                revreads=expand(config["raw_fastq_path"]+"/{sample}_R2{suffix}.fastq.gz",sample=SAMPLES,suffix=SFX),
            output:
                output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_raw_pe.manifest"
            params:
                SFX=SFX
            shell:
                """
                echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > {output}
                if [ {params.SFX} == "_001" ]; then
                    for f in {input.fwdreads}
                        do
                            echo -e "$(basename $f | sed 's/_R1.*\.fastq.gz//')\t$f\t$(echo $f | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/g')" >> {output}
                        done
                else
                    for f in {input.fwdreads}
                        do
                            echo -e "$(basename $f | sed 's/_R1.*\.fastq.gz//')\t$f\t$(echo $f | sed 's/_R1.fastq.gz/_R2.fastq.gz/g')" >> {output}
                        done
                fi
                
            
                """
        rule import_raw_fastq_demux_pe:
            input:
                output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_raw_pe.manifest",
            output:
                output_dir+config["run_name"]+"-samples/raw_pe_fastq.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                qiime tools import \
                --type 'SampleData[PairedEndSequencesWithQuality]' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format PairedEndFastqManifestPhred33V2      
                """
        ruleorder: cutadapt_pe > import_fastq_demux_pe
    elif config["trimmed_fastq_path"] != None:
        print("no manifest, not trimming reads\n")
        SAMPLES, =glob_wildcards(config["trimmed_fastq_path"]+"/{sample}.1.fastq.gz")
        ruleorder: import_fastq_demux_pe > cutadapt_pe
        ruleorder: make_manifest_pe > make_manifest_pe_file




primerF=config["fwd_primer"]
primerR=config["rev_primer"]
revcomp_primerF=Seq(config["fwd_primer"]).reverse_complement()
revcomp_primerR=Seq(config["rev_primer"]).reverse_complement()

rule trim_pe_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        output_dir+config["run_name"]+"-samples/raw_pe_fastq.qza",
        output_dir+config["run_name"]+"-samples/stats/raw_fastq_summary.qzv",
        output_dir+config["run_name"]+"-samples/stats/fastq_summary.qzv",
        #config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule no_trim_pe_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        output_dir+config["run_name"]+"-samples/stats/fastq_summary.qzv",
        #config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",



rule raw_fastq_summary:
    input:
        output_dir+config["run_name"]+"-samples/raw_pe_fastq.qza"
    output:
        output_dir+config["run_name"]+"-samples/stats/raw_fastq_summary.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"

rule cutadapt_pe:
    input:
        #rules.import_raw_fastq_demux_pe.output
        output_dir+config["run_name"]+"-samples/raw_pe_fastq.qza"
    output:
        output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza",
        output_dir+config["run_name"]+"-samples/stats/cutadapt_summary.txt",
    params:
        discard=config['discard_untrimmed'],
        min_len=config['minimum_length']
    conda:
        "qiime2-amplicon-2024.10"
    threads: workflow.cores
    shell:
        """
        # Gotta figure this out$( (( "{params.discard}" == "True" )) && printf %s '--p-discard-untrimmed' ) 
        # https://stackoverflow.com/questions/28678505/add-command-arguments-using-inline-if-statement-in-bash
        qiime cutadapt trim-paired \
        --p-cores {threads} \
        --i-demultiplexed-sequences {input} \
        --p-adapter-f {revcomp_primerR} \
        --p-adapter-r {revcomp_primerF} \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-minimum-length {params.min_len} \
        --verbose \
        --o-trimmed-sequences trimmed_1.qza 1> {output[1]}

        if [ {params.discard} = True ]; then
            qiime cutadapt trim-paired \
            --p-cores {threads} \
            --i-demultiplexed-sequences trimmed_1.qza \
            --p-front-f {primerF} \
            --p-front-r {primerR} \
            --p-match-read-wildcards \
            --p-match-adapter-wildcards \
            --p-discard-untrimmed \
            --p-minimum-length {params.min_len} \
            --verbose \
            --o-trimmed-sequences {output[0]} 1>> {output[1]}
        else
            qiime cutadapt trim-paired \
            --p-cores {threads} \
            --i-demultiplexed-sequences trimmed_1.qza \
            --p-front-f {primerF} \
            --p-front-r {primerR} \
            --p-match-read-wildcards \
            --p-match-adapter-wildcards \
            --p-minimum-length {params.min_len} \
            --verbose \
            --o-trimmed-sequences {output[0]} 1>> {output[1]}
        fi;

        rm trimmed_1.qza
        """

# make updated manifest file

rule make_manifest_pe:
    input:
        fwdreads=expand(output_dir+config["run_name"]+"-samples/trimmed/{sample}.1.fastq",sample=SAMPLES),
        revreads=expand(output_dir+config["run_name"]+"-samples/trimmed/{sample}.2.fastq",sample=SAMPLES),
    output:
        output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest"
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
        output_dir+"create_manifest.txt"
    output:
        output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest"
    params:
        file=config["sample_manifest_file"]
    run:
        if config["sample_manifest_file"] != None:
            os.makedirs(output_dir+config["run_name"]+"-samples", exist_ok=True)
            # copy the file
            shutil.copy(params.file, output[0])

rule import_fastq_demux_pe:
    input:
        output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_pe.manifest"
    output:
        output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input} "
        "--output-path {output} "
        "--input-format PairedEndFastqManifestPhred33V2"

rule summarize_fastq_demux_pe:
    input:
        output_dir+config["run_name"]+"-samples/"+config["run_name"]+"_fastq_pe.qza"
    output:
        output_dir+config["run_name"]+"-samples/stats/fastq_summary.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"


# make stats file with script

#not done
#rule trim_summary_stats:
    # input:
    #     config["run_name"]+"-samples/stats/fastq_summary.qzv",
    #     config["run_name"]+"-samples/stats/raw_fastq_summary.qzv"
    # output:
    #     config["run_name"]+"-samples/stats/"
    # conda:
    #     "qiime2-amplicon-2024.10"
    # threads: config["other_threads"]
    # shell:
    #     "unzip -qq -o {input} -d temp0; "
    #     "mv temp0/*/data/per-sample-fastq-counts.tsv {output}; "
    #     "/bin/rm -r temp0"


rule check_seq_qual_dropoff:
    input:
        R1=output_dir+config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        R2=output_dir+config["run_name"]+"-samples/trimmed/fastqc_R2/multiqc_report_R2.html",
    output:
        output_dir+config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",
    params:
        seq_cutoff = config["seq_quality_cutoff"],
        seq_qual_R1 = output_dir+config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt",
        seq_qual_R2 = output_dir+config["run_name"]+"-samples/trimmed/fastqc_R2/multiqc_report_R2_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt",
        stats_dir = output_dir+config["run_name"]+"-samples/stats/",
    conda:
        "qiime2-amplicon-2024.10"
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


