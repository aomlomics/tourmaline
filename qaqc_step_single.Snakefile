## Snakefile for qaqc step (single-end reads) of Tourmalne V2 pipeline

## NEED TO add manifest rules as in the paired end version

import os
import shutil

output_dir = config["output_dir"]+"/"

# HELPER FUNCTIONS

def check_file_separator(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline()
        if '\t' in first_line:
            return 'tab'
        elif ',' in first_line:
            return 'comma'
        else:
            return 'unknown'

def check_files_for_substring(directory, substring):
    for file_name in os.listdir(directory):
        if substring in file_name:
            return True
    return False


if config["sample_manifest_file"] != None:
    separator = check_file_separator(config["sample_manifest_file"])
    if config["to_trim"] == True:
        print("Manifest provided, trimming reads\n")
        rule raw_fastq_demux_pe:
            input:
                config["sample_manifest_file"],
            output:
                output_dir+config["run_name"]+"-qaqc/raw_se_fastq.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                if [ {separator} == "tab" ]; then
                    FORMAT='SingleEndFastqManifestPhred33V2'
                elif [ {separator} == "comma" ]; then
                    FORMAT='SingleEndFastqManifestPhred33'
                else
                    echo "The manifest file is neither tab-separated nor comma-separated."
                fi
                qiime tools import \
                --type 'SampleData[SequencesWithQuality]' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format $FORMAT     
                """
        ruleorder: cutadapt_se > import_fastq_demux_se
        SAMPLES=['nothing']
    else:
        print("Manifest provided, not trimming reads\n")
        rule notrim_fastq_demux_se:
            input:
                config["sample_manifest_file"],
            output:
                output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq_se.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                if [ {separator} == "tab" ]; then
                    FORMAT='SingleEndFastqManifestPhred33V2'
                elif [ {separator} == "comma" ]; then
                    FORMAT='SingleEndFastqManifestPhred33'
                else
                    echo "The manifest file is neither tab-separated nor comma-separated."
                fi
                qiime tools import \
                --type 'SampleData[SequencesWithQuality]' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format $FORMAT      
                """
        SAMPLES=['nothing']
else:
    if config["raw_fastq_path"] != None:
        print("no manifest, trimming reads\n")
        if check_files_for_substring(config["raw_fastq_path"], "_001.fastq"):
            SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1_001.fastq.gz")
            fwdreads=expand(config["raw_fastq_path"]+"/{sample}_R1_001.fastq.gz",sample=SAMPLES)
        else:
            SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz")
            fwdreads=expand(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz",sample=SAMPLES)
        rule make_raw_manifest_se_path:
            input:
                fwdreads
            output:
                output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_raw_se.manifest"
            shell:
                """
                echo -e "sample-id\tabsolute-filepath" > {output}
                for f in {input}
                    do
                        echo "$(basename $f | sed 's/_R1.*\.fastq.gz//')\t$f" >> {output}
                    done
            
                """
        
        rule import_raw_fastq_demux_se:
            input:
                output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_raw_se.manifest",
            output:
                output_dir+config["run_name"]+"-qaqc/raw_se_fastq.qza"
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                "qiime tools import "
                "--type 'SampleData[SequencesWithQuality]' "
                "--input-path {input[0]} "
                "--output-path {output} "
                "--input-format SingleEndFastqManifestPhred33V2"
    
        ruleorder: cutadapt_se > import_fastq_demux_se
    elif config["trimmed_fastq_path"] != None:
        print("no manifest, not trimming reads\n")
        SAMPLES, =glob_wildcards(config["trimmed_fastq_path"]+"/{sample}.1.fastq.gz")
        ruleorder: import_fastq_demux_se > cutadapt_se
        ruleorder: make_manifest_se > make_manifest_se_file


#ruleorder: make_raw_manifest_se_path > make_raw_manifest_se_file


rule trim_se_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        #output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_raw_se.manifest",
        output_dir+config["run_name"]+"-qaqc/raw_se_fastq.qza",
        output_dir+config["run_name"]+"-qaqc/stats/fastq_summary.qzv",
        output_dir+config["run_name"]+"-qaqc/stats/raw_fastq_summary.qzv",
        #config["run_name"]+"-qaqc/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule no_trim_se_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        #output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_se.manifest",
        output_dir+config["run_name"]+"-qaqc/stats/fastq_summary.qzv",
        #config["run_name"]+"-qaqc/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",



rule raw_fastq_summary:
    input:
        output_dir+config["run_name"]+"-qaqc/raw_se_fastq.qza"
    output:
        output_dir+config["run_name"]+"-qaqc/stats/raw_fastq_summary.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"

rule cutadapt_se:
    input:
        #rules.import_raw_fastq_demux_se.output
        output_dir+config["run_name"]+"-qaqc/raw_se_fastq.qza"
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq_se.qza",
        output_dir+config["run_name"]+"-qaqc/stats/cutadapt_summary.txt",
    params:
        discard=config['discard_untrimmed'],
        min_len=config['minimum_length'],
        primerF=config["fwd_primer"],
        primerR=config["rev_primer"],
        tempTrim=output_dir+config["run_name"]+"-qaqc/trimmed_1.qza"
    conda:
        "qiime2-amplicon-2024.10"
    threads: workflow.cores
    shell:
        """
        #revcomp_primerF=`echo {params.primerF} | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
        revcomp_primerR=`echo {params.primerR} | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`

        qiime cutadapt trim-single \
        --p-cores {threads} \
        --i-demultiplexed-sequences {input} \
        --p-adapter $revcomp_primerR \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --p-minimum-length {params.min_len} \
        --verbose \
        --o-trimmed-sequences {params.tempTrim} 1> {output[1]}

        if [ {params.discard} = True ]; then
            qiime cutadapt trim-single \
            --p-cores {threads} \
            --i-demultiplexed-sequences {params.tempTrim} \
            --p-front {params.primerF} \
            --p-match-read-wildcards \
            --p-match-adapter-wildcards \
            --p-discard-untrimmed \
            --p-minimum-length {params.min_len} \
            --verbose \
            --o-trimmed-sequences {output[0]} 1>> {output[1]}

        else
            qiime cutadapt trim-single \
            --p-cores {threads} \
            --i-demultiplexed-sequences {params.tempTrim} \
            --p-front {params.primerF} \
            --p-match-read-wildcards \
            --p-match-adapter-wildcards \
            --p-minimum-length {params.min_len} \
            --verbose \
            --o-trimmed-sequences {output[0]} 1>> {output[1]}
        fi;        
        rm {params.tempTrim}
        """

# make updated manifest file

rule make_manifest_se:
    input:
        fwdreads=expand(output_dir+config["run_name"]+"-qaqc/trimmed/{sample}.1.fastq",sample=SAMPLES),
        #lambda wildcards: os.path.join(output_dir+config["run_name"]+"-qaqc/trimmed/", f"{wildcards.sample}.1.fastq")
        #os.path.join(output_dir+config["run_name"]+"-qaqc/trimmed/{sample}.1.fastq")
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_se.manifest"
    shell:
        """
        echo -e "sample-id\tabsolute-filepath" > {output}
        for f in {input}
        do
            echo "$(basename $f | sed 's/.1.fastq//')\t$f" >> {output}
        done
       
        """

rule make_manifest_se_file:
    input:
        output_dir+"create_manifest.txt"
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_se.manifest"
    params:
        file=config["sample_manifest_file"]
    run:
        if config["sample_manifest_file"] != None:
            os.makedirs(config["run_name"]+"-qaqc", exist_ok=True)
            # copy the file
            shutil.copy(params.file, output[0])
        os.remove(output_dir+"create_manifest.txt")

rule import_fastq_demux_se:
    input:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_se.manifest"
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq_se.qza"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime tools import "
        "--type 'SampleData[SequencesWithQuality]' "
        "--input-path {input[0]} "
        "--output-path {output} "
        "--input-format SingleEndFastqManifestPhred33V2"

rule summarize_fastq_demux_se:
    input:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq_se.qza"
    output:
        output_dir+config["run_name"]+"-qaqc/stats/fastq_summary.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"


# make stats file with script

rule check_seq_qual_dropoff:
    input:
        R1=output_dir+config["run_name"]+"-qaqc/trimmed/fastqc_R1/multiqc_report_R1.html",
    output:
        output_dir+config["run_name"]+"-qaqc/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",
    params:
        seq_cutoff = config["seq_quality_cutoff"],
        seq_qual_R1 = output_dir+config["run_name"]+"-qaqc/trimmed/fastqc_R1/multiqc_report_R1_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt",
        stats_dir = output_dir+config["run_name"]+"-qaqc/stats/",
    conda:
        "qiime2-amplicon-2024.10"
    threads: 1
    shell:
        """
        mkdir -p {params.stats_dir}
        echo "R1 seq quality" > {output};
        python scripts/fastqc_per_base_sequence_quality_dropoff.py \
        --input {params.seq_qual_R1} --cutoff {params.seq_cutoff} >> {output};
        """









