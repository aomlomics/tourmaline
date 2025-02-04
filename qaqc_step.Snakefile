## Unified Snakefile for qaqc step of Tourmaline V2 pipeline
## Handles both single-end and paired-end reads

from Bio.Seq import Seq
import os
import shutil

output_dir = config["output_dir"]+"/"

# Copy config file to output directory
config_output_path = output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"-qaqc_config.yaml"
os.makedirs(os.path.dirname(config_output_path), exist_ok=True)
shutil.copy(workflow.configfiles[0], config_output_path)


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

# Set data type specific variables
IS_PAIRED = config.get("paired_end", False)
SEQUENCE_TYPE = "PairedEnd" if IS_PAIRED else "SingleEnd"
MANIFEST_FORMAT = f"{SEQUENCE_TYPE}FastqManifestPhred33V2"
IMPORT_TYPE = "SampleData[PairedEndSequencesWithQuality]" if IS_PAIRED else "SampleData[SequencesWithQuality]"
CUTADAPT_COMMAND = "trim-paired" if IS_PAIRED else "trim-single"

rule trim_all:
    """Trim all reads with all supplied trimming parameters"""
    input:
        output_dir+config["run_name"]+"-qaqc/raw_fastq.qza",
        output_dir+config["run_name"]+"-qaqc/stats/raw_fastq_summary.qzv",
        output_dir+config["run_name"]+"-qaqc/stats/fastq_summary.qzv"

rule no_trim_all:
    """Process all reads without trimming"""
    input:
        output_dir+config["run_name"]+"-qaqc/stats/fastq_summary.qzv",

if config["sample_manifest_file"] != None:
    separator = check_file_separator(config["sample_manifest_file"])
    if config["to_trim"] == True:
        print(f"Manifest provided, trimming {SEQUENCE_TYPE.lower()} reads\n")
        rule raw_fastq_demux:
            input:
                config["sample_manifest_file"],
            output:
                output_dir+config["run_name"]+"-qaqc/raw_fastq.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                if [ {separator} == "tab" ]; then
                    FORMAT='{SEQUENCE_TYPE}FastqManifestPhred33V2'
                elif [ {separator} == "comma" ]; then
                    FORMAT='{SEQUENCE_TYPE}FastqManifestPhred33'
                else
                    echo "The manifest file is neither tab-separated nor comma-separated."
                fi
                qiime tools import \
                --type '{IMPORT_TYPE}' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format $FORMAT     
                """
        ruleorder: cutadapt > import_fastq_demux
        SAMPLES=['nothing']
    else:
        print(f"Manifest provided, not trimming {SEQUENCE_TYPE.lower()} reads\n")
        rule notrim_fastq_demux:
            input:
                config["sample_manifest_file"],
            output:
                output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                if [ {separator} == "tab" ]; then
                    FORMAT='{SEQUENCE_TYPE}FastqManifestPhred33V2'
                elif [ {separator} == "comma" ]; then
                    FORMAT='{SEQUENCE_TYPE}FastqManifestPhred33'
                else
                    echo "The manifest file is neither tab-separated nor comma-separated."
                fi
                qiime tools import \
                --type '{IMPORT_TYPE}' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format $FORMAT      
                """
        SAMPLES=['nothing']
else:
    if config["raw_fastq_path"] != None:
        print(f"no manifest, trimming {SEQUENCE_TYPE.lower()} reads\n")
        if check_files_for_substring(config["raw_fastq_path"], "_001.fastq"):
            SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1_001.fastq.gz")
            fwdreads=expand(config["raw_fastq_path"]+"/{sample}_R1_001.fastq.gz",sample=SAMPLES)
            SUF="_001"
        else:
            SAMPLES, = glob_wildcards(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz")
            fwdreads=expand(config["raw_fastq_path"]+"/{sample}_R1.fastq.gz",sample=SAMPLES)
            SUF="other"
        
        if IS_PAIRED:
            rule make_raw_manifest_pe_path:
                input:
                    fwdreads
                output:
                    output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_raw.manifest"
                params:
                    SUF=SUF
                shell:
                    """
                    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > {output}
                    if [ {params.SUF} == "_001" ]; then
                        for f in {input}
                            do
                                echo -e "$(basename $f | sed 's/_R1.*\.fastq.gz//')\t$f\t$(echo $f | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/g')" >> {output}
                            done
                    else
                        for f in {input}
                            do
                                echo -e "$(basename $f | sed 's/_R1.*\.fastq.gz//')\t$f\t$(echo $f | sed 's/_R1.fastq.gz/_R2.fastq.gz/g')" >> {output}
                            done
                    fi
                    """
        else:
            rule make_raw_manifest_se_path:
                input:
                    fwdreads
                output:
                    output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_raw.manifest"
                shell:
                    """
                    echo -e "sample-id\tabsolute-filepath" > {output}
                    for f in {input}
                        do
                            echo "$(basename $f | sed 's/_R1.*\.fastq.gz//')\t$f" >> {output}
                        done
                    """

        rule import_raw_fastq_demux:
            input:
                output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_raw.manifest",
            output:
                output_dir+config["run_name"]+"-qaqc/raw_fastq.qza",
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                """
                qiime tools import \
                --type '{IMPORT_TYPE}' \
                --input-path {input[0]} \
                --output-path {output[0]} \
                --input-format {MANIFEST_FORMAT}      
                """
        ruleorder: cutadapt > import_fastq_demux
    elif config["trimmed_fastq_path"] != None:
        print(f"no manifest, not trimming {SEQUENCE_TYPE.lower()} reads\n")
        SAMPLES, =glob_wildcards(config["trimmed_fastq_path"]+"/{sample}.1.fastq.gz")
        ruleorder: import_fastq_demux > cutadapt
        ruleorder: make_manifest > make_manifest_file

#if IS_PAIRED:
#    primerF=config["fwd_primer"]
#    primerR=config["rev_primer"]
#    revcomp_primerF=Seq(config["fwd_primer"]).reverse_complement()
#    revcomp_primerR=Seq(config["rev_primer"]).reverse_complement()

rule raw_fastq_summary:
    input:
        output_dir+config["run_name"]+"-qaqc/raw_fastq.qza"
    output:
        output_dir+config["run_name"]+"-qaqc/stats/raw_fastq_summary.qzv"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"

rule cutadapt:
    input:
        output_dir+config["run_name"]+"-qaqc/raw_fastq.qza"
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq.qza",
        output_dir+config["run_name"]+"-qaqc/stats/cutadapt_summary.txt",
    params:
        discard=config['discard_untrimmed'],
        min_len=config['minimum_length'],
        tempTrim=output_dir+config["run_name"]+"-qaqc/trimmed_1.qza",
        primerF=config["fwd_primer"],
        primerR=config["rev_primer"],
    conda:
        "qiime2-amplicon-2024.10"
    threads: workflow.cores
    shell:
        """
        if [ "{IS_PAIRED}" = "True" ]; then
            revcomp_primerF=`echo {params.primerF} | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
            revcomp_primerR=`echo {params.primerR} | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
            qiime cutadapt trim-paired \
            --p-cores {threads} \
            --i-demultiplexed-sequences {input} \
            --p-adapter-f $revcomp_primerR \
            --p-adapter-r $revcomp_primerF \
            --p-match-read-wildcards \
            --p-match-adapter-wildcards \
            --p-minimum-length {params.min_len} \
            --verbose \
            --o-trimmed-sequences {params.tempTrim} 1> {output[1]}

            if [ {params.discard} = True ]; then
                qiime cutadapt trim-paired \
                --p-cores {threads} \
                --i-demultiplexed-sequences {params.tempTrim} \
                --p-front-f {params.primerF} \
                --p-front-r {params.primerR} \
                --p-match-read-wildcards \
                --p-match-adapter-wildcards \
                --p-discard-untrimmed \
                --p-minimum-length {params.min_len} \
                --verbose \
                --o-trimmed-sequences {output[0]} 1>> {output[1]}
            else
                qiime cutadapt trim-paired \
                --p-cores {threads} \
                --i-demultiplexed-sequences {params.tempTrim} \
                --p-front-f {params.primerF} \
                --p-front-r {params.primerR} \
                --p-match-read-wildcards \
                --p-match-adapter-wildcards \
                --p-minimum-length {params.min_len} \
                --verbose \
                --o-trimmed-sequences {output[0]} 1>> {output[1]}
            fi
        else
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
            fi
        fi
        rm {params.tempTrim}
        """

rule make_manifest:
    input:
        fwdreads=expand(output_dir+config["run_name"]+"-qaqc/trimmed/{sample}.1.fastq",sample=SAMPLES),
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_manifest"
    shell:
        """
        if [ "{IS_PAIRED}" = "True" ]; then
            echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > {output}
            for f in {input.fwdreads}
            do
                echo -e "$(basename $f | sed 's/.1.fastq//')\t$(pwd)/$f\t$(pwd)/$(echo $f | sed 's/.1.fastq/.2.fastq/g')" >> {output}
            done
        else
            echo -e "sample-id\tabsolute-filepath" > {output}
            for f in {input.fwdreads}
            do
                echo "$(basename $f | sed 's/.1.fastq//')\t$f" >> {output}
            done
        fi
        """

rule make_manifest_file:
    input:
        output_dir+"create_manifest.txt"
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_manifest"
    params:
        file=config["sample_manifest_file"]
    shell:
        """
        if [ ! -z "{params.file}" ]; then
            mkdir -p $(dirname {output})
            cp {params.file} {output}
        fi
        rm {input}
        """

rule import_fastq_demux:
    input:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_manifest"
    output:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq.qza"
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime tools import "
        f"--type '{IMPORT_TYPE}' "
        "--input-path {input} "
        "--output-path {output} "
        f"--input-format {MANIFEST_FORMAT}"

rule summarize_fastq_demux:
    input:
        output_dir+config["run_name"]+"-qaqc/"+config["run_name"]+"_fastq.qza"
    output:
        output_dir+config["run_name"]+"-qaqc/stats/fastq_summary.qzv",
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"
        #"cp

rule check_seq_qual_dropoff:
    input:
        R1=output_dir+config["run_name"]+"-qaqc/trimmed/fastqc_R1/multiqc_report_R1.html",
        R2=output_dir+config["run_name"]+"-qaqc/trimmed/fastqc_R2/multiqc_report_R2.html" if IS_PAIRED else [],
    output:
        output_dir+config["run_name"]+"-qaqc/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",
    params:
        seq_cutoff = config["seq_quality_cutoff"],
        seq_qual_R1 = output_dir+config["run_name"]+"-qaqc/trimmed/fastqc_R1/multiqc_report_R1_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt",
        seq_qual_R2 = output_dir+config["run_name"]+"-qaqc/trimmed/fastqc_R2/multiqc_report_R2_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt" if IS_PAIRED else None,
        stats_dir = output_dir+config["run_name"]+"-qaqc/stats/",
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        """
        mkdir -p {params.stats_dir}
        echo "R1 seq quality" > {output}
        python scripts/fastqc_per_base_sequence_quality_dropoff.py --input {params.seq_qual_R1} --cutoff {params.seq_cutoff} >> {output}
        if [ "{IS_PAIRED}" = "True" ]; then
            echo "R2 seq quality" >> {output}
            python scripts/fastqc_per_base_sequence_quality_dropoff.py --input {params.seq_qual_R2} --cutoff {params.seq_cutoff} >> {output}
        fi
        """
