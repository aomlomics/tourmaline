rule dada2_pe_denoise:
    """Run paired end dada2"""
    input:
        config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.fasta",
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest",
        config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt",

rule dada2_se_denoise:
    """Trim all reads with all supplied trimming parameters"""
    input:
        expand(config["run_name"]+"-samples/trimmed/fastqc_R1/{sample}.1_fastqc.html",sample=SAMPLES),
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest",
        config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt"

rule deblur_se_denoise:
    """Trim all reads with all supplied trimming parameters"""
    input:
        expand(config["run_name"]+"-samples/trimmed/fastqc_R1/{sample}.1_fastqc.html",sample=SAMPLES),
        config["run_name"]+"-samples/"+config["run_name"]+"_se.manifest",
        config["run_name"]+"-samples/trimmed/fastqc_R1/multiqc_report_R1.html",
        config["run_name"]+"-samples/stats/"+config["run_name"]+"-seq_qual_dropoff.txt"