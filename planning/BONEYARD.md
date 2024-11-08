# BONEYARD

## Trimming  
```
        #paste <(ls {input.fwdreads}) <(ls {input.revreads}) | while IFS=$'\t' read -r fwd rev
        #do
        #    echo "$(basename $fwd | sed 's/.1.fastq//'),$fwd,$(echo $rev | sed 's/.1.fastq/.2.fastq/g')" >> {output}
        #done
```

## Trying modules  

```
rule fastqc:
    input:
        config["run_name"]+"-trimmed/{sample}.{read}.fastq"
    output:
        html=config["run_name"]+"-trimmed/fastqc_R{read}/{sample}.{read}.html",
        #zip=config["run_name"]+"-trimmed/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet",
        outdir = config["run_name"]+"-trimmed/fastqc_R{read}"
    conda:
        "multiqc"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir};
        fastqc --threads {threads} {input} --outdir {params.outdir}
        """
```

```
rule fastqc:
    input:
        config["run_name"]+"-trimmed/{sample}.fastq"
    output:
        html=config["run_name"]+"-trimmed/qc/fastqc/{sample}.html",
        zip=config["run_name"]+"-trimmed/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v3.12.1/bio/fastqc"
```
Error:  snakemake --use-conda -s sample_test_paired.Snakefile --configfile $CONFIG --cores $CORES --conda-frontend conda

```
(snakemake) lutjanus:tourmaline katherine.silliman$ ./tourmaline.sh -s samples -c samp_config_test.yaml -n 5
Building DAG of jobs...
Traceback (most recent call last):
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/__init__.py", line 771, in snakemake
    success = workflow.execute(
              ^^^^^^^^^^^^^^^^^
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/workflow.py", line 924, in execute
    dag.create_conda_envs(
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/dag.py", line 349, in create_conda_envs
    env.create(dryrun)
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/deployment/conda.py", line 385, in create
    if self.pin_file:
       ^^^^^^^^^^^^^
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/common/__init__.py", line 217, in __get__
    value = self.method(instance)
            ^^^^^^^^^^^^^^^^^^^^^
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/deployment/conda.py", line 103, in pin_file
    f".{self.conda.platform}.pin.txt"
        ^^^^^^^^^^
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/common/__init__.py", line 217, in __get__
    value = self.method(instance)
            ^^^^^^^^^^^^^^^^^^^^^
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/deployment/conda.py", line 96, in conda
    return Conda(
           ^^^^^^
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/deployment/conda.py", line 667, in __init__
    self._check()
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/deployment/conda.py", line 736, in _check
    self._check_version()
  File "/Users/katherine.silliman/miniconda3/envs/snakemake/lib/python3.11/site-packages/snakemake/deployment/conda.py", line 745, in _check_version
    from packaging.version import Version
ModuleNotFoundError: No module named 'packaging'
```

–p-front-f is the forward primer (5’-3’), –p-front-r the reverse. I’ve also included –p-adapter-f, which is the reverse complement of the reverse primer, and –p-adapter-r, which is the reverse complement of the forward primer.

Trimming

fwd_primer, p-front-f:  GCCGGTAAAACTCGTGCCAGC;required
rev_primer, p-front-r:  CATAGTGGGGTATCTAATCCCAGTTTG;required
revcomp_rev, p-adapter-f: CAAACTGGGATTAGATACCCCACTATG;optional
revcomp_fwd,-p-adapter-r: GCTGGCACGAGTTTTACCGGC;optional

a_primer = f"{primerF};required...{revcomp_primerR};optional"
A_primer = f"{primerR};required...{revcomp_primerF};optional"

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../test-fastqs --input-format SingleLanePerSamplePairedEndFastqDirFmt --output-path demux-paired-end.qza

CasavaOneEightSingleLanePerSampleDirFmt
CasavaOneEightLanelessPerSampleDirFmt

qiime cutadapt trim-paired \
--p-cores 5 \
--i-demultiplexed-sequences demux-paired-end.qza \
--p-adapter-f CAAACTGGGATTAGATACCCCACTATG \
--p-adapter-r GCTGGCACGAGTTTTACCGGC \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 50 \
--verbose \
--o-trimmed-sequences trimmed_remove_primers_test_1.qza 1> ca_summary1

qiime cutadapt trim-paired \
--p-cores 5 \
--i-demultiplexed-sequences trimmed_remove_primers_test_1.qza \
--p-front-f GCCGGTAAAACTCGTGCCAGC \
--p-front-r CATAGTGGGGTATCTAATCCCAGTTTG \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-discard-untrimmed \
--p-minimum-length 50 \
--verbose \
--o-trimmed-sequences trimmed_remove_primers_test_2.qza 1> ca_summary2


    shell:
        "qiime metadata tabulate "
        "--m-input-file {input.repseqs} "
        "--m-input-file {input.taxonomy} "
        "--o-visualization merged-data.qzv; "
        "qiime tools export "
        "--input-path merged-data.qzv "
        "--output-path {output}; "
        "mv {output}/metadata.tsv temp; "
        "rm -r {output}; "
        "sed -e '2d' temp | sed '1 s|id\\t|featureid\\t|' | sed '1 s|Taxon|taxonomy|' | sed '1 s|Sequence|sequence|' > {output}; "
        "/bin/rm -r temp merged-data.qzv"
