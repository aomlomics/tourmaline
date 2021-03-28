import pandas as pd
import numpy as np
from qiime2 import Artifact
import matplotlib.pyplot as plt
import seaborn as sns
from tabulate import tabulate

# GLOBALS ----------------------------------------------------------------------

configfile: "config.yaml"

# RULEORDER DIRECTIVES ---------------------------------------------------------

ruleorder: summarize_fastq_demux_pe > summarize_fastq_demux_se
ruleorder: feature_classifier > import_taxonomy_to_qza
ruleorder: filter_taxonomy > unzip_taxonomy_to_tsv

# PSEUDO-RULES: DADA2 PAIRED-END -----------------------------------------------

rule dada2_pe_denoise:
    input:
        "01-imported/fastq_summary.qzv",
        "01-imported/fastq_counts_describe.md",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/table.qzv",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs_lengths_describe.md"

rule dada2_pe_taxonomy_unfiltered:
    input:
        "02-output-dada2-pe-unfiltered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-pe-unfiltered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-pe-unfiltered/01-taxonomy/taxa_barplot.qzv"

rule dada2_pe_diversity_unfiltered:
    input:
        "02-output-dada2-pe-unfiltered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-pe-unfiltered/02-alignment-tree/aligned_repseqs_gaps_describe.md",
        "02-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_properties.pdf",
        "02-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        "02-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_unassigned.tsv",
        "02-output-dada2-pe-unfiltered/03-alpha-diversity/alpha_rarefaction.qzv",
        "02-output-dada2-pe-unfiltered/03-alpha-diversity/faith_pd_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/03-alpha-diversity/observed_features_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/03-alpha-diversity/shannon_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/03-alpha-diversity/evenness_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/bray_curtis_emperor.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/bray_curtis_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/jaccard_emperor.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/jaccard_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/unweighted_unifrac_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/weighted_unifrac_emperor.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/weighted_unifrac_group_significance.qzv"

rule dada2_pe_report_unfiltered:
    input:
        "03-reports/report_dada2-pe_unfiltered.html"

rule dada2_pe_taxonomy_filtered:
    input:
        "02-output-dada2-pe-filtered/00-table-repseqs/table.qzv",
        "02-output-dada2-pe-filtered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-pe-filtered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-pe-filtered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-pe-filtered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-pe-filtered/00-table-repseqs/repseqs_lengths_describe.md",
        "02-output-dada2-pe-filtered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-pe-filtered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-pe-filtered/01-taxonomy/taxa_barplot.qzv"

rule dada2_pe_diversity_filtered:
    input:
        "02-output-dada2-pe-filtered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-pe-filtered/02-alignment-tree/aligned_repseqs_gaps_describe.md",
        "02-output-dada2-pe-filtered/02-alignment-tree/repseqs_properties.pdf",
        "02-output-dada2-pe-filtered/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        "02-output-dada2-pe-filtered/02-alignment-tree/repseqs_to_filter_unassigned.tsv",
        "02-output-dada2-pe-filtered/03-alpha-diversity/alpha_rarefaction.qzv",
        "02-output-dada2-pe-filtered/03-alpha-diversity/faith_pd_group_significance.qzv",
        "02-output-dada2-pe-filtered/03-alpha-diversity/observed_features_group_significance.qzv",
        "02-output-dada2-pe-filtered/03-alpha-diversity/shannon_group_significance.qzv",
        "02-output-dada2-pe-filtered/03-alpha-diversity/evenness_group_significance.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/bray_curtis_emperor.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/bray_curtis_group_significance.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/jaccard_emperor.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/jaccard_group_significance.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/unweighted_unifrac_group_significance.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/weighted_unifrac_emperor.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/weighted_unifrac_group_significance.qzv"

rule dada2_pe_report_filtered:
    input:
        "03-reports/report_dada2-pe_filtered.html"

# PSEUDO-RULES: DADA2 SINGLE-END -----------------------------------------------

rule dada2_se_denoise:
    input:
        "01-imported/fastq_summary.qzv",
        "01-imported/fastq_counts_describe.md",
        "02-output-dada2-se-unfiltered/00-table-repseqs/table.qzv",
        "02-output-dada2-se-unfiltered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-se-unfiltered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-se-unfiltered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-se-unfiltered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-se-unfiltered/00-table-repseqs/repseqs_lengths_describe.md"

rule dada2_se_taxonomy_unfiltered:
    input:
        "02-output-dada2-se-unfiltered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-se-unfiltered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-se-unfiltered/01-taxonomy/taxa_barplot.qzv"

rule dada2_se_diversity_unfiltered:
    input:
        "02-output-dada2-se-unfiltered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-se-unfiltered/02-alignment-tree/aligned_repseqs_gaps_describe.md",
        "02-output-dada2-se-unfiltered/02-alignment-tree/repseqs_properties.pdf",
        "02-output-dada2-se-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        "02-output-dada2-se-unfiltered/02-alignment-tree/repseqs_to_filter_unassigned.tsv",
        "02-output-dada2-se-unfiltered/03-alpha-diversity/alpha_rarefaction.qzv",
        "02-output-dada2-se-unfiltered/03-alpha-diversity/faith_pd_group_significance.qzv",
        "02-output-dada2-se-unfiltered/03-alpha-diversity/observed_features_group_significance.qzv",
        "02-output-dada2-se-unfiltered/03-alpha-diversity/shannon_group_significance.qzv",
        "02-output-dada2-se-unfiltered/03-alpha-diversity/evenness_group_significance.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/bray_curtis_emperor.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/bray_curtis_group_significance.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/jaccard_emperor.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/jaccard_group_significance.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/unweighted_unifrac_group_significance.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/weighted_unifrac_emperor.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/weighted_unifrac_group_significance.qzv"

rule dada2_se_report_unfiltered:
    input:
        "03-reports/report_dada2-se_unfiltered.html"

rule dada2_se_taxonomy_filtered:
    input:
        "02-output-dada2-se-filtered/00-table-repseqs/table.qzv",
        "02-output-dada2-se-filtered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-se-filtered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-se-filtered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-se-filtered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-se-filtered/00-table-repseqs/repseqs_lengths_describe.md",
        "02-output-dada2-se-filtered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-se-filtered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-se-filtered/01-taxonomy/taxa_barplot.qzv"

rule dada2_se_diversity_filtered:
    input:
        "02-output-dada2-se-filtered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-se-filtered/02-alignment-tree/aligned_repseqs_gaps_describe.md",
        "02-output-dada2-se-filtered/02-alignment-tree/repseqs_properties.pdf",
        "02-output-dada2-se-filtered/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        "02-output-dada2-se-filtered/02-alignment-tree/repseqs_to_filter_unassigned.tsv",
        "02-output-dada2-se-filtered/03-alpha-diversity/alpha_rarefaction.qzv",
        "02-output-dada2-se-filtered/03-alpha-diversity/faith_pd_group_significance.qzv",
        "02-output-dada2-se-filtered/03-alpha-diversity/observed_features_group_significance.qzv",
        "02-output-dada2-se-filtered/03-alpha-diversity/shannon_group_significance.qzv",
        "02-output-dada2-se-filtered/03-alpha-diversity/evenness_group_significance.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/bray_curtis_emperor.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/bray_curtis_group_significance.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/jaccard_emperor.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/jaccard_group_significance.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/unweighted_unifrac_group_significance.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/weighted_unifrac_emperor.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/weighted_unifrac_group_significance.qzv"

rule dada2_se_report_filtered:
    input:
        "03-reports/report_dada2-se_filtered.html"

# PSEUDO-RULES: DEBLUR SINGLE-END ----------------------------------------------

rule deblur_se_denoise:
    input:
        "01-imported/fastq_summary.qzv",
        "01-imported/fastq_counts_describe.md",
        "02-output-deblur-se-unfiltered/00-table-repseqs/table.qzv",
        "02-output-deblur-se-unfiltered/00-table-repseqs/table_summary_samples.txt",
        "02-output-deblur-se-unfiltered/00-table-repseqs/table_summary_features.txt",
        "02-output-deblur-se-unfiltered/00-table-repseqs/repseqs.qzv",
        "02-output-deblur-se-unfiltered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-deblur-se-unfiltered/00-table-repseqs/repseqs_lengths_describe.md"

rule deblur_se_taxonomy_unfiltered:
    input:
        "02-output-deblur-se-unfiltered/01-taxonomy/taxonomy.tsv",
        "02-output-deblur-se-unfiltered/01-taxonomy/taxonomy.qzv",
        "02-output-deblur-se-unfiltered/01-taxonomy/taxa_barplot.qzv"

rule deblur_se_diversity_unfiltered:
    input:
        "02-output-deblur-se-unfiltered/02-alignment-tree/rooted_tree.qzv",
        "02-output-deblur-se-unfiltered/02-alignment-tree/aligned_repseqs_gaps_describe.md",
        "02-output-deblur-se-unfiltered/02-alignment-tree/repseqs_properties.pdf",
        "02-output-deblur-se-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        "02-output-deblur-se-unfiltered/02-alignment-tree/repseqs_to_filter_unassigned.tsv",
        "02-output-deblur-se-unfiltered/03-alpha-diversity/alpha_rarefaction.qzv",
        "02-output-deblur-se-unfiltered/03-alpha-diversity/faith_pd_group_significance.qzv",
        "02-output-deblur-se-unfiltered/03-alpha-diversity/observed_features_group_significance.qzv",
        "02-output-deblur-se-unfiltered/03-alpha-diversity/shannon_group_significance.qzv",
        "02-output-deblur-se-unfiltered/03-alpha-diversity/evenness_group_significance.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/bray_curtis_emperor.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/bray_curtis_group_significance.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/jaccard_emperor.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/jaccard_group_significance.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/unweighted_unifrac_group_significance.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/weighted_unifrac_emperor.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/weighted_unifrac_group_significance.qzv"

rule deblur_se_report_unfiltered:
    input:
        "03-reports/report_deblur-se_unfiltered.html"

rule deblur_se_taxonomy_filtered:
    input:
        "02-output-deblur-se-filtered/00-table-repseqs/table.qzv",
        "02-output-deblur-se-filtered/00-table-repseqs/table_summary_samples.txt",
        "02-output-deblur-se-filtered/00-table-repseqs/table_summary_features.txt",
        "02-output-deblur-se-filtered/00-table-repseqs/repseqs.qzv",
        "02-output-deblur-se-filtered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-deblur-se-filtered/00-table-repseqs/repseqs_lengths_describe.md",
        "02-output-deblur-se-filtered/01-taxonomy/taxonomy.tsv",
        "02-output-deblur-se-filtered/01-taxonomy/taxonomy.qzv",
        "02-output-deblur-se-filtered/01-taxonomy/taxa_barplot.qzv"

rule deblur_se_diversity_filtered:
    input:
        "02-output-deblur-se-filtered/02-alignment-tree/rooted_tree.qzv",
        "02-output-deblur-se-filtered/02-alignment-tree/aligned_repseqs_gaps_describe.md",
        "02-output-deblur-se-filtered/02-alignment-tree/repseqs_properties.pdf",
        "02-output-deblur-se-filtered/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        "02-output-deblur-se-filtered/02-alignment-tree/repseqs_to_filter_unassigned.tsv",
        "02-output-deblur-se-filtered/03-alpha-diversity/alpha_rarefaction.qzv",
        "02-output-deblur-se-filtered/03-alpha-diversity/faith_pd_group_significance.qzv",
        "02-output-deblur-se-filtered/03-alpha-diversity/observed_features_group_significance.qzv",
        "02-output-deblur-se-filtered/03-alpha-diversity/shannon_group_significance.qzv",
        "02-output-deblur-se-filtered/03-alpha-diversity/evenness_group_significance.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/bray_curtis_emperor.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/bray_curtis_group_significance.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/jaccard_emperor.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/jaccard_group_significance.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/unweighted_unifrac_group_significance.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/weighted_unifrac_emperor.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/weighted_unifrac_group_significance.qzv"

rule deblur_se_report_filtered:
    input:
        "03-reports/report_deblur-se_filtered.html"

# RULES: IMPORT ----------------------------------------------------------------

rule import_ref_seqs:
    input:
        config["refseqs_fna"]
    output:
        config["refseqs_qza"]
    shell:
        "qiime tools import "
        "--type 'FeatureData[Sequence]' "
        "--input-path {input} "
        "--output-path {output}"

rule import_ref_tax:
    input:
        config["reftax_tsv"]
    output:
        config["reftax_qza"]
    shell:
        "qiime tools import "
        "--type 'FeatureData[Taxonomy]' "
        "--input-format HeaderlessTSVTaxonomyFormat "
        "--input-path {input} "
        "--output-path {output}"

rule import_fastq_demux_pe:
    input:
        config["manifest_pe"]
    output:
        config["fastq_pe_qza"]
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input} "
        "--output-path {output} "
        "--input-format PairedEndFastqManifestPhred33"

rule import_fastq_demux_se:
    input:
        config["manifest_se"]
    output:
        config["fastq_se_qza"]
    shell:
        "qiime tools import "
        "--type 'SampleData[SequencesWithQuality]' "
        "--input-path {input} "
        "--output-path {output} "
        "--input-format SingleEndFastqManifestPhred33"

# RULES: SUMMARIZE FASTQ SEQUENCES ---------------------------------------------

rule summarize_fastq_demux_pe:
    input:
        config["fastq_pe_qza"]
    output:
        "01-imported/fastq_summary.qzv"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"

rule summarize_fastq_demux_se:
    input:
        config["fastq_se_qza"]
    output:
        "01-imported/fastq_summary.qzv"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"

rule unzip_fastq_summary:
    input:
        "01-imported/fastq_summary.qzv"
    output:
        "01-imported/fastq_counts.tsv"
    shell:
        "unzip -qq -o {input} -d temp0; "
        "mv temp0/*/data/per-sample-fastq-counts.tsv {output}; "
        "/bin/rm -r temp0"

rule describe_fastq_counts:
    input:
        "01-imported/fastq_counts.tsv"
    output:
        "01-imported/fastq_counts_describe.md"
    run:
        s = pd.read_csv(input[0], sep='\t', index_col=0)
        t = s.describe()
        outstr = tabulate(pd.DataFrame(t.iloc[1:,0]), tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0,0].astype(int), 'Fastq sequences per sample'])
        with open(output[0], 'w') as target:
            target.write(outstr)
            target.write('\n')

# RULES: DENOISE ---------------------------------------------------------------

rule denoise_dada2_pe:
    input:
        config["fastq_pe_qza"]
    params:
        trunclenf=config["dada2pe_trunc_len_f"],
        trunclenr=config["dada2pe_trunc_len_r"],
        trimleftf=config["dada2pe_trim_left_f"],
        trimleftr=config["dada2pe_trim_left_r"],
        maxeef=config["dada2pe_max_ee_f"],
        maxeer=config["dada2pe_max_ee_r"],
        truncq=config["dada2pe_trunc_q"],
        poolingmethod=config["dada2pe_pooling_method"],        
        chimeramethod=config["dada2pe_chimera_method"],
        minfoldparentoverabundance=config["dada2pe_min_fold_parent_over_abundance"],
        nreadslearn=config["dada2pe_n_reads_learn"],
        hashedfeatureids=config["dada2pe_hashed_feature_ids"]
    output:
        table="02-output-dada2-pe-unfiltered/00-table-repseqs/table.qza",
        repseqs="02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs.qza",
        stats="02-output-dada2-pe-unfiltered/00-table-repseqs/dada2_stats.qza"
    threads: config["dada2pe_threads"]
    shell:
        "qiime dada2 denoise-paired "
        "--i-demultiplexed-seqs {input} "
        "--p-trunc-len-f {params.trunclenf} "
        "--p-trunc-len-r {params.trunclenr} "
        "--p-trim-left-f {params.trimleftf} "
        "--p-trim-left-r {params.trimleftr} "
        "--p-max-ee-f {params.maxeef} "
        "--p-max-ee-r {params.maxeer} "
        "--p-trunc-q {params.truncq} "
        "--p-pooling-method {params.poolingmethod} "
        "--p-chimera-method {params.chimeramethod} "
        "--p-min-fold-parent-over-abundance {params.minfoldparentoverabundance} "
        "--p-n-reads-learn {params.nreadslearn} "
        "--p-n-threads {threads} "
        "{params.hashedfeatureids} "
        "--o-table {output.table} "
        "--o-representative-sequences {output.repseqs} "
        "--o-denoising-stats {output.stats} "
        "--verbose"

rule denoise_dada2_se:
    input:
        config["fastq_se_qza"]
    params:
        trunclen=config["dada2se_trunc_len"],
        trimleft=config["dada2se_trim_left"],
        maxee=config["dada2se_max_ee"],
        truncq=config["dada2se_trunc_q"],
        poolingmethod=config["dada2se_pooling_method"],        
        chimeramethod=config["dada2se_chimera_method"],
        minfoldparentoverabundance=config["dada2se_min_fold_parent_over_abundance"],
        nreadslearn=config["dada2se_n_reads_learn"],
        hashedfeatureids=config["dada2se_hashed_feature_ids"]
    output:
        table="02-output-dada2-se-unfiltered/00-table-repseqs/table.qza",
        repseqs="02-output-dada2-se-unfiltered/00-table-repseqs/repseqs.qza",
        stats="02-output-dada2-se-unfiltered/00-table-repseqs/dada2_stats.qza"
    threads: config["dada2se_threads"]
    shell:
        "qiime dada2 denoise-single "
        "--i-demultiplexed-seqs {input} "
        "--p-trunc-len {params.trunclen} "
        "--p-trim-left {params.trimleft} "
        "--p-max-ee {params.maxee} "
        "--p-trunc-q {params.truncq} "
        "--p-pooling-method {params.poolingmethod} "
        "--p-chimera-method {params.chimeramethod} "
        "--p-min-fold-parent-over-abundance {params.minfoldparentoverabundance} "
        "--p-n-reads-learn {params.nreadslearn} "
        "--p-n-threads {threads} "
        "{params.hashedfeatureids} "
        "--o-table {output.table} "
        "--o-representative-sequences {output.repseqs} "
        "--o-denoising-stats {output.stats} "
        "--verbose"

rule denoise_deblur_se:
    input:
        seqs=config["fastq_se_qza"],
        refseqs=config["refseqs_qza"]
    params:
        trimlength=config["deblur_trim_length"],
        samplestats=config["deblur_sample_stats"],
        meanerror=config["deblur_mean_error"],
        indelprob=config["deblur_indel_prob"],
        indelmax=config["deblur_indel_max"],
        minreads=config["deblur_min_reads"],
        minsize=config["deblur_min_size"],
        hashedfeatureids=config["deblur_hashed_feature_ids"]
    output:
        table="02-output-deblur-se-unfiltered/00-table-repseqs/table.qza",
        repseqs="02-output-deblur-se-unfiltered/00-table-repseqs/repseqs.qza",
        stats="02-output-deblur-se-unfiltered/00-table-repseqs/deblur_stats.qza"
    threads: config["deblur_threads"]
    shell:
        "qiime deblur denoise-other "
        "--i-demultiplexed-seqs {input.seqs} "
        "--i-reference-seqs {input.refseqs} "
        "--p-trim-length {params.trimlength} "
        "{params.samplestats} "
        "--p-mean-error {params.meanerror} "
        "--p-indel-prob {params.indelprob} "
        "--p-indel-max {params.indelmax} "
        "--p-min-reads {params.minreads} "
        "--p-min-size {params.minsize} "
        "--p-jobs-to-start {threads} "
        "{params.hashedfeatureids} "
        "--o-table {output.table} "
        "--o-representative-sequences {output.repseqs} "
        "--o-stats {output.stats} "
        "--verbose"

# RULES: SUMMARIZE FEATURE TABLE -----------------------------------------------

rule summarize_feature_table:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        metadata=config["metadata"]
    output:
        "02-output-{method}-{filter}/00-table-repseqs/table.qzv"
    shell:
        "qiime feature-table summarize "
        "--i-table {input.table} "
        "--m-sample-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule unzip_table_to_biom:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/table.qza"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/table.biom"
    shell:
        "unzip -qq -o {input} -d temp1; "
        "mv temp1/*/data/feature-table.biom {output}; "
        "/bin/rm -r temp1"

rule summarize_biom_samples:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/table.biom"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/table_summary_samples.txt"
    shell:
        "biom summarize-table "
        "--input-fp {input} "
        "--output-fp {output}; "
        "cat {output} | sed 's/observation/feature/g' | sed 's/.000$//' > temp; "
        "mv temp {output}"

rule summarize_biom_features:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/table.biom"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/table_summary_features.txt"
    shell:
        "biom summarize-table "
        "--observations "
        "--input-fp {input} "
        "--output-fp {output}; "
        "cat {output} | sed 's/observation/feature/g' | sed 's/Counts\/sample/Counts\/feature/g' | sed 's/.000$//' > temp; "
        "mv temp {output}"

rule visualize_repseqs:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.qza"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.qzv"
    shell:
        "qiime feature-table tabulate-seqs "
        "--i-data {input} "
        "--o-visualization {output}"

rule unzip_repseqs_to_fasta:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.qza"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
    shell:
        "unzip -qq -o {input} -d temp2; "
        "mv temp2/*/data/dna-sequences.fasta {output}; "
        "/bin/rm -r temp2"

rule repseqs_detect_amplicon_locus:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_amplicon_type.txt"
    shell:
        "python scripts/detect_amplicon_locus.py -i {input} > {output}"

rule repseqs_lengths:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths.txt"
    shell:
        "perl scripts/fastaLengths.pl {input} > {output}"

rule repseqs_lengths_describe:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths.txt"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths_describe.md"
    run:
        s = pd.read_csv(input[0], header=None, index_col=0)
        t = s.describe()
        outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0].values[0].astype(int), 'Sequence length'])
        with open(output[0], 'w') as target:
            target.write(outstr)
            target.write('\n')

# RULES: TAXONOMY --------------------------------------------------------------

rule feature_classifier:
    input:
        refseqs=config["refseqs_qza"],
        reftax=config["reftax_qza"],
        repseqs="02-output-{method}-unfiltered/00-table-repseqs/repseqs.qza"
    output:
        "02-output-{method}-unfiltered/01-taxonomy/taxonomy.qza"
    params:
        classifymethod=config["classify_method"]
    threads: config["feature_classifier_threads"]
    shell:
        "echo classify_method: {params.classifymethod}; "
        "if [ {params.classifymethod} = 'naive-bayes' ]; then "
        "if [ ! -f 01-imported/classifier.qza ]; then "
        "qiime feature-classifier fit-classifier-naive-bayes "
        "--i-reference-reads {input.refseqs} "
        "--i-reference-taxonomy {input.reftax} "
        "--o-classifier 01-imported/classifier.qza; "
        "fi; "
        "qiime feature-classifier classify-sklearn "
        "--i-classifier 01-imported/classifier.qza "
        "--i-reads {input.repseqs} "
        "--o-classification {output} "
        "--p-n-jobs {threads}; "
        "elif [ {params.classifymethod} = 'consensus-blast' ]; then "
        "qiime feature-classifier classify-consensus-blast "
        "--i-reference-reads {input.refseqs} "
        "--i-reference-taxonomy {input.reftax} "
        "--i-query {input.repseqs} "
        "--o-classification {output}; "
        "fi"

rule visualize_taxonomy:
    input:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qza"
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qzv"
    shell:
        "qiime metadata tabulate "
        "--m-input-file {input} "
        "--o-visualization {output}"

rule taxa_barplot:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
        metadata=config["metadata"]
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxa_barplot.qzv"
    shell:
        "qiime taxa barplot "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--m-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule unzip_taxonomy_to_tsv:
    input:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qza"
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.tsv"
    shell:
        "unzip -qq -o {input} -d temp3; "
        "mv temp3/*/data/taxonomy.tsv {output}; "
        "/bin/rm -r temp3"

rule import_taxonomy_to_qza:
    input:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.tsv"
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qza"
    shell:
        "qiime tools import "
        "--type 'FeatureData[Taxonomy]' "
        "--input-format TSVTaxonomyFormat "
        "--input-path {input} "
        "--output-path {output}"

# RULES: ALIGNMENT & TREE ------------------------------------------------------

# OPTION 1: MUSCLE - default (leave lines 653-666 uncommented and lines 670-683 and 687-717 commented)

rule alignment_muscle:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
    output:
        aln_fasta="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta",
        aln_qza="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza"
    shell:
        "muscle -maxiters 2 -diags "
        "-in {input} "
        "-out {output.aln_fasta}; "
        "qiime tools import "
        "--type 'FeatureData[AlignedSequence]' "
        "--input-path {output.aln_fasta} "
        "--output-path {output.aln_qza}"

# OPTION 2: Clustal Omega - comment lines 653-666 and uncomment lines 670-683 (leave lines 687-717 commented)

# rule alignment_clustalo:
#     input:
#         "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
#     output:
#         aln_fasta="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta",
#         aln_qza="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza"
#     shell:
#         "clustalo --verbose --force "
#         "--in {input} "
#         "--out {output.aln_fasta}; "
#         "qiime tools import "
#         "--type 'FeatureData[AlignedSequence]' "
#         "--input-path {output.aln_fasta} "
#         "--output-path {output.aln_qza}"

# OPTION 3: MAFFT with masking - comment lines 653-666 and uncomment lines 687-717 (leave lines 670-683 commented)

# rule alignment_mafft:
#     input:
#         "02-output-{method}-{filter}/00-table-repseqs/repseqs.qza"
#     output:
#         "02-output-{method}-{filter}/02-alignment-tree/unmasked_aligned_repseqs.qza"
#     threads: config["alignment_mafft_threads"]
#     shell:
#         "qiime alignment mafft "
#         "--i-sequences {input} "
#         "--o-alignment {output} "
#         "--p-n-threads {threads}"

# rule alignment_mask:
#     input:
#         "02-output-{method}-{filter}/02-alignment-tree/unmasked_aligned_repseqs.qza"
#     output:
#         "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza"
#     shell:
#         "qiime alimvgnment mask "
#         "--i-alignment {input} "
#         "--o-masked-alignment {output}"

# rule unzip_alignment_to_fasta:
#     input:
#         "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza"
#     output:
#         "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta"
#     shell:
#         "unzip -qq -o {input} -d temp4; "
#         "mv temp4/*/data/aligned-dna-sequences.fasta {output}; "
#         "/bin/rm -r temp4"

rule phylogeny_fasttree:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/unrooted_tree.qza"
    threads: config["phylogeny_fasttree_threads"]
    shell:
        "qiime phylogeny fasttree "
        "--i-alignment {input} "
        "--o-tree {output} "
        "--p-n-threads {threads}"

rule phylogeny_midpoint_root:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/unrooted_tree.qza"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qza"
    shell:
        "qiime phylogeny midpoint-root "
        "--i-tree {input} "
        "--o-rooted-tree {output}"

rule visualize_tree:
    input:
        tree="02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qza",
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        metadata=config["metadata"],
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
        outliers="02-output-{method}-{filter}/02-alignment-tree/outliers.qza"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qzv"
    shell:
        "qiime empress community-plot "
        "--i-tree {input.tree} "
        "--i-feature-table {input.table} "
        "--m-sample-metadata-file {input.metadata} "
        "--m-feature-metadata-file {input.taxonomy} "
        "--m-feature-metadata-file {input.outliers} "
        "--o-visualization {output}"

rule alignment_count_gaps:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps.txt"
    shell:
        "cat {input} | grep -v '>' | sed 's/[^-]//g' | awk '{{ print length }}' > {output}"

rule alignment_gaps_describe:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps.txt"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps_describe.md"
    run:
        gaps = pd.read_csv(input[0], header=None)
        t = gaps.describe()
        outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0].values[0].astype(int), 'Alignment gaps per sequence'])
        with open(output[0], 'w') as target:
            target.write(outstr)
            target.write('\n')

rule alignment_detect_outliers:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta"
    params:
        metric = config["odseq_distance_metric"],
        replicates = config["odseq_bootstrap_replicates"],
        threshold = config["odseq_threshold"]
    output:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_outliers.tsv"
    shell:
        "Rscript --vanilla scripts/run_odseq.R {input} {params.metric} {params.replicates} {params.threshold} {output}"

rule tabulate_plot_repseq_properties:
    input:
        lengths="02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths.txt",
        gaps="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps.txt",
        outliers="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_outliers.tsv",
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza"
    output:
        proptsv="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties.tsv",
        propdescribe="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties_describe.md",
        proppdf="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties.pdf",
        outliersforqza="02-output-{method}-{filter}/02-alignment-tree/outliers.tsv"
    run:
        lengths = pd.read_csv(input['lengths'], header=None)
        gaps = pd.read_csv(input['gaps'], header=None)
        outliers = pd.read_csv(input['outliers'], header=None, sep='\t')
        taxonomy = Artifact.load(input['taxonomy'])
        taxonomydf = taxonomy.view(view_type=pd.DataFrame)
        taxonomydf['level_1'] = [x.split(';')[0] for x in taxonomydf['Taxon']]
        table = Artifact.load(input['table'])
        tabledf = table.view(view_type=pd.DataFrame)
        merged = pd.concat([lengths, gaps, outliers[1], taxonomydf['Taxon'].reset_index(drop=True), taxonomydf['level_1'].reset_index(drop=True), tabledf.sum().reset_index(drop=True)],
                           axis=1, ignore_index=True)
        merged.columns = ['featureid', 'length', 'gaps', 'outlier', 'taxonomy', 'taxonomy_level_1', 'observations']
        merged['log10(observations)'] = [np.log10(x) for x in merged['observations']]
        merged.sort_values('log10(observations)', ascending=False, inplace=True)
        merged.to_csv(output['proptsv'], index=False, sep='\t')
        t = merged.describe()
        tcolumns = t.columns
        tcolumns = tcolumns.insert(0, 'Statistic (n=%s)' % t.iloc[0].values[0].astype(int))
        outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=tcolumns)
        with open(output['propdescribe'], 'w') as target:
            target.write(outstr)
            target.write('\n')
        g = sns.relplot(data=merged, x='length', y='gaps', col='outlier', hue='taxonomy_level_1', size='log10(observations)', sizes=(1,500), edgecolor = 'none', alpha=0.7)
        g.set_axis_labels('length (bp) not including gaps', 'gaps (bp) in multiple sequence alignment')
        plt.savefig(output['proppdf'], bbox_inches='tight')
        outliers.columns = ['Feature ID', 'Outlier']
        outliers = outliers*1
        outliers.to_csv(output['outliersforqza'], index=False, sep='\t')

rule import_outliers_to_qza:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/outliers.tsv"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/outliers.qza"
    shell:
        "qiime tools import "
        "--type 'FeatureData[Importance]' "
        "--input-path {input} "
        "--output-path {output}"

rule tabulate_repseqs_to_filter:
    input:
        proptsv="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties.tsv"
    output:
        outliers="02-output-{method}-{filter}/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        unassigned="02-output-{method}-{filter}/02-alignment-tree/repseqs_to_filter_unassigned.tsv"
    shell:
        "cat {input.proptsv} | grep -i 'outlier\|true' | cut -f1,4 > {output.outliers}; "
        "cat {input.proptsv} | grep -i 'taxonomy\|unassigned' | cut -f1,5 > {output.unassigned}"

# RULES: FILTER ----------------------------------------------------------------

rule filter_sequences_by_taxonomy:
    input:
        repseqs="02-output-{method}-unfiltered/00-table-repseqs/repseqs.qza",
        taxonomy="02-output-{method}-unfiltered/01-taxonomy/taxonomy.qza"
    params:
        excludeterms=config["exclude_terms"]
    output:
        "02-output-{method}-filtered/00-table-repseqs/repseqs_pre_id_filtering.qza"
    shell:
        "qiime taxa filter-seqs "
        "--i-sequences {input.repseqs} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude {params.excludeterms} "
        "--o-filtered-sequences {output}"

rule filter_sequences_by_id:
    input:
        repseqs="02-output-{method}-filtered/00-table-repseqs/repseqs_pre_id_filtering.qza",
        repseqstofilter="00-data/repseqs_to_filter_{method}.tsv"
    output:
        "02-output-{method}-filtered/00-table-repseqs/repseqs.qza"
    shell:
        "qiime feature-table filter-seqs "
        "--i-data {input.repseqs} "
        "--m-metadata-file {input.repseqstofilter} "
        "--p-exclude-ids "
        "--o-filtered-data {output}"

rule filter_table:
    input:
        table="02-output-{method}-unfiltered/00-table-repseqs/table.qza",
        filteredseqs="02-output-{method}-filtered/00-table-repseqs/repseqs.qza"
    output:
        "02-output-{method}-filtered/00-table-repseqs/table.qza",
    shell:
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--m-metadata-file {input.filteredseqs} "
        "--o-filtered-table {output}"

rule filter_taxonomy:
    input:
        taxonomy="02-output-{method}-unfiltered/01-taxonomy/taxonomy.tsv",
        repseqstofilter="00-data/repseqs_to_filter_{method}.tsv"
    params:
        excludeterms=config["exclude_terms"]
    output:
        taxonomy="02-output-{method}-filtered/01-taxonomy/taxonomy.tsv"
    run:
        df_taxonomy = pd.read_csv(input['taxonomy'], sep='\t', index_col=0)
        exclude_terms = params['excludeterms'].split(',')
        exclude_ids_taxa = df_taxonomy.index[[any(x.lower() in y.lower() for x in exclude_terms) for y in df_taxonomy['Taxon']]]
        df_repseqs_to_filter = pd.read_csv(input['repseqstofilter'], sep='\t')
        exclude_ids_direct = df_repseqs_to_filter['featureid'].values
        keep_ids = set(df_taxonomy.index) - set(exclude_ids_taxa) - set(exclude_ids_direct)
        df_taxonomy_filtered = df_taxonomy.loc[list(keep_ids)]
        df_taxonomy_filtered.to_csv(output['taxonomy'], sep='\t')

# RULES: DIVERSITY -------------------------------------------------------------

rule diversity_alpha_rarefaction:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        phylogeny="02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qza",
        metadata=config["metadata"]
    params:
        maxdepth=config["alpha_max_depth"]
    output:
        "02-output-{method}-{filter}/03-alpha-diversity/alpha_rarefaction.qzv"
    shell:
        "qiime diversity alpha-rarefaction "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-max-depth {params.maxdepth} "
        "--p-metrics faith_pd "
        "--p-metrics observed_features "
        "--p-metrics shannon "
        "--p-metrics pielou_e "
        "--m-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule diversity_core_metrics_phylogenetic:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        phylogeny="02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qza",
        metadata=config["metadata"]
    params:
        samplingdepth=config["core_sampling_depth"]
    output:
        rarefiedtable="02-output-{method}-{filter}/03-alpha-diversity/rarefied_table.qza",
        faithpdvector="02-output-{method}-{filter}/03-alpha-diversity/faith_pd_vector.qza",
        observedfeaturesvector="02-output-{method}-{filter}/03-alpha-diversity/observed_features_vector.qza",
        shannonvector="02-output-{method}-{filter}/03-alpha-diversity/shannon_vector.qza",
        evennessvector="02-output-{method}-{filter}/03-alpha-diversity/evenness_vector.qza",
        unweightedunifracdistancematrix="02-output-{method}-{filter}/04-beta-diversity/unweighted_unifrac_distance_matrix.qza",
        weightedunifracdistancematrix="02-output-{method}-{filter}/04-beta-diversity/weighted_unifrac_distance_matrix.qza",
        jaccarddistancematrix="02-output-{method}-{filter}/04-beta-diversity/jaccard_distance_matrix.qza",
        braycurtisdistancematrix="02-output-{method}-{filter}/04-beta-diversity/bray_curtis_distance_matrix.qza",
        unweightedunifracpcoaresults="02-output-{method}-{filter}/04-beta-diversity/unweighted_unifrac_pcoa_results.qza",
        weightedunifracpcoaresults="02-output-{method}-{filter}/04-beta-diversity/weighted_unifrac_pcoa_results.qza",
        jaccardpcoaresults="02-output-{method}-{filter}/04-beta-diversity/jaccard_pcoa_results.qza",
        braycurtispcoaresults="02-output-{method}-{filter}/04-beta-diversity/bray_curtis_pcoa_results.qza",
        unweightedunifracemperor="02-output-{method}-{filter}/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        weightedunifracemperor="02-output-{method}-{filter}/04-beta-diversity/weighted_unifrac_emperor.qzv",
        jaccardemperor="02-output-{method}-{filter}/04-beta-diversity/jaccard_emperor.qzv",
        braycurtisemperor="02-output-{method}-{filter}/04-beta-diversity/bray_curtis_emperor.qzv"
    threads: config["diversity_core_metrics_phylogenetic_threads"]
    shell:
        "qiime diversity core-metrics-phylogenetic "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-sampling-depth {params.samplingdepth} "
        "--m-metadata-file {input.metadata} "
        "--o-rarefied-table {output.rarefiedtable} "
        "--o-faith-pd-vector {output.faithpdvector} "
        "--o-observed-features-vector {output.observedfeaturesvector} "
        "--o-shannon-vector {output.shannonvector} "
        "--o-evenness-vector {output.evennessvector} "
        "--o-unweighted-unifrac-distance-matrix {output.unweightedunifracdistancematrix} "
        "--o-weighted-unifrac-distance-matrix {output.weightedunifracdistancematrix} "
        "--o-jaccard-distance-matrix {output.jaccarddistancematrix} "
        "--o-bray-curtis-distance-matrix {output.braycurtisdistancematrix} "
        "--o-unweighted-unifrac-pcoa-results {output.unweightedunifracpcoaresults} "
        "--o-weighted-unifrac-pcoa-results {output.weightedunifracpcoaresults} "
        "--o-jaccard-pcoa-results {output.jaccardpcoaresults} "
        "--o-bray-curtis-pcoa-results {output.braycurtispcoaresults} "
        "--o-unweighted-unifrac-emperor {output.unweightedunifracemperor} "
        "--o-weighted-unifrac-emperor {output.weightedunifracemperor} "
        "--o-jaccard-emperor {output.jaccardemperor} "
        "--o-bray-curtis-emperor {output.braycurtisemperor} "
        "--p-n-jobs-or-threads {threads}"

rule diversity_alpha_group_significance:
    input:
        alphadiversity="02-output-{method}-{filter}/03-alpha-diversity/{metric}_vector.qza",
        metadata=config["metadata"]
    output:
        "02-output-{method}-{filter}/03-alpha-diversity/{metric}_group_significance.qzv"
    shell:
        "qiime diversity alpha-group-significance "
        "--i-alpha-diversity {input.alphadiversity} "
        "--m-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule diversity_beta_group_significance:
    input:
        distancematrix="02-output-{method}-{filter}/04-beta-diversity/{metric}_distance_matrix.qza",
        metadata=config["metadata"]
    params:
        column=config["beta_group_column"],
        method=config["beta_group_method"],
        pairwise=config["beta_group_pairwise"]
    output:
        "02-output-{method}-{filter}/04-beta-diversity/{metric}_group_significance.qzv"
    shell:
        "qiime diversity beta-group-significance "
        "--i-distance-matrix {input.distancematrix} "
        "--m-metadata-file {input.metadata} "
        "--m-metadata-column {params.column} "
        "--p-method {params.method} "
        "{params.pairwise} "
        "--o-visualization {output}"

# RULES: REPORT ----------------------------------------------------------------

rule summarize_metadata:
    input:
        metadata=config["metadata"]
    output:
        "03-reports/metadata_summary.md"
    run:
        df = pd.read_csv(input.metadata, sep='\t')
        cols = df.columns
        df2 = pd.DataFrame(columns =[0,1], index=cols)
        for col in cols:
            if col in df.columns:
                vc = df[col].value_counts()
                if vc.index.shape == (0,):
                    df2.loc[col, 0] = '(no values in column)'
                    df2.loc[col, 1] = '--'
                else:
                    df2.loc[col, 0] = vc.index[0]
                    df2.loc[col, 1] = vc.values[0]
            else:
                df2.loc[col, 0] = '(column not provided)'
                df2.loc[col, 1] = '--'
        df2.columns = ['Most common value', 'Count']
        df2.index.name = 'Column name'
        outstr = tabulate(df2, tablefmt="pipe", headers="keys")
        with open(output[0], 'w') as target:
            target.write(outstr)
            target.write('\n')

rule generate_report_md:
    input:
        configfile="config.yaml",
        mdsummary="03-reports/metadata_summary.md",
        fastq="01-imported/fastq_counts_describe.md",
        amplicontype="02-output-{method}-{filter}/00-table-repseqs/repseqs_amplicon_type.txt",
        repseqstsv="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties.tsv",
        repseqsdescribe="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties_describe.md",
        repseqspdf="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties.pdf",
        repseqstree="02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qzv",
        repseqsoutliers="02-output-{method}-{filter}/02-alignment-tree/repseqs_to_filter_outliers.tsv",
        repseqsunassigned="02-output-{method}-{filter}/02-alignment-tree/repseqs_to_filter_unassigned.tsv",
        samples="02-output-{method}-{filter}/00-table-repseqs/table_summary_samples.txt",
        features="02-output-{method}-{filter}/00-table-repseqs/table_summary_features.txt",
        visfastq="01-imported/fastq_summary.qzv",
        visrepseqs="02-output-{method}-{filter}/00-table-repseqs/repseqs.qzv",
        tsvtaxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.tsv",
        vistaxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qzv",
        vistable="02-output-{method}-{filter}/00-table-repseqs/table.qzv",
        vistaxbar="02-output-{method}-{filter}/01-taxonomy/taxa_barplot.qzv",
        visalpharare="02-output-{method}-{filter}/03-alpha-diversity/alpha_rarefaction.qzv",
        visevengs="02-output-{method}-{filter}/03-alpha-diversity/evenness_group_significance.qzv",
        visfaithgs="02-output-{method}-{filter}/03-alpha-diversity/faith_pd_group_significance.qzv",
        visobsfeaturesgs="02-output-{method}-{filter}/03-alpha-diversity/observed_features_group_significance.qzv",
        visshannongs="02-output-{method}-{filter}/03-alpha-diversity/shannon_group_significance.qzv",
        visbcemp="02-output-{method}-{filter}/04-beta-diversity/bray_curtis_emperor.qzv",
        visbcgs="02-output-{method}-{filter}/04-beta-diversity/bray_curtis_group_significance.qzv",
        visjacemp="02-output-{method}-{filter}/04-beta-diversity/jaccard_emperor.qzv",
        visjacgs="02-output-{method}-{filter}/04-beta-diversity/jaccard_group_significance.qzv",
        viswuemp="02-output-{method}-{filter}/04-beta-diversity/weighted_unifrac_emperor.qzv",
        viswugs="02-output-{method}-{filter}/04-beta-diversity/weighted_unifrac_group_significance.qzv",
        visuwuemp="02-output-{method}-{filter}/04-beta-diversity/unweighted_unifrac_emperor.qzv",
        visuwugs="02-output-{method}-{filter}/04-beta-diversity/unweighted_unifrac_group_significance.qzv"
    output:
        "03-reports/report_{method}_{filter}.md"
    shell:
        "echo '# Tourmaline Report' > {output};"
        "echo '' >> {output};"
        "echo 'View this HTML report with [Chrome](https://www.google.com/chrome/){{target=\"_blank\"}} or [Firefox](https://www.mozilla.org/en-US/firefox/new/){{target=\"_blank\"}} for best results.' >> {output};"
        "echo '' >> {output};"
        "echo 'To view the linked files below: ' >> {output};"
        "echo '' >> {output};"
        "echo '* QZV (QIIME 2 visualization): click to download, then drag and drop in [https://view.qiime2.org](https://view.qiime2.org){{target=\"_blank\"}}.' >> {output};"
        "echo '* TSV (tab-separated values): click to download, then open in Microsoft Excel or Tabview (command line tool that comes with Tourmaline).' >> {output};"
        "echo '* PDF (portable document format): click to open and view in new tab.' >> {output};"
        "echo '* Markdown and text: click to open and view in new tab.' >> {output};"
        "echo '' >> {output};"
        "echo 'Note: Downloaded files can be deleted after viewing because they are already stored in your Tourmaline directory.' >> {output};"
        "echo '' >> {output};"
        "echo 'For information on Tourmaline outputs, visit [https://github.com/aomlomics/tourmaline](https://github.com/aomlomics/tourmaline){{target=\"_blank\"}}.' >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Metadata Summary' >> {output};"
        "echo '' >> {output};"
        "echo Markdown: \[{input.mdsummary}\]\(../{input.mdsummary}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "cat {input.mdsummary} >> {output};"
        "echo '' >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Fastq Sequences Information' >> {output};"
        "echo '' >> {output};"
        "echo '### Fastq Sequences per Sample' >> {output};"
        "echo '' >> {output};"
        "echo Markdown: \[{input.fastq}\]\(../{input.fastq}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "cat {input.fastq} >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Fastq Sequences' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.visfastq}\]\(../{input.visfastq}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Representative Sequences Information' >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Properties Summary' >> {output};"
        "echo '' >> {output};"
        "echo Markdown: \[{input.repseqsdescribe}\]\(../{input.repseqsdescribe}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "cat {input.repseqsdescribe} >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Properties Table' >> {output};"
        "echo '' >> {output};"
        "echo TSV: \[{input.repseqstsv}\]\(../{input.repseqstsv}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo 'Columns:' >> {output};"
        "echo '' >> {output};"
        "echo '* featureid' >> {output};"
        "echo '* length - length (bp) not including gaps' >> {output};"
        "echo '* gaps - gaps (bp) in multiple sequence alignment' >> {output};"
        "echo '* outlier - outlier (True/False) determined by OD-seq' >> {output};"
        "echo '* taxonomy - domain level' >> {output};"
        "echo '* observations - total observations summed across all samples (unrarefied)' >> {output};"
        "echo '* log10(observations) - log base-10 of total observations' >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Properties Plot' >> {output};"
        "echo '' >> {output};"
        "echo PDF: \[{input.repseqspdf}\]\(../{input.repseqspdf}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo 'Plot elements:' >> {output};"
        "echo '' >> {output};"
        "echo '* x: length (bp) not including gaps' >> {output};"
        "echo '* y: gaps (bp) in multiple sequence alignment' >> {output};"
        "echo '* color: taxonomy (domain)' >> {output};"
        "echo '* size: log10(observations)' >> {output};"
        "echo '* facets: outlier (True/False) determined by OD-seq' >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Rooted Tree' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.repseqstree}\]\(../{input.repseqstree}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences to Filter' >> {output};"
        "echo '' >> {output};"
        "echo 'To filter, go to 02-output-{{method}}-{{filter}}/02-alignment-tree, merge or copy repseqs_to_filter_outliers.tsv and/or repseqs_to_filter_unassigned.tsv, rename to 00-data/repseqs_to_filter_{{method}}.tsv, then run Snakemake in filtered mode.' >> {output};"
        "echo '' >> {output};"
        "echo Outliers TSV: \[{input.repseqsoutliers}\]\(../{input.repseqsoutliers}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Unassigned TSV: \[{input.repseqsunassigned}\]\(../{input.repseqsunassigned}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Representative Sequences' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.visrepseqs}\]\(../{input.visrepseqs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Taxonomy Table' >> {output};"
        "echo '' >> {output};"
        "echo TSV: \[{input.tsvtaxonomy}\]\(../{input.tsvtaxonomy}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Taxonomy' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.vistaxonomy}\]\(../{input.vistaxonomy}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Predicted Amplicon Type' >> {output};"
        "echo '' >> {output};"
        "echo Text: \[{input.amplicontype}\]\(../{input.amplicontype}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "head -n 5 {input.amplicontype} >> {output};"
        "echo '```' >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Observation Table Information' >> {output};"
        "echo '' >> {output};"
        "echo '### Table Summary: Samples' >> {output};"
        "echo '' >> {output};"
        "echo Text: \[{input.samples}\]\(../{input.samples}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "head -n 13 {input.samples} >> {output};"
        "echo '```' >> {output};"
        "echo '' >> {output};"
        "echo '### Table Summary: Features' >> {output};"
        "echo '' >> {output};"
        "echo Text: \[{input.features}\]\(../{input.features}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "head -n 13 {input.features} >> {output};"
        "echo '```' >> {output};"
        "echo '' >> {output};"
        "echo '### Visualization of Table Summary' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.vistable}\]\(../{input.vistable}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Taxonomic Diversity Results' >> {output};"
        "echo '' >> {output};"
        "echo '### Taxonomy Barplot' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.vistaxbar}\]\(../{input.vistaxbar}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Alpha Diversity Results' >> {output};"
        "echo '' >> {output};"
        "echo '### Alpha Rarefaction' >> {output};"
        "echo '' >> {output};"
        "echo QZV: \[{input.visalpharare}\]\(../{input.visalpharare}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Alpha Group Significance' >> {output};"
        "echo '' >> {output};"
        "echo Evenness QZV: \[{input.visevengs}\]\(../{input.visevengs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Faith PD QZV: \[{input.visfaithgs}\]\(../{input.visfaithgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Observed features QZV: \[{input.visobsfeaturesgs}\]\(../{input.visobsfeaturesgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Shannon QZV: \[{input.visshannongs}\]\(../{input.visshannongs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Beta Diversity Results' >> {output};"
        "echo '' >> {output};"
        "echo '### PCoA Emperor Plots' >> {output};"
        "echo '' >> {output};"
        "echo Bray-Curtis QZV: \[{input.visbcemp}\]\(../{input.visbcemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Jaccard QZV: \[{input.visjacemp}\]\(../{input.visjacemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Weighted UniFrac QZV: \[{input.viswuemp}\]\(../{input.viswuemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Unweighted UniFrac QZV: \[{input.visuwuemp}\]\(../{input.visuwuemp}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '### Beta Group Significance' >> {output};"
        "echo '' >> {output};"
        "echo Bray-Curtis QZV: \[{input.visbcgs}\]\(../{input.visbcgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Jaccard QZV: \[{input.visjacgs}\]\(../{input.visjacgs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Weighted UniFrac QZV: \[{input.viswugs}\]\(../{input.viswugs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo Unweighted UniFrac QZV: \[{input.visuwugs}\]\(../{input.visuwugs}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '---' >> {output};"
        "echo '' >> {output};"
        "echo '## Tourmaline Config File' >> {output};"
        "echo '' >> {output};"
        "echo YAML: \[{input.configfile}\]\(../{input.configfile}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "echo '```' >> {output};"
        "cat {input.configfile} >> {output};"
        "echo '```' >> {output};"

rule generate_report_html:
    input:
        markdown="03-reports/report_{method}_{filter}.md"
    params:
        theme=config["report_theme"]
    output:
        "03-reports/report_{method}_{filter}.html"
    shell:
        "pandoc -i {input} -o {output};"
        "echo '<!DOCTYPE html>' > header.html;"
        "echo '<html>' >> header.html;"
        "echo '<head>' >> header.html;"
        "echo '<link rel=\"stylesheet\" type=\"text/css\" href=\"../css/{params.theme}.css\">' >> header.html;"
        "echo '</head>' >> header.html;"
        "echo '<body>' >> header.html;"
        "echo '' >> header.html;"
        "cat header.html {output} > temp.html;"
        "echo '' >> temp.html;"
        "echo '</body>' >> temp.html;"
        "echo '</html>' >> temp.html;"
        "mv temp.html {output};"
        "/bin/rm header.html"
