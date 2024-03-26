# Requires qiime2-2023.5 conda environment 

# GLOBALS ----------------------------------------------------------------------

configfile: "config.yaml"

# RULEORDER DIRECTIVES ---------------------------------------------------------

ruleorder: check_inputs_params_pe > import_fastq_demux_pe
ruleorder: check_inputs_params_se > import_fastq_demux_se
ruleorder: summarize_fastq_demux_pe > summarize_fastq_demux_se
ruleorder: feature_classifier > import_taxonomy_to_qza
ruleorder: filter_taxonomy > import_taxonomy_to_qza > export_taxonomy_to_tsv

# PSEUDO-RULES: DADA2 PAIRED-END -----------------------------------------------

rule dada2_pe_denoise:
    input:
        "01-imported/check_inputs_params_pe.done",
        "01-imported/fastq_summary.qzv",
        "01-imported/fastq_counts_describe.md",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/table_summary.qzv",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs_lengths_describe.md"

rule dada2_pe_taxonomy_unfiltered:
    input:
        "01-imported/check_inputs_params_pe.done",
        "02-output-dada2-pe-unfiltered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-pe-unfiltered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-pe-unfiltered/01-taxonomy/taxa_barplot.qzv",
        "02-output-dada2-pe-unfiltered/01-taxonomy/taxa_sample_table.tsv",
        "02-output-dada2-pe-unfiltered/01-taxonomy/asv_taxa_sample_table.tsv"

rule dada2_pe_diversity_unfiltered:
    input:
        "01-imported/check_inputs_params_pe.done",
        "02-output-dada2-pe-unfiltered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_properties_describe.md",
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
        "02-output-dada2-pe-unfiltered/04-beta-diversity/weighted_unifrac_group_significance.qzv",
        "02-output-dada2-pe-unfiltered/04-beta-diversity/deicode_biplot_emperor.qzv"

rule dada2_pe_report_unfiltered:
    input:
        "01-imported/check_inputs_params_pe.done",
        "03-reports/report_dada2-pe_unfiltered.html"

rule dada2_pe_taxonomy_filtered:
    input:
        "01-imported/check_inputs_params_pe.done",
        "02-output-dada2-pe-filtered/00-table-repseqs/table_summary.qzv",
        "02-output-dada2-pe-filtered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-pe-filtered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-pe-filtered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-pe-filtered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-pe-filtered/00-table-repseqs/repseqs_lengths_describe.md",
        "02-output-dada2-pe-filtered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-pe-filtered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-pe-filtered/01-taxonomy/taxa_barplot.qzv",
        "02-output-dada2-pe-filtered/01-taxonomy/taxa_sample_table.tsv",
        "02-output-dada2-pe-filtered/01-taxonomy/asv_taxa_sample_table.tsv"

rule dada2_pe_diversity_filtered:
    input:
        "01-imported/check_inputs_params_pe.done",
        "02-output-dada2-pe-filtered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-pe-filtered/02-alignment-tree/repseqs_properties_describe.md",
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
        "02-output-dada2-pe-filtered/04-beta-diversity/weighted_unifrac_group_significance.qzv",
        "02-output-dada2-pe-filtered/04-beta-diversity/deicode_biplot_emperor.qzv"

rule dada2_pe_report_filtered:
    input:
        "01-imported/check_inputs_params_pe.done",
        "03-reports/report_dada2-pe_filtered.html"

# PSEUDO-RULES: DADA2 SINGLE-END -----------------------------------------------

rule dada2_se_denoise:
    input:
        "01-imported/check_inputs_params_se.done",
        "01-imported/fastq_summary.qzv",
        "01-imported/fastq_counts_describe.md",
        "02-output-dada2-se-unfiltered/00-table-repseqs/table_summary.qzv",
        "02-output-dada2-se-unfiltered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-se-unfiltered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-se-unfiltered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-se-unfiltered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-se-unfiltered/00-table-repseqs/repseqs_lengths_describe.md"

rule dada2_se_taxonomy_unfiltered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-dada2-se-unfiltered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-se-unfiltered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-se-unfiltered/01-taxonomy/taxa_barplot.qzv",
        "02-output-dada2-se-unfiltered/01-taxonomy/taxa_sample_table.tsv",
        "02-output-dada2-se-unfiltered/01-taxonomy/asv_taxa_sample_table.tsv"

rule dada2_se_diversity_unfiltered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-dada2-se-unfiltered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-se-unfiltered/02-alignment-tree/repseqs_properties_describe.md",
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
        "02-output-dada2-se-unfiltered/04-beta-diversity/weighted_unifrac_group_significance.qzv",
        "02-output-dada2-se-unfiltered/04-beta-diversity/deicode_biplot_emperor.qzv"

rule dada2_se_report_unfiltered:
    input:
        "01-imported/check_inputs_params_se.done",
        "03-reports/report_dada2-se_unfiltered.html"

rule dada2_se_taxonomy_filtered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-dada2-se-filtered/00-table-repseqs/table_summary.qzv",
        "02-output-dada2-se-filtered/00-table-repseqs/table_summary_samples.txt",
        "02-output-dada2-se-filtered/00-table-repseqs/table_summary_features.txt",
        "02-output-dada2-se-filtered/00-table-repseqs/repseqs.qzv",
        "02-output-dada2-se-filtered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-dada2-se-filtered/00-table-repseqs/repseqs_lengths_describe.md",
        "02-output-dada2-se-filtered/01-taxonomy/taxonomy.tsv",
        "02-output-dada2-se-filtered/01-taxonomy/taxonomy.qzv",
        "02-output-dada2-se-filtered/01-taxonomy/taxa_barplot.qzv",
        "02-output-dada2-se-filtered/01-taxonomy/taxa_sample_table.tsv",
        "02-output-dada2-se-filtered/01-taxonomy/asv_taxa_sample_table.tsv"

rule dada2_se_diversity_filtered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-dada2-se-filtered/02-alignment-tree/rooted_tree.qzv",
        "02-output-dada2-se-filtered/02-alignment-tree/repseqs_properties_describe.md",
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
        "02-output-dada2-se-filtered/04-beta-diversity/weighted_unifrac_group_significance.qzv",
        "02-output-dada2-se-filtered/04-beta-diversity/deicode_biplot_emperor.qzv"

rule dada2_se_report_filtered:
    input:
        "01-imported/check_inputs_params_se.done",
        "03-reports/report_dada2-se_filtered.html"

# PSEUDO-RULES: DEBLUR SINGLE-END ----------------------------------------------

rule deblur_se_denoise:
    input:
        "01-imported/check_inputs_params_se.done",
        "01-imported/fastq_summary.qzv",
        "01-imported/fastq_counts_describe.md",
        "02-output-deblur-se-unfiltered/00-table-repseqs/table_summary.qzv",
        "02-output-deblur-se-unfiltered/00-table-repseqs/table_summary_samples.txt",
        "02-output-deblur-se-unfiltered/00-table-repseqs/table_summary_features.txt",
        "02-output-deblur-se-unfiltered/00-table-repseqs/repseqs.qzv",
        "02-output-deblur-se-unfiltered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-deblur-se-unfiltered/00-table-repseqs/repseqs_lengths_describe.md"

rule deblur_se_taxonomy_unfiltered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-deblur-se-unfiltered/01-taxonomy/taxonomy.tsv",
        "02-output-deblur-se-unfiltered/01-taxonomy/taxonomy.qzv",
        "02-output-deblur-se-unfiltered/01-taxonomy/taxa_barplot.qzv",
        "02-output-deblur-se-unfiltered/01-taxonomy/taxa_sample_table.tsv",
        "02-output-deblur-se-unfiltered/01-taxonomy/asv_taxa_sample_table.tsv"

rule deblur_se_diversity_unfiltered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-deblur-se-unfiltered/02-alignment-tree/rooted_tree.qzv",
        "02-output-deblur-se-unfiltered/02-alignment-tree/repseqs_properties_describe.md",
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
        "02-output-deblur-se-unfiltered/04-beta-diversity/weighted_unifrac_group_significance.qzv",
        "02-output-deblur-se-unfiltered/04-beta-diversity/deicode_biplot_emperor.qzv"

rule deblur_se_report_unfiltered:
    input:
        "01-imported/check_inputs_params_se.done",
        "03-reports/report_deblur-se_unfiltered.html"

rule deblur_se_taxonomy_filtered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-deblur-se-filtered/00-table-repseqs/table_summary.qzv",
        "02-output-deblur-se-filtered/00-table-repseqs/table_summary_samples.txt",
        "02-output-deblur-se-filtered/00-table-repseqs/table_summary_features.txt",
        "02-output-deblur-se-filtered/00-table-repseqs/repseqs.qzv",
        "02-output-deblur-se-filtered/00-table-repseqs/repseqs_amplicon_type.txt",
        "02-output-deblur-se-filtered/00-table-repseqs/repseqs_lengths_describe.md",
        "02-output-deblur-se-filtered/01-taxonomy/taxonomy.tsv",
        "02-output-deblur-se-filtered/01-taxonomy/taxonomy.qzv",
        "02-output-deblur-se-filtered/01-taxonomy/taxa_barplot.qzv",
        "02-output-deblur-se-filtered/01-taxonomy/taxa_sample_table.tsv",
        "02-output-deblur-se-filtered/01-taxonomy/asv_taxa_sample_table.tsv"

rule deblur_se_diversity_filtered:
    input:
        "01-imported/check_inputs_params_se.done",
        "02-output-deblur-se-filtered/02-alignment-tree/rooted_tree.qzv",
        "02-output-deblur-se-filtered/02-alignment-tree/repseqs_properties_describe.md",
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
        "02-output-deblur-se-filtered/04-beta-diversity/weighted_unifrac_group_significance.qzv",
        "02-output-deblur-se-filtered/04-beta-diversity/deicode_biplot_emperor.qzv"

rule deblur_se_report_filtered:
    input:
        "01-imported/check_inputs_params_se.done",
        "03-reports/report_deblur-se_filtered.html"

# RULES: CHECKS ----------------------------------------------------------------

rule check_metadata:
    output:
        touch("01-imported/check_metadata.done")
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "if [ -r '00-data/metadata.tsv' ]; then "
        "    echo 'OK: Metadata file 00-data/metadata.tsv found.'; "
        "else echo 'Error: Metadata file 00-data/metadata.tsv not found.' && exit 1; "
        "fi; "

rule summarize_metadata:
    input:
        metadata="00-data/metadata.tsv",
        check="01-imported/check_metadata.done"
    output:
        "01-imported/metadata_summary.md",
        "01-imported/metadata_columns.txt"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "python scripts/summarize_metadata.py {input.metadata} {output[0]} {output[1]}"

rule check_inputs_params_pe:
    input:
        "01-imported/metadata_columns.txt",
        check="01-imported/check_metadata.done"
    params:
        column=config["beta_group_column"],
        classifymethod=config["classify_method"],
    output:
        touch("01-imported/check_inputs_params_pe.done")
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "if [ -r '01-imported/fastq_pe.qza' ]; then "
        "    echo 'OK: FASTQ archive 01-imported/fastq_pe.qza found; FASTQ manifest file 00-data/manifest_pe.csv not required.'; "
        "elif [ -r '00-data/manifest_pe.csv' ]; then "
        "    echo 'OK: FASTQ manifest file 00-data/manifest_pe.csv found; it will be used to create FASTQ archive 01-imported/fastq_pe.qza.'; "
        "else echo 'Error: FASTQ sequence data not found; either 00-data/manifest_pe.csv or 01-imported/fastq_pe.qza is required.' && exit 1; "
        "fi; "
        "if [ {params.classifymethod} = naive-bayes ]; then "
        "    if [ -r '01-imported/classifier.qza' ]; then "
        "        echo 'OK: Reference sequences classifier 01-imported/classifier.qza found; reference sequences FASTA file 00-data/refseqs.fna not required.'; "
        "    elif [ -r '01-imported/refseqs.qza' ]; then "
        "        echo 'OK: Reference sequences archive 01-imported/refseqs.qza found; reference sequences FASTA file 00-data/refseqs.fna not required.'; "
        "    elif [ -r '00-data/refseqs.fna' ]; then "
        "        echo 'OK: Reference sequences FASTA file 00-data/refseqs.fna found; it will be used to create reference sequences archive 01-imported/refseqs.qza.'; "
        "        qiime tools import --type 'FeatureData[Sequence]' --input-path 00-data/refseqs.fna --output-path 01-imported/refseqs.qza; "
        "        if [ -r '00-data/reftax.tsv' ]; then "
        "             qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 00-data/reftax.tsv --output-path 01-imported/reftax.qza; "
        "    else echo 'Error: Reference sequences not found; either 01-imported/classifier.qza or 00-data/refseqs.fna or 01-imported/refseqs.qza is required.' && exit 1; "
        "    fi; "
        "elif [ -r '01-imported/refseqs.qza' ]; then "
        "    echo 'OK: Reference sequences archive 01-imported/refseqs.qza found; reference sequences FASTA file 00-data/refseqs.fna not required.'; "
        "elif [ -r '00-data/refseqs.fna' ]; then "
        "    echo 'OK: Reference sequences FASTA file 00-data/refseqs.fna found; it will be used to create reference sequences archive 01-imported/refseqs.qza.'; "
        "else echo 'Error: Reference sequences not found; either 00-data/refseqs.fna or 01-imported/refseqs.qza is required.' && exit 1; "
        "fi; "
        "if [ {params.classifymethod} != naive-bayes ]; then "
        "    if [ -r '01-imported/reftax.qza' ]; then "
        "        echo 'OK: Reference taxonomy archive 01-imported/reftax.qza found; reference taxonomy file 00-data/reftax.tsv not required.'; "
        "    elif [ -r '00-data/reftax.tsv' ]; then "
        "        echo 'OK: Reference taxonomy file 00-data/reftax.tsv found; it will be used to create reference taxonomy archive 01-imported/reftax.qza.'; "
        "        qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 00-data/reftax.tsv --output-path 01-imported/reftax.qza;"
        "    else echo 'Error: Reference taxonomy not found; either 00-data/reftax.tsv or 01-imported/reftax.qza is required.' && exit 1; "
        "    fi; "
        "fi; "
        "if grep -q ^{params.column}$ {input}; then "
        "    echo 'OK: Metadata contains the column \"{params.column}\" that is specified as beta_group_column in config.yaml.'; "
        "else echo 'Error: Metadata does not contain the column \"{params.column}\" that is specified as beta_group_column in config.yaml.' && exit 1; "
        "fi"

rule check_inputs_params_se:
    input:
        "01-imported/metadata_columns.txt"
    params:
        column=config["beta_group_column"],
        classifymethod=config["classify_method"],
    output:
        touch("01-imported/check_inputs_params_se.done")
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "if [ -r '01-imported/fastq_se.qza' ]; then "
        "    echo 'OK: FASTQ archive 01-imported/fastq_se.qza found; FASTQ manifest file 00-data/manifest_se.csv not required.'; "
        "elif [ -r '00-data/manifest_se.csv' ]; then "
        "    echo 'OK: FASTQ manifest file 00-data/manifest_se.csv found; it will be used to create FASTQ archive 01-imported/fastq_se.qza.'; "
        "else echo 'Error: FASTQ sequence data not found; either 00-data/manifest_se.csv or 01-imported/fastq_se.qza is required.' && exit 1; "
        "fi; "
        "if [ {params.classifymethod} = naive-bayes ]; then "
        "    if [ -r '01-imported/classifier.qza' ]; then "
        "        echo 'OK: Reference sequences classifier 01-imported/classifier.qza found; reference sequences FASTA file 00-data/refseqs.fna not required.'; "
        "    elif [ -r '01-imported/refseqs.qza' ]; then "
        "        echo 'OK: Reference sequences archive 01-imported/refseqs.qza found; reference sequences FASTA file 00-data/refseqs.fna not required.'; "
        "    elif [ -r '00-data/refseqs.fna' ]; then "
        "        echo 'OK: Reference sequences FASTA file 00-data/refseqs.fna found; it will be used to create reference sequences archive 01-imported/refseqs.qza.'; "
        "        qiime tools import --type 'FeatureData[Sequence]' --input-path 00-data/refseqs.fna --output-path 01-imported/refseqs.qza; "
        "        if [ -r '00-data/reftax.tsv' ]; then "
        "             qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 00-data/reftax.tsv --output-path 01-imported/reftax.qza; "
        "    else echo 'Error: Reference sequences not found; either 01-imported/classifier.qza or 00-data/refseqs.fna or 01-imported/refseqs.qza is required.' && exit 1; "
        "    fi; "
        "elif [ -r '01-imported/refseqs.qza' ]; then "
        "    echo 'OK: Reference sequences archive 01-imported/refseqs.qza found; reference sequences FASTA file 00-data/refseqs.fna not required.'; "
        "elif [ -r '00-data/refseqs.fna' ]; then "
        "    echo 'OK: Reference sequences FASTA file 00-data/refseqs.fna found; it will be used to create reference sequences archive 01-imported/refseqs.qza.'; "
        "    qiime tools import --type 'FeatureData[Sequence]' --input-path 00-data/refseqs.fna --output-path 01-imported/refseqs.qza; "
        "else echo 'Error: Reference sequences not found; either 00-data/refseqs.fna or 01-imported/refseqs.qza is required.' && exit 1; "
        "fi; "
        "if [ {params.classifymethod} != naive-bayes ]; then "
        "    if [ -r '01-imported/reftax.qza' ]; then "
        "        echo 'OK: Reference taxonomy archive 01-imported/reftax.qza found; reference taxonomy file 00-data/reftax.tsv not required.'; "
        "    elif [ -r '00-data/reftax.tsv' ]; then "
        "        echo 'OK: Reference taxonomy file 00-data/reftax.tsv found; it will be used to create reference taxonomy archive 01-imported/reftax.qza.'; "
        "        qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path 00-data/reftax.tsv --output-path 01-imported/reftax.qza;"
        "    else echo 'Error: Reference taxonomy not found; either 00-data/reftax.tsv or 01-imported/reftax.qza is required.' && exit 1; "
        "    fi; "
        "fi; "
        "if grep -q ^{params.column}$ {input}; then "
        "    echo 'OK: Metadata contains the column \"{params.column}\" that is specified as beta_group_column in config.yaml.'; "
        "else echo 'Error: Metadata does not contain the column \"{params.column}\" that is specified as beta_group_column in config.yaml.' && exit 1; "
        "fi"

# RULES: IMPORT ----------------------------------------------------------------

rule import_ref_seqs:
    input:
        "00-data/refseqs.fna"
    output:
        "01-imported/refseqs.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools import "
        "--type 'FeatureData[Sequence]' "
        "--input-path {input} "
        "--output-path {output}"

rule import_ref_tax:
    input:
        "00-data/reftax.tsv"
    output:
        "01-imported/reftax.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools import "
        "--type 'FeatureData[Taxonomy]' "
        "--input-format HeaderlessTSVTaxonomyFormat "
        "--input-path {input} "
        "--output-path {output}"

rule import_fastq_demux_pe:
    input:
        "00-data/manifest_pe.csv",
        "01-imported/check_inputs_params_pe.done"
    output:
        "01-imported/fastq_pe.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input[0]} "
        "--output-path {output} "
        "--input-format PairedEndFastqManifestPhred33"

rule import_fastq_demux_se:
    input:
        "00-data/manifest_se.csv",
        "01-imported/check_inputs_params_se.done"
    output:
        "01-imported/fastq_se.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools import "
        "--type 'SampleData[SequencesWithQuality]' "
        "--input-path {input[0]} "
        "--output-path {output} "
        "--input-format SingleEndFastqManifestPhred33"

# RULES: SUMMARIZE FASTQ SEQUENCES ---------------------------------------------

rule summarize_fastq_demux_pe:
    input:
        "01-imported/fastq_pe.qza",
        "01-imported/check_inputs_params_pe.done"
    output:
        "01-imported/fastq_summary.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime demux summarize "
        "--i-data {input[0]} "
        "--o-visualization {output}"

rule summarize_fastq_demux_se:
    input:
        "01-imported/fastq_se.qza",
        "01-imported/check_inputs_params_se.done"
    output:
        "01-imported/fastq_summary.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime demux summarize "
        "--i-data {input[0]} "
        "--o-visualization {output}"

rule export_fastq_summary_to_counts:
    input:
        "01-imported/fastq_summary.qzv"
    output:
        "01-imported/fastq_counts.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "unzip -qq -o {input} -d temp0; "
        "mv temp0/*/data/per-sample-fastq-counts.tsv {output}; "
        "/bin/rm -r temp0"

rule describe_fastq_counts:
    input:
        "01-imported/fastq_counts.tsv"
    output:
        "01-imported/fastq_counts_describe.md"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "python scripts/describe_fastq_counts.py {input} {output}"

# RULES: DENOISE ---------------------------------------------------------------

rule denoise_dada2_pe:
    input:
        "01-imported/fastq_pe.qza",
        "01-imported/check_inputs_params_pe.done"
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
    conda:
        "qiime2-2023.5"
    threads: config["dada2pe_threads"]
    shell:
        "qiime dada2 denoise-paired "
        "--i-demultiplexed-seqs {input[0]} "
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
        "01-imported/fastq_se.qza",
        "01-imported/check_inputs_params_se.done"
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
    conda:
        "qiime2-2023.5"
    conda:
        "qiime2-2023.5"
    threads: config["dada2se_threads"]
    shell:
        "qiime dada2 denoise-single "
        "--i-demultiplexed-seqs {input[0]} "
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
        seqs="01-imported/fastq_se.qza",
        refseqs="01-imported/refseqs.qza",
        check="01-imported/check_inputs_params_se.done"
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
    conda:
        "qiime2-2023.5"
    conda:
        "qiime2-2023.5"
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
        metadata="00-data/metadata.tsv"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/table_summary.qzv"
    conda:
        "qiime2-2023.5"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime feature-table summarize "
        "--i-table {input.table} "
        "--m-sample-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule summarize_repseqs:
    input:
        stats="02-output-{method}-{filter}/00-table-repseqs/dada2_stats.qza"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/dada2_stats.qzv"
    conda:
        "qiime2-2023.5"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime metadata tabulate "
        "--m-input-file {input.stats} "
        "--o-visualization {output}"


rule export_table_to_biom:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/table.qza"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/table.biom"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "--output-format BIOMV210Format"

rule summarize_biom_samples:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/table.biom"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/table_summary_samples.txt"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
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
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
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
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime feature-table tabulate-seqs "
        "--i-data {input} "
        "--o-visualization {output}"

rule export_repseqs_to_fasta:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.qza"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "--output-format DNAFASTAFormat"

rule repseqs_detect_amplicon_locus:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_amplicon_type.txt"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "python scripts/detect_amplicon_locus.py -i {input} > {output}"

rule repseqs_lengths:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "perl scripts/fastaLengths.pl {input} > {output}"

rule repseqs_lengths_describe:
    input:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths.tsv"
    output:
        "02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths_describe.md"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "python scripts/repseqs_lengths_describe.py {input} {output}"

# RULES: TAXONOMY --------------------------------------------------------------

rule feature_classifier:
    input:
        repseqs="02-output-{method}-unfiltered/00-table-repseqs/repseqs.qza"
    output:
        "02-output-{method}-unfiltered/01-taxonomy/taxonomy.qza"
    params:
        refseqs="01-imported/refseqs.qza",
        reftax="01-imported/reftax.qza",
        classifymethod=config["classify_method"],
        classifyparams=config["classify_parameters"],
        searchout="02-output-{method}-unfiltered/01-taxonomy/search_results.qza"
    conda:
        "qiime2-2023.5"
    threads: config["feature_classifier_threads"]
    shell:
        "echo classify_method: {params.classifymethod}; "
        "if [ {params.classifymethod} = naive-bayes ]; then "
        "    if [ ! -r 01-imported/classifier.qza ]; then "
        "        qiime feature-classifier fit-classifier-naive-bayes "
        "        --i-reference-reads {params.refseqs} "
        "        --i-reference-taxonomy {params.reftax} "
        "        --o-classifier 01-imported/classifier.qza; "
        "    fi; "
        "    qiime feature-classifier classify-sklearn "
        "    --i-classifier 01-imported/classifier.qza "
        "    --i-reads {input.repseqs} "
        "    --o-classification {output} "
        "    --p-n-jobs {threads} "
        "    {params.classifyparams}; "
        "elif [ {params.classifymethod} = consensus-blast ]; then "
        "    qiime feature-classifier classify-consensus-blast "
        "    --i-reference-reads {params.refseqs} "
        "    --i-reference-taxonomy {params.reftax} "
        "    --i-query {input.repseqs} "
        "    --o-classification {output} "
        "    --o-search-results {params.searchout} "
        "    {params.classifyparams}; "
        "elif [ {params.classifymethod} = consensus-vsearch ]; then "
        "    qiime feature-classifier classify-consensus-vsearch "
        "    --i-reference-reads {params.refseqs} "
        "    --i-reference-taxonomy {params.reftax} "
        "    --i-query {input.repseqs} "
        "    --o-classification {output} "
        "    --o-search-results {params.searchout} "
        "    --p-threads {threads} "
        "    {params.classifyparams}; "
        "fi"

rule visualize_taxonomy:
    input:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qza"
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime metadata tabulate "
        "--m-input-file {input} "
        "--o-visualization {output}"

rule taxa_barplot:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
        metadata="00-data/metadata.tsv"
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxa_barplot.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime taxa barplot "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--m-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule export_taxa_biom:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxa_sample_table.tsv"
    params:
        taxalevel=config["classify_taxalevel"]
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime taxa collapse "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--p-level {params.taxalevel} "
        "--o-collapsed-table tempfile_collapsed.qza;"
        "qiime tools export "
        "--input-path tempfile_collapsed.qza "
        "--output-path temp_export;"
        "biom convert "
        "-i temp_export/feature-table.biom "
        "-o {output} "
        "--to-tsv;"
        "/bin/rm -r tempfile_collapsed.qza temp_export/"

rule export_asv_seq_taxa_obis:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
        repseqs="02-output-{method}-{filter}/00-table-repseqs/repseqs.qza"
    output:
        "02-output-{method}-{filter}/01-taxonomy/asv_taxa_sample_table.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime feature-table transpose "
        "--i-table {input.table} "
        "--o-transposed-feature-table transposed-table.qza; "
        "qiime metadata tabulate "
        "--m-input-file {input.repseqs} "
        "--m-input-file {input.taxonomy} "
        "--m-input-file transposed-table.qza "
        "--o-visualization merged-data.qzv; "
        "qiime tools export "
        "--input-path merged-data.qzv "
        "--output-path {output}; "
        "mv {output}/metadata.tsv temp; "
        "rm -r {output}; "
        "sed -e '2d' temp | sed '1 s|id\\t|featureid\\t|' | sed '1 s|Taxon|taxonomy|' | sed '1 s|Sequence|sequence|' > {output}; "
        "/bin/rm -r temp transposed-table.qza merged-data.qzv"

rule export_taxonomy_to_tsv:
    input:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qza"
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "--output-format TSVTaxonomyFormat"

rule import_taxonomy_to_qza:
    input:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.tsv"
    output:
        "02-output-{method}-{filter}/01-taxonomy/taxonomy.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime tools import "
        "--type 'FeatureData[Taxonomy]' "
        "--input-format TSVTaxonomyFormat "
        "--input-path {input} "
        "--output-path {output}"

# RULES: ALIGNMENT & TREE ------------------------------------------------------

rule align_repseqs:
    input:
        repseqsfasta="02-output-{method}-{filter}/00-table-repseqs/repseqs.fasta",
        repseqsqza="02-output-{method}-{filter}/00-table-repseqs/repseqs.qza"
    output:
        alnfasta="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta",
        alnqza="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza"
    params:
        method=config["alignment_method"],
        muscle_iters=config["alignment_muscle_iters"]
    conda:
        "qiime2-2023.5"
    threads: config["alignment_threads"],
    shell:
        "if [ {params.method} = muscle ]; then "
        "    echo 'Multiple sequence alignment method: MUSCLE ...'; "
        "    muscle "
        "    -super5 {input.repseqsfasta} "
        "    -threads {threads} "
        "    -refineiters {params.muscle_iters} "
        "    -output temp_aligned_repseqs.fasta; "
        "    perl scripts/cleanupMultiFastaNoBreaks.pl temp_aligned_repseqs.fasta > {output.alnfasta}; "
        "    echo 'Line breaks removed to generate {output.alnfasta}'; "
        "    /bin/rm temp_aligned_repseqs.fasta; "
        "    qiime tools import "
        "    --type 'FeatureData[AlignedSequence]' "
        "    --input-path {output.alnfasta} "
        "    --output-path {output.alnqza}; "
        "elif [ {params.method} = clustalo ]; then "
        "    echo 'Multiple sequence alignment method: Clustal Omega ...'; "
        "    clustalo --verbose --force "
        "    --in {input.repseqsfasta} "
        "    --out temp_aligned_repseqs.fasta "
        "    --threads={threads}; "
        "    perl scripts/cleanupMultiFastaNoBreaks.pl temp_aligned_repseqs.fasta > {output.alnfasta}; "
        "    echo 'Line breaks removed to generate {output.alnfasta}'; "
        "    /bin/rm temp_aligned_repseqs.fasta; "
        "    qiime tools import "
        "    --type 'FeatureData[AlignedSequence]' "
        "    --input-path {output.alnfasta} "
        "    --output-path {output.alnqza}; "
        "elif [ {params.method} = mafft ]; then "
        "    echo 'Multiple sequence alignment method: MAFFT ...'; "
        "    qiime alignment mafft "
        "    --i-sequences {input.repseqsqza} "
        "    --o-alignment tempfile_unmasked_aligned_repseqs.qza "
        "    --p-n-threads {threads}; "
        "    qiime alignment mask "
        "    --i-alignment tempfile_unmasked_aligned_repseqs.qza "
        "    --o-masked-alignment {output.alnqza}; "
        "    /bin/rm tempfile_unmasked_aligned_repseqs.qza; "
        "    qiime tools export "
        "    --input-path {output.alnqza} "
        "    --output-path {output.alnfasta} "
        "    --output-format AlignedDNAFASTAFormat; "
        "else "
        "    echo 'Multiple sequence alignment method: MUSCLE ...'; "
        "fi"

rule phylogeny_fasttree:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.qza"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/unrooted_tree.qza"
    conda:
        "qiime2-2023.5"
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
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime phylogeny midpoint-root "
        "--i-tree {input} "
        "--o-rooted-tree {output}"

rule visualize_tree:
    input:
        tree="02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qza",
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        metadata="00-data/metadata.tsv",
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
        outliers="02-output-{method}-{filter}/02-alignment-tree/outliers.qza"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
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
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "bash scripts/alignment_count_gaps.sh < {input} > {output}"

# rule alignment_gaps_describe:
#     input:
#         "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps.tsv"
#     output:
#         "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps_describe.md"
#     run:
#         gaps = pd.read_csv(input[0], header=None, sep='\t')
#         t = gaps.describe()
#         outstr = tabulate(t.iloc[1:], tablefmt="pipe", headers=['Statistic (n=%s)' % t.iloc[0].values[0].astype(int), 'Alignment gaps per sequence'])
#         with open(output[0], 'w') as target:
#             target.write(outstr)
#             target.write('\n')

rule alignment_detect_outliers:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs.fasta"
    params:
        metric = config["odseq_distance_metric"],
        replicates = config["odseq_bootstrap_replicates"],
        threshold = config["odseq_threshold"]
    output:
        "02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_outliers.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "Rscript --vanilla scripts/run_odseq.R {input} {params.metric} {params.replicates} {params.threshold} temp_odseq; "
        "cat temp_odseq | sed 's/^X//' > {output}; "
        "/bin/rm temp_odseq"

rule tabulate_plot_repseq_properties:
    input:
        lengths="02-output-{method}-{filter}/00-table-repseqs/repseqs_lengths.tsv",
        gaps="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_gaps.tsv",
        outliers="02-output-{method}-{filter}/02-alignment-tree/aligned_repseqs_outliers.tsv",
        taxonomy="02-output-{method}-{filter}/01-taxonomy/taxonomy.qza",
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza"
    output:
        proptsv="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties.tsv",
        propdescribe="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties_describe.md",
        proppdf="02-output-{method}-{filter}/02-alignment-tree/repseqs_properties.pdf",
        outliersforqza="02-output-{method}-{filter}/02-alignment-tree/outliers.tsv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "python scripts/plot_repseq_properties.py {input.lengths} {input.gaps} {input.outliers} {input.taxonomy} "
        "{input.table} {output.proptsv} {output.propdescribe} {output.proppdf} {output.outliersforqza}"

rule import_outliers_to_qza:
    input:
        "02-output-{method}-{filter}/02-alignment-tree/outliers.tsv"
    output:
        "02-output-{method}-{filter}/02-alignment-tree/outliers.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
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
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "cat {input.proptsv} | grep -i 'outlier\|true' | cut -f1,4 > {output.outliers}; "
        "cat {input.proptsv} | grep -i 'taxonomy\|unassigned' | cut -f1,5 > {output.unassigned}"

# RULES: FILTER ----------------------------------------------------------------

# NEW RULES

rule filter_sequences_table:
    input:
        table="02-output-{method}-unfiltered/00-table-repseqs/table.qza",
        taxonomy="02-output-{method}-unfiltered/01-taxonomy/taxonomy.qza",
        repseqs="02-output-{method}-unfiltered/00-table-repseqs/repseqs.qza",
        repseqstofilter="00-data/repseqs_to_filter_{method}.tsv",
        samplestofilter="00-data/samples_to_filter_{method}.tsv",
        metadata="00-data/metadata.tsv"
    params:
        excludeterms=config["exclude_terms"],
        minlength=config["repseq_min_length"],
        maxlength=config["repseq_max_length"],
        minabund=config["repseq_min_abundance"],
        minprev=config["repseq_min_prevalence"],
        metadatafilt=config["metadata_filter"]
    output:
        repseqs="02-output-{method}-filtered/00-table-repseqs/repseqs.qza",
        table="02-output-{method}-filtered/00-table-repseqs/table.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        # FILTER SEQUENCES BY TAXONOMY
        "qiime taxa filter-seqs "
        "--i-sequences {input.repseqs} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude {params.excludeterms} "
        "--o-filtered-sequences temp_repseqs1.qza; "
        # FILTER SEQUENCES BY LENGTH
        "qiime feature-table filter-seqs "
        "--i-data temp_repseqs1.qza "
        "--m-metadata-file temp_repseqs1.qza "
        "--p-where 'length(sequence) >= {params.minlength} AND length(sequence) <= {params.maxlength}' "
        "--o-filtered-data temp_repseqs2.qza; "
        # FILTER SEQUENCES BY ID
        "qiime feature-table filter-seqs "
        "--i-data temp_repseqs2.qza "
        "--m-metadata-file {input.repseqstofilter} "
        "--p-exclude-ids "
        "--o-filtered-data temp_repseqs3.qza; "
        # FILTER TABLE BY FEATURE IDS
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--m-metadata-file temp_repseqs3.qza "
        "--o-filtered-table temp_table.qza; "
        # FILTER TABLE BY SAMPLE IDS
        "qiime feature-table filter-samples "
        "--i-table temp_table.qza "
        "--m-metadata-file {input.samplestofilter} "
        "--p-exclude-ids "
        "--o-filtered-table temp_table2.qza; "
        # FILTER TABLE BY SAMPLE METADATA
        "if [ {params.metadatafilt} != none ]; then "
        "    qiime feature-table filter-samples "
        "    --i-table temp_table2.qza "
        "    --m-metadata-file {input.metadata} "
        "    --p-where \"{params.metadatafilt}\" "
        "    --o-filtered-table temp_table3.qza; "
        "else "
        "    cp temp_table2.qza temp_table3.qza; "
        "fi; "
        # FILTER TABLE BY ABUNDANCE & PREVALENCE
        "qiime feature-table filter-features-conditionally "
        "--i-table temp_table3.qza "
        "--p-abundance {params.minabund} "
        "--p-prevalence {params.minprev} "
        "--o-filtered-table {output.table}; "
        # FILTER SEQUENCES USING TABLE
        "qiime feature-table filter-seqs "
        "--i-data temp_repseqs3.qza "
        "--i-table {output.table} "
        "--p-no-exclude-ids "
        "--o-filtered-data {output.repseqs}; "
        # REMOVE TEMP FILES
        "/bin/rm temp_repseqs1.qza; "
        "/bin/rm temp_repseqs2.qza; "
        "/bin/rm temp_repseqs3.qza; "
        "/bin/rm temp_table.qza; "
        "/bin/rm temp_table2.qza; "
        "/bin/rm temp_table3.qza;"

rule filter_taxonomy:
    input:
        taxonomy="02-output-{method}-unfiltered/01-taxonomy/taxonomy.tsv",
        repseqs="02-output-{method}-filtered/00-table-repseqs/repseqs_lengths.tsv"
    output:
        taxonomytsv="02-output-{method}-filtered/01-taxonomy/taxonomy.tsv",
        taxonomyqza="02-output-{method}-filtered/01-taxonomy/taxonomy.qza"        
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "python scripts/filter_taxonomy.py {input.taxonomy} {input.repseqs} {output.taxonomytsv} {output.taxonomyqza}"


# RULES: DIVERSITY -------------------------------------------------------------

rule diversity_alpha_rarefaction:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza",
        phylogeny="02-output-{method}-{filter}/02-alignment-tree/rooted_tree.qza",
        metadata="00-data/metadata.tsv"
    params:
        maxdepth=config["alpha_max_depth"]
    output:
        "02-output-{method}-{filter}/03-alpha-diversity/alpha_rarefaction.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
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
        metadata="00-data/metadata.tsv"
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
    conda:
        "qiime2-2023.5"
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
        metadata="00-data/metadata.tsv"
    output:
        "02-output-{method}-{filter}/03-alpha-diversity/{metric}_group_significance.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime diversity alpha-group-significance "
        "--i-alpha-diversity {input.alphadiversity} "
        "--m-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule diversity_beta_group_significance:
    input:
        distancematrix="02-output-{method}-{filter}/04-beta-diversity/{metric}_distance_matrix.qza",
        metadata="00-data/metadata.tsv"
    params:
        column=config["beta_group_column"],
        method=config["beta_group_method"],
        pairwise=config["beta_group_pairwise"]
    output:
        "02-output-{method}-{filter}/04-beta-diversity/{metric}_group_significance.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime diversity beta-group-significance "
        "--i-distance-matrix {input.distancematrix} "
        "--m-metadata-file {input.metadata} "
        "--m-metadata-column {params.column} "
        "--p-method {params.method} "
        "{params.pairwise} "
        "--o-visualization {output}"

rule deicode_auto_rpca:
    input:
        table="02-output-{method}-{filter}/00-table-repseqs/table.qza"
    params:
        minsamplecount=config["deicode_min_sample_count"],
        minfeaturecount=config["deicode_min_feature_count"],
        minfeaturefrequency=config["deicode_min_feature_frequency"],
        maxiterations=config["deicode_max_iterations"]
    output:
        biplot="02-output-{method}-{filter}/04-beta-diversity/deicode_biplot.qza",
        distancematrix="02-output-{method}-{filter}/04-beta-diversity/deicode_distance_matrix.qza"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime deicode auto-rpca "
        "--i-table {input.table} "
        "--p-min-sample-count {params.minsamplecount} "
        "--p-min-feature-count {params.minfeaturecount} "
        "--p-min-feature-frequency {params.minfeaturefrequency} "
        "--p-max-iterations {params.maxiterations} "
        "--o-biplot {output.biplot} "
        "--o-distance-matrix {output.distancematrix}"

rule emperor_biplot:
    input:
        biplot="02-output-{method}-{filter}/04-beta-diversity/deicode_biplot.qza",
        metadata="00-data/metadata.tsv"
    params:
        numfeatures=config["deicode_num_features"]
    output:
        emperor="02-output-{method}-{filter}/04-beta-diversity/deicode_biplot_emperor.qzv"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "qiime emperor biplot "
        "--i-biplot {input.biplot} "
        "--m-sample-metadata-file {input.metadata} "
        "--o-visualization {output.emperor} "
        "--p-number-of-features {params.numfeatures}"

# RULES: REPORT ----------------------------------------------------------------

rule generate_report_md:
    input:
        configfile="config.yaml",
        mdsummary="01-imported/metadata_summary.md",
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
        vistable="02-output-{method}-{filter}/00-table-repseqs/table_summary.qzv",
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
        visuwugs="02-output-{method}-{filter}/04-beta-diversity/unweighted_unifrac_group_significance.qzv",
        visbiplotemp="02-output-{method}-{filter}/04-beta-diversity/deicode_biplot_emperor.qzv"
    params:
        refdatabase=config["database_name"]
    output:
        "03-reports/report_{method}_{filter}.md"
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
    shell:
        "echo '# Tourmaline Report' > {output};"
        "echo '' >> {output};"
        "echo 'View this HTML report with [Chrome](https://www.google.com/chrome/){{target=\"_blank\"}} or [Firefox](https://www.mozilla.org/en-US/firefox/new/){{target=\"_blank\"}} for best results.' >> {output};"
        "echo '' >> {output};"
        "echo 'To view the linked files below: ' >> {output};"
        "echo '' >> {output};"
        "echo '* QZV (QIIME 2 visualization): click to download, then drag and drop in [https://view.qiime2.org](https://view.qiime2.org){{target=\"_blank\"}}.' >> {output};"
        "echo '* TSV (tab-separated values): click to download, then open in Microsoft Excel or Tabview (command line tool).' >> {output};"
        "echo '* PDF (portable document format): click to open and view in new tab.' >> {output};"
        "echo '* Markdown and text: click to open and view in new tab.' >> {output};"
        "echo '' >> {output};"
        "echo 'Note: Downloaded files can be deleted after viewing, as they are already stored in your Tourmaline directory.' >> {output};"
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
        "echo '### Representative Sequences Properties Summary' >> {output};"
        "echo '' >> {output};"
        "echo Markdown: \[{input.repseqsdescribe}\]\(../{input.repseqsdescribe}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        "cat {input.repseqsdescribe} >> {output};"
        "echo '' >> {output};"
        "echo '### Representative Sequences Properties Plot' >> {output};"
        "echo '' >> {output};"
        "echo PDF: \[{input.repseqspdf}\]\(../{input.repseqspdf}\){{target=\"_blank\"}} >> {output};"
        "echo '' >> {output};"
        #"echo '<img src=\"../{input.repseqspdf}\" alt=\"../{input.repseqspdf}\" width=\"70%\"/>' >> {output};"
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
        "echo 'To filter, copy 02-output-{{method}}-{{filter}}/02-alignment-tree/repseqs_to_filter_outliers.tsv to 00-data/repseqs_to_filter_{{method}}.tsv, then run Snakemake in filtered mode. A file with unassigned repseqs is provided for reference; these sequences can be filtered by adding the term \"unassigned\" to \"exclude_terms\".' >> {output};"
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
        "echo '### Taxonomy Reference Database' >> {output};"
        "echo '' >> {output};"
        "echo '{params.refdatabase}' >> {output};"
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
        "echo DEICODE Biplot QZV: \[{input.visbiplotemp}\]\(../{input.visbiplotemp}\){{target=\"_blank\"}} >> {output};"
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
    conda:
        "qiime2-2023.5"
    threads: config["other_threads"]
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
