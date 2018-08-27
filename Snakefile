import pandas as pd

configfile: "config.yaml"

# pseudo-rules ("stats" input is not included in "all")

rule deblur_se_all:
    input:
        "01-imported/fastq_pe_count_describe.tsv",
        "02-denoised/deblur-se/table.qzv",
        "02-denoised/deblur-se/table_summary_samples.txt",
        "02-denoised/deblur-se/table_summary_features.txt",
        "02-denoised/deblur-se/representative_sequences.qzv",
        "02-denoised/deblur-se/representative_sequences_lengths_describe.tsv",
        "03-repseqs/deblur-se/taxonomy.qzv",
        "03-repseqs/deblur-se/rooted_tree.qza",
        "03-repseqs/deblur-se/aligned_dna_sequences_gaps_describe.tsv",
        "04-diversity/deblur-se/alpha_rarefaction.qzv",
        "04-diversity/deblur-se/unweighted_unifrac_emperor.qzv",
        "04-diversity/deblur-se/taxa_barplot.qzv",
        "05-reports/deblur-se/report.txt"

rule deblur_se_denoise:
    input:
        "01-imported/fastq_pe_count_describe.tsv",
        "02-denoised/deblur-se/table.qzv",
        "02-denoised/deblur-se/table_summary_samples.txt",
        "02-denoised/deblur-se/table_summary_features.txt",
        "02-denoised/deblur-se/representative_sequences.qzv",
        "02-denoised/deblur-se/representative_sequences_lengths_describe.tsv"

rule deblur_se_stats:
    input:
        "04-diversity/deblur-se/unweighted_unifrac_group_significance.qzv"

rule dada2_se_all:
    input:
        "01-imported/fastq_pe_count_describe.tsv",
        "02-denoised/dada2-se/table.qzv",
        "02-denoised/dada2-se/table_summary_samples.txt",
        "02-denoised/dada2-se/table_summary_features.txt",
        "02-denoised/dada2-se/representative_sequences.qzv",
        "02-denoised/dada2-se/representative_sequences_lengths_describe.tsv",
        "03-repseqs/dada2-se/taxonomy.qzv",
        "03-repseqs/dada2-se/rooted_tree.qza",
        "03-repseqs/dada2-se/aligned_dna_sequences_gaps_describe.tsv",
        "04-diversity/dada2-se/alpha_rarefaction.qzv",
        "04-diversity/dada2-se/unweighted_unifrac_emperor.qzv",
        "04-diversity/dada2-se/taxa_barplot.qzv",
        "05-reports/dada2-se/report.txt"

rule dada2_se_denoise:
    input:
        "01-imported/fastq_pe_count_describe.tsv",
        "02-denoised/dada2-se/table.qzv",
        "02-denoised/dada2-se/table_summary_samples.txt",
        "02-denoised/dada2-se/table_summary_features.txt",
        "02-denoised/dada2-se/representative_sequences.qzv",
        "02-denoised/dada2-se/representative_sequences_lengths_describe.tsv"

rule dada2_se_stats:
    input:
        "04-diversity/dada2-se/unweighted_unifrac_group_significance.qzv"

rule dada2_pe_all:
    input:
        "01-imported/fastq_pe_count_describe.tsv",
        "02-denoised/dada2-pe/table.qzv",
        "02-denoised/dada2-pe/table_summary_samples.txt",
        "02-denoised/dada2-pe/table_summary_features.txt",
        "02-denoised/dada2-pe/representative_sequences.qzv",
        "02-denoised/dada2-pe/representative_sequences_lengths_describe.tsv",
        "03-repseqs/dada2-pe/taxonomy.qzv",
        "03-repseqs/dada2-pe/rooted_tree.qza",
        "03-repseqs/dada2-pe/aligned_dna_sequences_gaps_describe.tsv",
        "04-diversity/dada2-pe/alpha_rarefaction.qzv",
        "04-diversity/dada2-pe/unweighted_unifrac_emperor.qzv",
        "04-diversity/dada2-pe/taxa_barplot.qzv",
        "05-reports/dada2-pe/report.txt"

rule dada2_pe_denoise:
    input:
        "01-imported/fastq_pe_count_describe.tsv",
        "02-denoised/dada2-pe/table.qzv",
        "02-denoised/dada2-pe/table_summary_samples.txt",
        "02-denoised/dada2-pe/table_summary_features.txt",
        "02-denoised/dada2-pe/representative_sequences.qzv",
        "02-denoised/dada2-pe/representative_sequences_lengths_describe.tsv"

rule dada2_pe_stats:
    input:
        "04-diversity/dada2-pe/unweighted_unifrac_group_significance.qzv"

# rules

rule import_fastq_se_demux:
    input:
        config["manifest_se"]
    output:
        "01-imported/fastq_se.qza"
    shell:
        "qiime tools import "
        "--type 'SampleData[SequencesWithQuality]' "
        "--input-path {input} "
        "--output-path {output} "
        "--source-format SingleEndFastqManifestPhred33"

rule import_fastq_pe_demux:
    input:
        config["manifest_pe"]
    output:
        "01-imported/fastq_pe.qza"
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {input} "
        "--output-path {output} "
        "--source-format PairedEndFastqManifestPhred33"

# change zcat to cat if fastq files are not gzipped (but they should be)
rule count_fastq_pe_demux:
    input:
        config["manifest_pe"]
    output:
        "01-imported/fastq_pe_count.csv"
    shell:
        "for line in `tail -n +2 {input} | cut -d',' -f2`; "
        "do echo -n $line; "
        "echo -n ','; "
        "cat $line | echo $((`wc -l`/4)); "
        "done > {output}"

rule fastq_pe_count_describe:
    input:
        "01-imported/fastq_pe_count.csv"
    output:
        "01-imported/fastq_pe_count_describe.tsv"
    run:
        s = pd.read_csv(input[0], header=None)
        t = s.describe()
        t.to_csv(output[0], sep='\t')

# NOTE: refseqs.qza, reftax.qza, classifier.qza should be in external directory

rule import_ref_seqs:
    input:
        config["refseqs"]
    output:
        "01-imported/refseqs.qza"
    shell:
        "qiime tools import "
        "--type 'FeatureData[Sequence]' "
        "--input-path {input} "
        "--output-path {output}"

rule import_ref_tax:
    input:
        config["reftax"]
    output:
        "01-imported/reftax.qza"
    shell:
        "qiime tools import "
        "--type 'FeatureData[Taxonomy]' "
        '--source-format HeaderlessTSVTaxonomyFormat '
        "--input-path {input} "
        "--output-path {output}"

rule denoise_deblur_se:
    input:
        seqs="01-imported/fastq_se.qza",
        refseqs="01-imported/refseqs.qza"
    params:
        trim=config["deblur_trim"]
    output:
        "02-denoised/deblur-se"
    shell:
        "qiime deblur denoise-other "
        "--i-demultiplexed-seqs {input.seqs} "
        "--i-reference-seqs {input.refseqs} "
        "--p-trim-length {params.trim} "
        "--output-dir {output} "
        "--verbose"

rule denoise_dada2_se:
    input:
        "01-imported/fastq_se.qza"
    params:
        trunclen=config["dada_trunc_len"],
        nthreads=config["dada_n_threads"]
    output:
        "02-denoised/dada2-se"
    shell:
        "qiime dada2 denoise-single "
        "--i-demultiplexed-seqs {input} "
        "--p-trunc-len {params.trunclen} "
        "--p-n-threads {params.nthreads} "
        "--output-dir {output} "
        "--verbose"

rule denoise_dada2_pe:
    input:
        "01-imported/fastq_pe.qza"
    params:
        trunclenf=config["dada_trunc_len_f"],
        trunclenr=config["dada_trunc_len_r"],
        truncq=config["dada_trunc_q"],
        trimleftf=config["dada_trim_left_f"],
        nthreads=config["dada_n_threads"]
    output:
        "02-denoised/dada2-pe"
    shell:
        "qiime dada2 denoise-paired "
        "--i-demultiplexed-seqs {input} "
        "--p-trunc-len-f {params.trunclenf} "
        "--p-trunc-len-r {params.trunclenr} "
        "--p-trunc-q {params.truncq} "
        "--p-trim-left-f {params.trimleftf} "
        "--p-n-threads {params.nthreads} "
        "--output-dir {output} "
        "--verbose"

rule visualize_table:
    input:
        table="02-denoised/{method}/table.qza",
        metadata=config["metadata"]
    output:
        "02-denoised/{method}/table.qzv"
    shell:
        "qiime feature-table summarize "
        "--i-table {input.table} "
        "--m-sample-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule unzip_table_to_biom:
    input:
        "02-denoised/{method}/table.qza"
    output:
        "02-denoised/{method}/table.biom"
    shell:
        "unzip -o {input} -d temp; "
        "mv temp/*/data/feature-table.biom {output}; "
        "/bin/rm -r temp"

rule summarize_biom_samples:
    input:
        "02-denoised/{method}/table.biom"
    output:
        "02-denoised/{method}/table_summary_samples.txt"
    shell:
        "biom summarize-table "
        "--input-fp {input} "
        "--output-fp {output}; "
        "cat {output} | sed 's/observation/feature/g' > temp; "
        "mv temp {output}"

rule summarize_biom_features:
    input:
        "02-denoised/{method}/table.biom"
    output:
        "02-denoised/{method}/table_summary_features.txt"
    shell:
        "biom summarize-table "
        "--observations "
        "--input-fp {input} "
        "--output-fp {output}; "
        "cat {output} | sed 's/observation/feature/g' | sed 's/Counts\/sample/Counts\/feature/g' > temp; "
        "mv temp {output}"

rule visualize_repseqs:
    input:
        "02-denoised/{method}/representative_sequences.qza"
    output:
        "02-denoised/{method}/representative_sequences.qzv"
    shell:
        "qiime feature-table tabulate-seqs "
        "--i-data {input} "
        "--o-visualization {output}"

rule unzip_repseq_to_fasta:
    input:
        "02-denoised/{method}/representative_sequences.qza"
    output:
        "02-denoised/{method}/representative_sequences.fasta"
    shell:
        "unzip -o {input} -d temp; "
        "mv temp/*/data/dna-sequences.fasta {output}; "
        "/bin/rm -r temp"

rule repseq_length_distribution:
    input:
        "02-denoised/{method}/representative_sequences.fasta"
    output:
        "02-denoised/{method}/representative_sequences_lengths.txt"
    shell:
        "fastaLengthDist.pl {input} | sort > {output}"

rule repseq_length_distribution_describe:
    input:
        "02-denoised/{method}/representative_sequences_lengths.txt"
    output:
        "02-denoised/{method}/representative_sequences_lengths_describe.tsv"
    run:
        s = pd.read_csv(input[0], header=None)
        t = s.describe()
        t.to_csv(output[0], sep='\t')

# skipping this step for now because it is slow with degenerate primers
# rule feature_classifier_extract_reads 
#     input:
#         "01-imported/refseqs.qza"
#     output:
#         "01-imported/refseqs_extracted.qza"
#     shell:
#         "qiime feature-classifier extract-reads "
#         "--i-sequences {input} "
#         "--p-f-primer {params.fprimer} "
#         "--p-r-primer {params.rprimer} "
#         "--o-reads {output}"

# this step should be done only once for each amplicon, outside of snakemake
# then specify classifer (eg "classifier_18s_fprimer_rprimer.qza") in configfile
rule feature_classifier_fit_classifier_naive_bayes:
    input:
        seqs="01-imported/refseqs.qza",
        tax="01-imported/reftax.qza"
    output:
        "01-imported/classifier.qza"
    shell:
        "qiime feature-classifier fit-classifier-naive-bayes "
        "--i-reference-reads {input.seqs} "
        "--i-reference-taxonomy {input.tax} "
        "--o-classifier {output}"

rule feature_classifier_classify_sklearn:
    input:
        classifier="01-imported/classifier.qza",
        repseqs="02-denoised/{method}/representative_sequences.qza"
    output:
        "03-repseqs/{method}/taxonomy.qza"
    shell:
        "qiime feature-classifier classify-sklearn "
        "--i-classifier {input.classifier} "
        "--i-reads {input.repseqs} "
        "--o-classification {output}"

rule visualize_taxonomy:
    input:
        "03-repseqs/{method}/taxonomy.qza"
    output:
        "03-repseqs/{method}/taxonomy.qzv"
    shell:
        "qiime metadata tabulate "
        "--m-input-file {input} "
        "--o-visualization {output}"

rule alignment_mafft:
    input:
        "02-denoised/{method}/representative_sequences.qza"
    output:
        "03-repseqs/{method}/aligned_representative_sequences.qza"
    shell:
        "qiime alignment mafft "
        "--i-sequences {input} "
        "--o-alignment {output}"

rule alignment_mask:
    input:
        "03-repseqs/{method}/aligned_representative_sequences.qza"
    output:
        "03-repseqs/{method}/masked_aligned_representative_sequences.qza"
    shell:
        "qiime alignment mask "
        "--i-alignment {input} "
        "--o-masked-alignment {output}"

rule phylogeny_fasttree:
    input:
        "03-repseqs/{method}/masked_aligned_representative_sequences.qza"
    output:
        "03-repseqs/{method}/unrooted_tree.qza"
    shell:
        "qiime phylogeny fasttree "
        "--i-alignment {input} "
        "--o-tree {output}"

rule phylogeny_midpoint_root:
    input:
        "03-repseqs/{method}/unrooted_tree.qza"
    output:
        "03-repseqs/{method}/rooted_tree.qza"
    shell:
        "qiime phylogeny midpoint-root "
        "--i-tree {input} "
        "--o-rooted-tree {output}"

# skipping this step because tree visualization in q2 is still under construction
# rule visualize_tree:
#     input:
#         "03-repseqs/{method}/rooted_tree.qza"
#     output:
#         "03-repseqs/{method}/rooted_tree.qzv"
#     shell:

rule unzip_alignment_to_fasta:
    input:
        "03-repseqs/{method}/masked_aligned_representative_sequences.qza"
    output:
        "03-repseqs/{method}/aligned_dna_sequences.fasta"
    shell:
        "unzip -o {input} -d temp; "
        "mv temp/*/data/aligned-dna-sequences.fasta {output}; "
        "/bin/rm -r temp"

rule alignment_count_gaps:
    input:
        "03-repseqs/{method}/aligned_dna_sequences.fasta"
    output:
        "03-repseqs/{method}/aligned_dna_sequences_gaps.txt"
    shell:
        "cat {input} | grep -v '>' | sed 's/[^-]//g' | awk '{{ print length }}' > {output}"

rule alignment_gaps_distribution_describe:
    input:
        "03-repseqs/{method}/aligned_dna_sequences_gaps.txt"
    output:
        "03-repseqs/{method}/aligned_dna_sequences_gaps_describe.tsv"
    run:
        s = pd.read_csv(input[0], header=None)
        t = s.describe()
        t.to_csv(output[0], sep='\t')

rule diversity_alpha_rarefaction:
    input:
        table="02-denoised/{method}/table.qza",
        phylogeny="03-repseqs/{method}/rooted_tree.qza",
        metadata=config["metadata"]
    params:
        maxdepth=config["alpha_max_depth"]
    output:
        "04-diversity/{method}/alpha_rarefaction.qzv"
    shell:
        "qiime diversity alpha-rarefaction "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-max-depth {params.maxdepth} "
        "--p-metrics faith_pd "
        "--p-metrics shannon "
        "--p-metrics observed_otus "
        "--p-metrics chao1 "
        "--m-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule diversity_core_metrics_phylogenetic:
    input:
        table="02-denoised/{method}/table.qza",
        phylogeny="03-repseqs/{method}/rooted_tree.qza",
        metadata=config["metadata"]
    params:
        samplingdepth=config["core_sampling_depth"]
    output:
        rarefiedtable="04-diversity/{method}/rarefied_table.qza",
        faithpdvector="04-diversity/{method}/faith_pd_vector.qza",        
        observedotusvector="04-diversity/{method}/observed_otus_vector.qza",
        shannonvector="04-diversity/{method}/shannon_vector.qza",
        evennessvector="04-diversity/{method}/evenness_vector.qza",
        unweightedunifracdistancematrix="04-diversity/{method}/unweighted_unifrac_distance_matrix.qza",
        weightedunifracdistancematrix="04-diversity/{method}/weighted_unifrac_distance_matrix.qza",
        jaccarddistancematrix="04-diversity/{method}/jaccard_distance_matrix.qza",
        braycurtisdistancematrix="04-diversity/{method}/bray_curtis_distance_matrix.qza",
        unweightedunifracpcoaresults="04-diversity/{method}/unweighted_unifrac_pcoa_results.qza",
        weightedunifracpcoaresults="04-diversity/{method}/weighted_unifrac_pcoa_results.qza",
        jaccardpcoaresults="04-diversity/{method}/jaccard_pcoa_results.qza",
        braycurtispcoaresults="04-diversity/{method}/bray_curtis_pcoa_results.qza",
        unweightedunifracemperor="04-diversity/{method}/unweighted_unifrac_emperor.qzv",
        weightedunifracemperor="04-diversity/{method}/weighted_unifrac_emperor.qzv",
        jaccardemperor="04-diversity/{method}/jaccard_emperor.qzv",
        braycurtisemperor="04-diversity/{method}/bray_curtis_emperor.qzv"
    shell:
        "qiime diversity core-metrics-phylogenetic "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-sampling-depth {params.samplingdepth} "
        "--m-metadata-file {input.metadata} "
        "--o-rarefied-table {output.rarefiedtable} "
        "--o-faith-pd-vector {output.faithpdvector} "
        "--o-observed-otus-vector {output.observedotusvector} "
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
        "--o-bray-curtis-emperor {output.braycurtisemperor}"

rule diversity_beta_group_significance:
    input:
        distancematrix="04-diversity/{method}/unweighted_unifrac_distance_matrix.qza",
        metadata=config["metadata"]
    params:
        column=config["beta_group_column"]
    output:
        "04-diversity/{method}/unweighted_unifrac_group_significance.qzv"
    shell:
        "qiime diversity beta-group-significance "
        "--i-distance-matrix {input.distancematrix} "
        "--m-metadata-file {input.metadata} "
        "--m-metadata-column {params.column} "
        "--o-visualization {output} "
        "--p-pairwise"

rule taxa_barplot:
    input:
        table="02-denoised/{method}/table.qza",
        taxonomy="03-repseqs/{method}/taxonomy.qza",
        metadata=config["metadata"]
    output:
        "04-diversity/{method}/taxa_barplot.qzv"
    shell:
        "qiime taxa barplot "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--m-metadata-file {input.metadata} "
        "--o-visualization {output}"

rule report:
    input:
        fastq="01-imported/fastq_pe_count_describe.tsv",
        samples="02-denoised/{method}/table_summary_samples.txt",
        features="02-denoised/{method}/table_summary_features.txt",
        repseqs="02-denoised/{method}/representative_sequences_lengths_describe.tsv",
        alignment="03-repseqs/{method}/aligned_dna_sequences_gaps_describe.tsv"
    output:
        "05-reports/{method}/report.txt"
    shell:
        "head -n 13 {input.fastq} {input.samples} {input.features} {input.repseqs} {input.alignment} > {output}"

