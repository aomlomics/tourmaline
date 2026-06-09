## Tourmaline taxonomy Snakemake workflow.
## Invoked via `tourmaline.sh --step taxonomy`; classifies representative
## sequences, exports per-feature taxonomy tables, and produces visual summaries.
import os
import shutil

output_dir = config["output_dir"]+"/"
taxonomy_dir = output_dir + config["run_name"] + "-taxonomy/"
taxonomy_qza = taxonomy_dir + config["run_name"] + "-taxonomy.qza"
taxonomy_tsv = taxonomy_dir + config["run_name"] + "-taxonomy.tsv"
classifier_qza = taxonomy_dir + "classifier.qza"

# Copy config file to output directory
config_output_path = taxonomy_dir + config["run_name"] + "-taxonomy_config.yaml"
os.makedirs(os.path.dirname(config_output_path), exist_ok=True)
shutil.copy(workflow.configfiles[0], config_output_path)

def has_fa_suffix(file, suffixes):
    return any(file.endswith(suffix) for suffix in suffixes)

def change_suffix(file, new_suffix):
    base_name = os.path.basename(file)
    file_name = os.path.splitext(base_name)[0]
    return file_name + new_suffix

if config["sample_metadata_file"] != None:
    use_metadata="yes"
else:
    use_metadata="no"

input_table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza"
if config["repseqs_run_name"] != None:
    repseqs_run_name=config["repseqs_run_name"]
    input_repseqs=output_dir+repseqs_run_name+"-repseqs/"+repseqs_run_name+"-repseqs.qza"
    input_table=output_dir+repseqs_run_name+"-repseqs/"+repseqs_run_name+"-table.qza"
elif config["repseqs_qza_file"] != None:
    input_repseqs=config["repseqs_qza_file"]
    input_table=config["table_qza_file"]
else:
    input_repseqs=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-repseqs.qza"
    input_table=output_dir+config["run_name"]+"-repseqs/"+config["run_name"]+"-table.qza"

fasta_suffixes = [".fna", ".fa",".fasta"]

if config["pretrained_classifier"] != None:
    use_classifier="yes"
else:
    use_classifier="no"
    print(f"No pretrained classifier provided, using refseqs and reftax files.\n")
    if config["refseqs_file"] == None or config["taxa_file"] == None:
        print(f"ERROR: refseqs_file and taxa_file must be provided if pretrained_classifier is not used.\n")

rule run_taxonomy:
    input:
        taxonomy_tsv,
        taxonomy_dir + "figures/" + config["run_name"] + "-taxa_barplot.qzv",
        taxonomy_dir + config["run_name"] + "-taxa_sample_table_" + "l" + str(config["collapse_taxalevel"]) + ".tsv",
        taxonomy_dir + config["run_name"] + "-asv_taxa_features.tsv"

if config["classify_method"] == "bt2-blca":
    if has_fa_suffix(input_repseqs, [".qza"]):
        fasta_repseqs=taxonomy_dir+change_suffix(input_repseqs, ".fasta")
        rule export_repseqs_to_fasta:
            input:
                input_repseqs
            output:
                fasta_repseqs
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                "qiime tools export "
                "--input-path {input} "
                "--output-path {output} "
                "--output-format DNAFASTAFormat; "
    else:
        fasta_repseqs = input_repseqs

if has_fa_suffix(config["refseqs_file"], fasta_suffixes) and config["classify_method"] != "bt2-blca":
    output_seq = taxonomy_dir+change_suffix(config["refseqs_file"], ".qza")
    output_tax = taxonomy_dir+change_suffix(config["taxa_file"], ".qza")
    rule import_ref_seqs:
        input:
            config["refseqs_file"]
        output:
            output_seq
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "qiime tools import "
            "--type 'FeatureData[Sequence]' "
            "--input-path {input} "
            "--output-path {output}"

    rule import_ref_tax:
        input:
            config["taxa_file"]
        output:
            output_tax
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "qiime tools import "
            "--type 'FeatureData[Taxonomy]' "
            "--input-format HeaderlessTSVTaxonomyFormat "
            "--input-path {input} "
            "--output-path {output}"
elif has_fa_suffix(config["refseqs_file"], [".qza"]) and config["classify_method"] != "bt2-blca":
    output_seq = config["refseqs_file"]
    output_tax = config["taxa_file"]
elif config["classify_method"] == "bt2-blca":
    if has_fa_suffix(config["refseqs_file"], [".qza"]):
        output_seq=taxonomy_dir+change_suffix(config["refseqs_file"], ".fasta")
        output_tax = taxonomy_dir+change_suffix(config["taxa_file"], ".txt")
        rule export_refseqs_to_fasta:
            input:
                config["refseqs_file"]
            output:
                output_seq
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                "qiime tools export "
                "--input-path {input} "
                "--output-path {output} "
                "--output-format DNAFASTAFormat"
        rule export_ref_taxonomy_to_tsv:
            input:
                config["taxa_file"]
            output:
                output_tax
            conda:
                "qiime2-amplicon-2024.10"
            shell:
                "qiime tools export "
                "--input-path {input} "
                "--output-path {output} "
                "--output-format TSVTaxonomyFormat"
    else:
        output_seq = config["refseqs_file"]
        output_tax = config["taxa_file"]
elif use_classifier == "no":
    raise ValueError("refseqs_file must have one of the following extensions: .qza, .fna, .fa, .fasta")
else:
    if config["pretrained_classifier"] != None:
        print(f"Using pretrained classifier.\n")
    elif config["refseqs_file"] == None or config["taxa_file"] == None:
        print(f"ERROR: refseqs_file and taxa_file must be provided if pretrained_classifier is not used.\n")

include: "rules/taxonomy_assignment.smk"

rule export_taxa_biom:
    input:
        table=input_table,
        taxonomy=taxonomy_qza,
    output:
        taxa_table=taxonomy_dir + config["run_name"] + "-taxa_sample_table_" + "l" + str(config["collapse_taxalevel"]) + ".tsv",
    params:
        taxalevel=config["collapse_taxalevel"]
    conda:
        "qiime2-amplicon-2024.10"
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
        "-o TEMP.tsv "
        "--to-tsv "
        "&& cat TEMP.tsv | tail -n +2 | sed 's/^#OTU ID/taxonomy/' > {output.taxa_table} "
        "&& /bin/rm -r tempfile_collapsed.qza temp_export/ TEMP.tsv"

rule export_asv_taxa_features:
    input:
        taxonomy=taxonomy_qza,
        repseqs=input_repseqs
    output:
        taxonomy_dir + config["run_name"] + "-asv_taxa_features.tsv"
    params:
        taxaranks=config["taxa_ranks"]
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        "python scripts/create_asv_seq_taxa_output.py --input_repseqs {input.repseqs} --input_taxonomy {input.taxonomy} --output {output} --taxaranks {params.taxaranks}"

rule taxa_barplot:
    input:
        table=input_table,
        taxonomy=taxonomy_qza,
    output:
        taxonomy_dir + "figures/" + config["run_name"] + "-taxa_barplot.qzv"
    params:
        metadata=config["sample_metadata_file"]
    conda:
        "qiime2-amplicon-2024.10"
    shell:
        """
        if [ {use_metadata} == "yes" ]; then
            qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --m-metadata-file {params.metadata} \
            --o-visualization {output};
        else
            qiime taxa barplot \
            --i-table {input.table} \
            --i-taxonomy {input.taxonomy} \
            --o-visualization {output};
        fi;
        """
