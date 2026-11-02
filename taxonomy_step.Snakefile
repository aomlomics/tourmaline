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
    if file is None or file == "":
        return False
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

if config["classify_method"] == "revamp":
    # REVAMP classifies against a local NCBI nt BLAST database, not a QIIME classifier
    # or reference artifacts; refseqs_file and taxa_file are unused.
    use_classifier="no"
elif config["pretrained_classifier"] != None:
    use_classifier="yes"
else:
    use_classifier="no"
    print(f"No pretrained classifier provided, using refseqs and reftax files.\n")
    if config["refseqs_file"] == None or config["taxa_file"] == None:
        print(f"ERROR: refseqs_file and taxa_file must be provided if pretrained_classifier is not used.\n")

# consensus-blast: blastn only honours --p-num-threads against a pre-indexed
# database. Given --i-reference-reads, q2-feature-classifier runs
# `blastn -subject`, which is single-threaded whatever classify_threads says.
# Three modes, in precedence order:
#   blast_database set           import that database       (threaded)
#   build_blast_database: true   makeblastdb from refseqs   (threaded)
#   neither                      -subject, exactly as before (the default)
# Read with .get so configs written before these options still parse.
BLASTDB_SUFFIXES = (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto", ".njs")
blastdb_qza = taxonomy_dir + "blastdb.qza"


def validate_blast_database(path):
    """Fail at parse time if *path* is not a QIIME-importable BLAST database.

    QIIME 2's BLASTDBDirFmtV5 requires all of BLASTDB_SUFFIXES under a single
    basename, with none optional, so the common cases that cannot work are
    worth naming here rather than surfacing as an import traceback much later.
    """
    if not os.path.isdir(path):
        raise ValueError(
            f"blast_database is not a directory: {path}\n"
            "Point it at the directory holding the database files (the ones "
            "makeblastdb wrote), not at a file or a database basename."
        )
    names = os.listdir(path)
    volumes = sorted(n for n in names if n.endswith((".nal", ".pal")))
    if volumes:
        raise ValueError(
            f"blast_database is a multi-volume BLAST database ({volumes[0]}): {path}\n"
            "QIIME 2's BLASTDB type cannot represent these. Either rebuild it as "
            "a single volume (makeblastdb without -max_file_sz splitting), or "
            "leave blast_database unset and set build_blast_database: true to build one "
            "from refseqs_file."
        )
    basenames = {
        n[: -len(suffix)]
        for n in names
        for suffix in BLASTDB_SUFFIXES
        if n.endswith(suffix)
    }
    if not basenames:
        raise ValueError(
            f"blast_database contains no BLAST database files: {path}\n"
            f"Expected one set of {', '.join(BLASTDB_SUFFIXES)} files."
        )
    complete = [
        b for b in sorted(basenames)
        if all(b + suffix in names for suffix in BLASTDB_SUFFIXES)
    ]
    if len(complete) > 1:
        raise ValueError(
            f"blast_database holds more than one BLAST database ({', '.join(complete)}): "
            f"{path}\nQIIME 2 imports a directory, not a basename, so give each "
            "database its own directory."
        )
    if not complete:
        base = sorted(basenames)[0]
        missing = [sfx for sfx in BLASTDB_SUFFIXES if base + sfx not in names]
        hint = ""
        if missing == [".njs"]:
            hint = (
                "\nOnly .njs is missing, which BLAST writes from 2.13.0 onward: "
                "this database predates it. Rebuild with a newer makeblastdb, or "
                "set build_blast_database: true."
            )
        elif ".ndb" in missing:
            hint = (
                "\nA missing .ndb usually means a version-4 database; rebuild "
                "with `makeblastdb -blastdb_version 5`."
            )
        raise ValueError(
            f"blast_database is missing required files for database {base!r}: "
            f"{', '.join(missing)}\n{path}\nQIIME 2's BLASTDB type requires all "
            f"of {', '.join(BLASTDB_SUFFIXES)}.{hint}"
        )
    return complete[0]


blast_database = config.get("blast_database") or None
build_blast_database = bool(config.get("build_blast_database", False))
blastdb_source = None
if config["classify_method"] == "consensus-blast":
    if blast_database:
        validate_blast_database(blast_database)
        blastdb_source = "import"
        if build_blast_database:
            print(
                "Both blast_database and build_blast_database are set; using blast_database "
                "and not building one.\n"
            )
    elif build_blast_database:
        blastdb_source = "build"
elif blast_database or build_blast_database:
    print(
        "blast_database / build_blast_database apply to consensus-blast only; ignoring "
        f"them for classify_method {config['classify_method']}.\n"
    )

# Krona plots are opt-in: they need the `krona` conda env, which most runs don't have.
# Read with .get so configs written before this option still parse.
make_krona = config.get("make_krona", False)
krona_dir = taxonomy_dir + "figures/krona_inputs/"
krona_manifest = krona_dir + "krona_datasets.tsv"
krona_html = taxonomy_dir + "figures/" + config["run_name"] + "-krona.html"

rule run_taxonomy:
    input:
        taxonomy_tsv,
        taxonomy_dir + "figures/" + config["run_name"] + "-taxa_barplot.qzv",
        taxonomy_dir + config["run_name"] + "-taxa_sample_table_" + "l" + str(config["collapse_taxalevel"]) + ".tsv",
        taxonomy_dir + config["run_name"] + "-asv_taxa_features.tsv",
        *([krona_html] if make_krona else [])

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

if config["classify_method"] == "revamp":
    output_seq = None
    output_tax = None
elif has_fa_suffix(config["refseqs_file"], fasta_suffixes) and config["classify_method"] != "bt2-blca":
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

if make_krona:
    rule krona_inputs:
        input:
            table=input_table,
            taxonomy=taxonomy_qza
        output:
            manifest=krona_manifest
        params:
            outdir=krona_dir,
            persample="yes" if config.get("krona_per_sample", True) else "no"
        conda:
            "qiime2-amplicon-2024.10"
        shell:
            "python scripts/taxonomy_to_krona.py "
            "--table {input.table} "
            "--taxonomy {input.taxonomy} "
            "--outdir {params.outdir} "
            "--manifest {output.manifest} "
            "--per-sample {params.persample}"

    rule krona_plot:
        input:
            manifest=krona_manifest
        output:
            krona_html
        conda:
            "krona"
        shell:
            # The manifest gives each dataset as `path<TAB>label`; Krona wants `path,label`.
            "ktImportText -o {output} "
            "$(awk -F'\\t' '{{print $1\",\"$2}}' {input.manifest} | tr '\\n' ' ')"

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
