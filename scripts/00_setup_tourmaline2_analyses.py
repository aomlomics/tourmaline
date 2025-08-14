#!/usr/bin/env python3

import os
import argparse
import shutil
import glob
import yaml
import pandas as pd

def create_output_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Creating output folder: {path}")
    else:
        print(f"Output folder already exists: {path}")
    return path

def copy_configs(tourmaline_folder, output_folder, config_to_copy=None, run_name=None, classifier_method=None, database_name=None):
    # Look for original config files without 00_ prefix
    if config_to_copy:
        # Remove 00_ prefix and convert underscores to hyphens to match original files
        original_name = config_to_copy.replace('00_', '').replace('_', '-')
        configs = [os.path.join(tourmaline_folder, original_name)]
    else:
        configs = glob.glob(os.path.join(tourmaline_folder, "config-*.yaml"))
    
    for config_file in configs:
        if not os.path.exists(config_file):
            print(f"Warning: Config file {config_file} not found. Skipping.")
            continue
            
        config_filename = os.path.basename(config_file)
        # Add 00_ prefix and replace hyphens with underscores
        new_config_filename = f"00_{config_filename.replace('-', '_')}"
        
        # Add run name and other identifiers to config filenames
        if "01_qaqc" in new_config_filename:
            new_config_filename = f"00_config_01_qaqc_{run_name}.yaml"
        elif "02_repseqs" in new_config_filename:
            new_config_filename = f"00_config_02_repseqs_{run_name}.yaml"
        elif "03_taxonomy" in new_config_filename:
            new_config_filename = f"00_config_03_taxonomy_{classifier_method}_{database_name}_{run_name}.yaml"
        
        destination_path = os.path.join(output_folder, new_config_filename)
        if not os.path.exists(destination_path):
            # Copy the original file and rename it with 00_ prefix and underscores
            shutil.copy(config_file, destination_path)
            print(f"Copying config: {config_filename} -> {new_config_filename}")
        else:
            print(f"Config file {new_config_filename} already exists in {output_folder}. Skipping copy.")

def generate_manifest(input_data_path, output_folder):
    fastq_files = sorted(glob.glob(os.path.join(input_data_path, "*.fastq.gz")))
    manifest_data = []

    for fq in fastq_files:
        file_name = os.path.basename(fq)
        direction = "forward" if "_R1" in file_name else "reverse" if "_R2" in file_name else None
        if direction:
            sample_id = file_name.split("_R")[0]
            manifest_data.append({
                "sample-id": sample_id,
                "absolute-filepath": os.path.abspath(fq),
                "direction": direction
            })

    manifest_df = pd.DataFrame(manifest_data)
    manifest_path = os.path.join(output_folder, "00_manifest_pe.csv")
    manifest_df.to_csv(manifest_path, index=False)
    print(f"Creating manifest: {manifest_path}")
    return os.path.abspath(manifest_path), manifest_df

def create_metadata_from_manifest(manifest_df, metadata_path):
    sample_ids = manifest_df['sample-id'].drop_duplicates()
    metadata_df = pd.DataFrame({'sample_name': sample_ids})
    metadata_df.to_csv(metadata_path, sep="\t", index=False)
    print(f"Creating metadata file: {metadata_path}")

def update_config_01_qaqc(config_path, manifest_path, primers, run_name, output_dir, threads):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    config["sample_manifest_file"] = os.path.abspath(manifest_path)
    config["run_name"] = f"01_{run_name}"
    config["output_dir"] = output_dir
    config["trimming_threads"] = threads

    if primers:
        fwd, rev = primers.split("/")
        config["fwd_primer"] = fwd
        config["rev_primer"] = rev

    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    print(f"Updating config: {config_path}")

def update_config_02_repseqs(config_path, run_name, output_dir, metadata_file, threads):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    config["run_name"] = f"02_{run_name}"
    config["output_dir"] = output_dir
    config["sample_metadata_file"] = os.path.abspath(metadata_file)
    config["sample_run_name"] = None
    config["asv_threads"] = threads
    
    # Update fastq_qza_file path using run name from config-01
    qaqc_run_name = f"01_{run_name}"
    config["fastq_qza_file"] = os.path.join(output_dir, f"{qaqc_run_name}-qaqc", f"{qaqc_run_name}_fastq.qza")

    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    print(f"Updating config: {config_path}")

def update_config_03_taxonomy(config_path, run_name, output_dir, metadata_file, classifier_method, database_path, output_subfolder, threads):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Get database name without path for naming purposes
    database_name = os.path.basename(database_path)
    if database_name.endswith('.qza'):
        database_name = database_name[:-4]

    # Update run name with classifier method and database, using output_subfolder
    config["run_name"] = f"03_{output_subfolder}_{classifier_method}_{database_name}_{run_name}"
    config["output_dir"] = output_dir
    config["sample_metadata_file"] = os.path.abspath(metadata_file)
    config["classify_threads"] = threads
    
    # Update repseqs and table file paths
    repseqs_run_name = f"02_{run_name}"
    config["repseqs_run_name"] = None
    config["repseqs_qza_file"] = os.path.join(output_dir, f"{repseqs_run_name}-repseqs", f"{repseqs_run_name}-repseqs.qza")
    config["table_qza_file"] = os.path.join(output_dir, f"{repseqs_run_name}-repseqs", f"{repseqs_run_name}-table.qza")
    
    # Update classifier method and database name
    config["classify_method"] = classifier_method
    config["database_name"] = database_name
    
    # Update database paths based on classifier method using absolute paths
    if classifier_method in ['consensus-blast', 'consensus-vsearch']:
        # For consensus methods, use the provided database path directly
        config["refseqs_file"] = os.path.abspath(f"{database_path}-seqs.qza")
        config["taxa_file"] = os.path.abspath(f"{database_path}-tax.qza")
    elif classifier_method == 'naive-bayes':
        # For naive-bayes, use the provided classifier directly and set refseqs and taxa to null
        config["pretrained_classifier"] = os.path.abspath(database_path)
        config["refseqs_file"] = None
        config["taxa_file"] = None
    elif classifier_method == 'bt2-blca':
        # For bt2-blca, use the provided database path directly
        config["bowtie_database"] = os.path.abspath(database_path)

    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    print(f"Updating config: {config_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--read-folder", required=True)
    parser.add_argument("--working-directory", required=True)
    parser.add_argument("--data-type", required=True)
    parser.add_argument("--tourmaline2-folder", required=True)
    parser.add_argument("--run-name", default="")
    parser.add_argument("--primers", default="")
    parser.add_argument("--metadata-folder", help="Optional: Path to existing metadata files. If not provided, metadata will be created from manifest.")
    parser.add_argument(
        "--config-file",
        default="all",
        choices=["00_config_01_qaqc.yaml", "00_config_02_repseqs.yaml", "00_config_03_taxonomy.yaml", "all"],
        help="Specify which config file to update. Default is 'all'."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force regeneration of manifest and metadata files if they already exist."
    )
    parser.add_argument(
        "--classifier-method",
        choices=['naive-bayes', 'consensus-blast', 'consensus-vsearch', 'bt2-blca'],
        required=True,
        help="Taxonomy classification method"
    )
    parser.add_argument(
        "--database",
        required=True,
        help="Path to the database for taxonomy classification"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=6,
        help="Number of threads to use for processing (default: 6)"
    )
    args = parser.parse_args()

    # Create metadata folder if provided
    if args.metadata_folder:
        os.makedirs(args.metadata_folder, exist_ok=True)

    for subdir in sorted(os.listdir(args.read_folder)):
        input_data_path = os.path.join(args.read_folder, subdir, args.data_type)
        if not os.path.exists(input_data_path):
            continue

        output_subfolder = f"{subdir}_{args.data_type}"
        full_output_path = os.path.join(args.working_directory, args.data_type, output_subfolder)
        create_output_folder(full_output_path)
        
        # Get database name without path
        database_name = os.path.basename(args.database)
        if database_name.endswith('.qza'):
            database_name = database_name[:-4]
        
        # Copy configs based on selection
        copy_configs(args.tourmaline2_folder, full_output_path, args.config_file if args.config_file != "all" else None,
                    args.run_name, args.classifier_method, database_name)

        # Manifest generation
        manifest_filename = "00_manifest_pe.csv"
        manifest_path_abs = os.path.join(full_output_path, manifest_filename)
        
        if os.path.exists(manifest_path_abs) and not args.force:
            print(f"Manifest file {manifest_path_abs} already exists. Using existing file.")
            manifest_df = pd.read_csv(manifest_path_abs)
            manifest_path_for_config = manifest_path_abs
        else:
            if args.force and os.path.exists(manifest_path_abs):
                print(f"Force regenerating manifest file: {manifest_path_abs}")
            else:
                print(f"Generating manifest file: {manifest_path_abs}")
            manifest_path_for_config, manifest_df = generate_manifest(input_data_path, full_output_path)

        # Metadata handling
        if args.metadata_folder:
            # Use existing metadata from provided folder
            metadata_filename = f"00_{output_subfolder}_metadata.tsv"
            metadata_path_abs = os.path.join(args.metadata_folder, metadata_filename)
            if not os.path.exists(metadata_path_abs):
                print(f"Warning: Metadata file not found at {metadata_path_abs}")
                print("Creating metadata from manifest in the working directory...")
                metadata_path_abs = os.path.join(full_output_path, "00_metadata.tsv")
                create_metadata_from_manifest(manifest_df, metadata_path_abs)
        else:
            # Create metadata in the same folder as manifest
            metadata_path_abs = os.path.join(full_output_path, "00_metadata.tsv")
            if os.path.exists(metadata_path_abs) and not args.force:
                print(f"Metadata file {metadata_path_abs} already exists. Using existing file.")
            else:
                if args.force and os.path.exists(metadata_path_abs):
                    print(f"Force regenerating metadata file: {metadata_path_abs}")
                else:
                    print(f"Generating metadata file: {metadata_path_abs}")
                create_metadata_from_manifest(manifest_df, metadata_path_abs)

        run_name_full = f"{output_subfolder}_{args.run_name}" if args.run_name else output_subfolder
        full_output_dir_abs = os.path.abspath(os.path.join(full_output_path, f"{output_subfolder}_Output"))

        config_files_to_update = []
        if args.config_file == "all":
            config_files_to_update = ["00_config_01_qaqc.yaml", "00_config_02_repseqs.yaml", "00_config_03_taxonomy.yaml"]
        else:
            config_files_to_update = [args.config_file]

        if "00_config_01_qaqc.yaml" in config_files_to_update:
            config_qaqc_path = os.path.join(full_output_path, f"00_config_01_qaqc_{args.run_name}.yaml")
            if os.path.exists(config_qaqc_path):
                update_config_01_qaqc(config_qaqc_path, manifest_path_for_config, args.primers, run_name_full, full_output_dir_abs, args.threads)
            else:
                print(f"Warning: Config file {config_qaqc_path} not found. Skipping update.")
        
        if "00_config_02_repseqs.yaml" in config_files_to_update:
            config_repseqs_path = os.path.join(full_output_path, f"00_config_02_repseqs_{args.run_name}.yaml")
            if os.path.exists(config_repseqs_path):
                update_config_02_repseqs(config_repseqs_path, run_name_full, full_output_dir_abs, metadata_path_abs, args.threads)
            else:
                print(f"Warning: Config file {config_repseqs_path} not found. Skipping update.")

        if "00_config_03_taxonomy.yaml" in config_files_to_update:
            config_taxonomy_path = os.path.join(full_output_path, f"00_config_03_taxonomy_{args.classifier_method}_{database_name}_{args.run_name}.yaml")
            if os.path.exists(config_taxonomy_path):
                update_config_03_taxonomy(config_taxonomy_path, run_name_full, full_output_dir_abs, metadata_path_abs, 
                                       args.classifier_method, args.database, output_subfolder, args.threads)
            else:
                print(f"Warning: Config file {config_taxonomy_path} not found. Skipping update.")

if __name__ == "__main__":
    main()


