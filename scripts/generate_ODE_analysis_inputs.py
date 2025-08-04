import os
import glob
import subprocess
import argparse
import shutil

parser = argparse.ArgumentParser(description="Generate analysis metadata commands.")
parser.add_argument('--working-directory', '--wd', required=True, help='Base directory to search for analyses.')
parser.add_argument('--classify-method', required=True, help='Classification method: consensus-blast, naive-bayes, consensus-vsearch, bt2-blca.')
parser.add_argument('--database', required=True, help='Database string for taxonomy config.')
parser.add_argument('--project-id', required=True, help='Project ID string.')
parser.add_argument('--assay-name', required=True, help='Assay name string.')
parser.add_argument('--output-folder-path', '-O', default=None, help='If provided, output files will be written here.')
parser.add_argument('--ode-inputs', default='FALSE', choices=['TRUE', 'FALSE'], help='If TRUE and -O is provided, copy ODE input files.')
parser.add_argument('--tourmaline-path', required=True, help='Path to the Tourmaline GitHub repository.')
args = parser.parse_args()

BASE_DIR = args.working_directory
CLASSIFY_METHOD_HYPHEN = args.classify_method
CLASSIFY_METHOD_UNDERSCORE = args.classify_method.replace('-', '_')
DATABASE = args.database
PROJECT_ID = args.project_id
ASSAY_NAME = args.assay_name
OUTPUT_FOLDER_PATH = args.output_folder_path
ODE_INPUTS = args.ode_inputs.upper() == 'TRUE'

TOURMALINE_GITHUB = os.path.join(TOURMALINE_PATH, "scripts/format_analysisMetadata.py")
TOURMALINE_METADATA = os.path.join(TOURMALINE_PATH, "00-data/tourmaline_metadata.yaml")

for folder in sorted(os.listdir(BASE_DIR)):
    folder_path = os.path.join(BASE_DIR, folder)
    if not os.path.isdir(folder_path):
        continue

    # -s: latest 00_config_01_qaqc*.yaml
    qaqc_files = sorted(glob.glob(os.path.join(folder_path, "00_config_01_qaqc*.yaml")))
    if not qaqc_files:
        print(f"Warning: No qaqc config found in {folder}")
        print()
        continue
    s_file = qaqc_files[-1]

    # -r: latest 00_config_02_repseqs*.yaml
    repseqs_files = sorted(
        glob.glob(os.path.join(folder_path, "00_config_02_repseqs*.yaml"))
    )
    if not repseqs_files:
        print(f"Warning: No repseqs config found in {folder}")
        print()
        continue
    r_file = repseqs_files[-1]

    # -t: latest taxonomy config with classify_method and database (underscores)
    taxonomy_pattern = f"00_config_03_taxonomy_{CLASSIFY_METHOD_UNDERSCORE}_{DATABASE}*.yaml"
    taxonomy_files = sorted(
        glob.glob(os.path.join(folder_path, taxonomy_pattern))
    )
    if not taxonomy_files:
        print(f"Warning: No taxonomy config found in {folder} for pattern {taxonomy_pattern}")
        print()
        continue
    t_file = taxonomy_files[-1]

    # -p, -a
    project_id = PROJECT_ID
    assay_name = ASSAY_NAME

    # -A: FOLDERNAME_LATESTTAXONOMY (strip 00_config_03_taxonomy_ from filename)
    taxonomy_base = os.path.basename(t_file)
    taxonomy_suffix = taxonomy_base.replace("00_config_03_taxonomy_", "")
    analysis_run_name = f"{folder}_{taxonomy_suffix.replace('.yaml','')}"

    # -o: analysisMetadata-ANALYSIS_RUN_NAME.tsv
    if OUTPUT_FOLDER_PATH:
        output_file = os.path.join(OUTPUT_FOLDER_PATH, f"analysisMetadata-{analysis_run_name}.tsv")
    else:
        output_file = os.path.join(folder_path, f"analysisMetadata-{analysis_run_name}.tsv")

    # Build the command as a list for subprocess
    cmd = [
        "python", TOURMALINE_GITHUB,
        "-s", s_file,
        "-r", r_file,
        "-t", t_file,
        "-p", project_id,
        "-a", assay_name,
        "-A", analysis_run_name,
        "-T", TOURMALINE_METADATA,
        "-o", output_file
    ]
    print()
    print(f"Running: {' '.join(cmd)}")
    print()
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in {folder}: {result.stderr}")
        print()
    else:
        print(f"Success for {folder}: {output_file} created.")
        print()

    # ODE input file copy logic
    if OUTPUT_FOLDER_PATH and ODE_INPUTS:
        # 1. Find the most recent repseqs folder
        repseqs_parent = os.path.join(folder_path, f"{folder}_Output")
        repseqs_glob = os.path.join(repseqs_parent, "02_*-repseqs")
        repseqs_dirs = sorted(glob.glob(repseqs_glob))
        if repseqs_dirs:
            latest_repseqs_dir = repseqs_dirs[-1]
            # Copy 02_*-table.tsv from latest repseqs dir
            repseqs_table_files = sorted(glob.glob(os.path.join(latest_repseqs_dir, "02_*-table.tsv")))
            if repseqs_table_files:
                repseqs_table_file = repseqs_table_files[-1]
                shutil.copy2(repseqs_table_file, OUTPUT_FOLDER_PATH)
                print(f"Copied {repseqs_table_file} to {OUTPUT_FOLDER_PATH}")
                # Remove the first row from the copied file
                copied_table_path = os.path.join(OUTPUT_FOLDER_PATH, os.path.basename(repseqs_table_file))
                with open(copied_table_path, 'r') as f:
                    lines = f.readlines()
                with open(copied_table_path, 'w') as f:
                    f.writelines(lines[1:])
            else:
                print(f"No repseqs table file found in {latest_repseqs_dir}")
        else:
            print(f"No repseqs directories found in {repseqs_parent}")

        # 2. Find the most recent taxonomy folder (hyphens)
        taxonomy_folder_glob = os.path.join(
            repseqs_parent,
            f"03_{folder}_{CLASSIFY_METHOD_HYPHEN}_{DATABASE}_*-taxonomy"
        )
        taxonomy_dirs = sorted(glob.glob(taxonomy_folder_glob))
        if taxonomy_dirs:
            latest_taxonomy_dir = taxonomy_dirs[-1]
            # Copy 03_*asv_taxa_features.tsv from latest taxonomy dir (hyphens)
            asv_taxa_files = sorted(glob.glob(os.path.join(latest_taxonomy_dir, f"03_*{CLASSIFY_METHOD_HYPHEN}*asv_taxa_features.tsv")))
            if asv_taxa_files:
                asv_taxa_file = asv_taxa_files[-1]
                shutil.copy2(asv_taxa_file, OUTPUT_FOLDER_PATH)
                print(f"Copied {asv_taxa_file} to {OUTPUT_FOLDER_PATH}")
            else:
                print(f"No asv_taxa_features file found in {latest_taxonomy_dir}")
        else:
            print(f"No taxonomy directories found in {repseqs_parent}") 
