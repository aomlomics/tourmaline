import argparse
import yaml
import pandas as pd
import os
import glob
import shutil

## ADD check for repeated run names, place to add project_id, assay_name, user provided terms

def load_yaml(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def find_config_file(working_dir, step_name, run_name):
    """
    Find config file for a given step and run name.
    Looks for files matching pattern: {working_dir}/{run_name}-{step_name}/*.yaml
    """
    step_dir = os.path.join(working_dir, f"{run_name}-{step_name}")
    if not os.path.exists(step_dir):
        raise FileNotFoundError(f"Step directory not found: {step_dir}")
    
    # Look for YAML files in the step directory
    yaml_files = glob.glob(os.path.join(step_dir, "*.yaml"))
    if not yaml_files:
        raise FileNotFoundError(f"No YAML config files found in: {step_dir}")
    
    # If multiple YAML files, prefer the one with the run name
    for yaml_file in yaml_files:
        if run_name in os.path.basename(yaml_file):
            return yaml_file
    
    # If no exact match, return the first YAML file found
    return yaml_files[0]

def find_and_copy_file(working_dir, step_name, run_name, pattern, output_dir, prefix, new_filename):
    """
    Find a file matching the pattern in the step directory and copy it to output directory.
    Returns the path to the copied file.
    """
    step_dir = os.path.join(working_dir, f"{run_name}-{step_name}")
    if not os.path.exists(step_dir):
        raise FileNotFoundError(f"Step directory not found: {step_dir}")
    
    # Look for files matching the pattern
    matching_files = glob.glob(os.path.join(step_dir, pattern))
    if not matching_files:
        raise FileNotFoundError(f"No files matching pattern '{pattern}' found in: {step_dir}")
    
    # Use the first matching file
    source_file = matching_files[0]
    dest_file = os.path.join(output_dir, f"{prefix}_{new_filename}")
    
    # Copy the file
    shutil.copy2(source_file, dest_file)
    print(f"Copied {source_file} to {dest_file}")
    
    return dest_file

def dict_to_tsv(data, file_path):
    df = pd.DataFrame(list(data.items()), columns=['term_name', 'values'])
    df.to_csv(file_path, sep='\t', index=False)

def save_yaml(data, file_path):
    with open(file_path, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)

# FORMATTING FUNCTIONS

def min_reads(repseqs):
    if repseqs['repseq_min_abundance'] > 0:
        return (repseqs['repseq_min_abundance'],"%")
    else:
        return (1,"reads")

def asv_tools(repseqs,tour):
    if repseqs['asv_method'] in ['dada2pe','dada2se']:
        return ";".join([tour['qiime2_version'], "DADA2 "+str(tour['dada2_version'])])
    else:
        return ";".join([tour['qiime2_version'], "deblur "+str(tour['deblur_version'])])
    
def assign_tools(taxa,tour):
    if taxa['classify_method'] == 'consensus-blast':
        software = ";".join([tour['qiime2_version'], "blast "+str(tour['blast_version'])])
        cat = "sequence similarity"
    elif taxa['classify_method'] == 'vsearch':
        software = ";".join([tour['qiime2_version'], "vsearch "+str(tour['vsearch_version'])])
        cat = "sequence similarity"
    elif taxa['classify_method'] == 'naive-bayes':
        software = ";".join([tour['qiime2_version'], "naive-bayes classifier; scikit-learn "+str(tour['scikit-learn_version'])])
        cat = "sequence composition"
    return (software,cat)

def trim_paramF(qaqc,repseqs,tour):
    output=""
    if qaqc['to_trim']:
        software = "; ".join([tour['qiime2_version'], "Cutadapt "+str(tour['cutadapt_version'])])
        if qaqc['paired_end']:
            output+=f"Trim forward reads of reverse complement of reverse primer, and reverse reads of reverse complement of forward primer. "
            #+=f"qiime cutadapt trim-paired --p-adapter-f {revcomp_primerR} --p-adapter-r {revcomp_primerF} --p-match-read-wildcards --p-match-adapter-wildcards --p-minimum-length {qaqc['minimum_length']} "
            if qaqc['discard_untrimmed']:
                output+=f"Then trim forward reads of forward primer, and reverse reads of reverse primer, discarding untrimmed reads. Minimum length of {str(qaqc['minimum_length'])}. "
                return (software,output)
            else:
                output+=f"Then trim forward reads of forward primer, and reverse reads of reverse primer. Minimum length of {str(qaqc['minimum_length'])} bp. "
                return (software,output)
        else:
            output+=f"Trim reverse complement of reverse primer. "
            if qaqc['discard_untrimmed']:
                output+=f"Then trim reads of forward primer, discarding untrimmed reads. Minimum length of {qaqc['minimum_length']} bp. "
                return (software,output)
            else:
                output+=f"Then trim reads of forward primer. Minimum length of {qaqc['minimum_length']}. "
                return (software,output)
    elif repseqs['asv_method'] == 'dada2pe' and (repseqs['dada2_trim_left_f'] > 0 or repseqs['dada2pe_trim_left_r'] > 0):
        software = asv_tools(repseqs,tour)
        output+=f"Trim {str(repseqs['dada2_trim_left_f'])} bp from forward reads, and {str(repseqs['dada2pe_trim_left_r'])} bp from reverse reads. "
        return (software,output)
    elif repseqs['asv_method'] == 'dada2se' and (repseqs['dada2_trim_left_f'] > 0):
        software = asv_tools(repseqs,tour)
        output+=f"Trim {str(repseqs['dada2_trim_left_f'])} bp from reads"
        return (software,output)
    else:
        return ("not applicable","not applicable")

def assign_collapse(taxa):
    if taxa['classify_method'] in ['consensus-blast','consensus-vsearch']:
        output = f"minimum consensus lowest common ancestor of {str(taxa['min_consensus'])}"
    elif taxa['classify_method'] == 'naive-bayes':
        output = f"confidence threshold of {str(taxa['skl_confidence'])}"
    elif taxa['classify_method'] == 'bt2-blca':
        output = f"minimum confidence threshold of {str(taxa['confidence_threshold'])}"
    return output



def main():
    parser = argparse.ArgumentParser(description="Generate a single TSV file from multiple YAML files.")
    parser.add_argument('-w','--working_dir', required=True, help='Working directory containing step folders')
    parser.add_argument('-s','--qaqc_run_name', required=True, help='Run name for qaqc step')
    parser.add_argument('-r','--repseqs_run_name', required=True, help='Run name for repseqs step')
    parser.add_argument('-t','--taxonomy_run_name', required=True, help='Run name for taxonomy step')
    parser.add_argument('-p','--project_id', required=True, help='Value for project_id')
    parser.add_argument('-a','--assay_name', help='Value for assay_name, otherwise use value in qaqc config')
    parser.add_argument('-A','--analysis_run_name', help='Value for analysis_run_name, otherwise use taxonomy run name')
    parser.add_argument('-T','--tourmaline_metadata',default="./00-data/tourmaline_metadata.yaml", help='Path to tourmaline metadata')
    parser.add_argument('-O','--output_folder', required=True, help='Output folder path where files will be saved')
    #parser.add_argument('--checklist', required=True, help='Path to the CSV file with metadata terms')

    args = parser.parse_args()

    # Validate working directory
    if not os.path.exists(args.working_dir):
        print(f"Error: Working directory does not exist: {args.working_dir}")
        return 1

    # Create output directory if it doesn't exist
    os.makedirs(args.output_folder, exist_ok=True)
    print(f"Output folder: {args.output_folder}")

    # Find and load the YAML config files
    try:
        qaqc_config_path = find_config_file(args.working_dir, "qaqc", args.qaqc_run_name)
        repseqs_config_path = find_config_file(args.working_dir, "repseqs", args.repseqs_run_name)
        taxonomy_config_path = find_config_file(args.working_dir, "taxonomy", args.taxonomy_run_name)
        
        print(f"Found qaqc config: {qaqc_config_path}")
        print(f"Found repseqs config: {repseqs_config_path}")
        print(f"Found taxonomy config: {taxonomy_config_path}")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error finding config files: {e}")
        return 1

    # Load the YAML files
    qaqc1 = load_yaml(qaqc_config_path)
    repseqs2 = load_yaml(repseqs_config_path)
    taxa3 = load_yaml(taxonomy_config_path)
    tour = load_yaml(args.tourmaline_metadata)
    project_id = args.project_id
    assay_name = args.assay_name if args.assay_name else qaqc1['amplicon_name']
    analysis_run_name = args.analysis_run_name if args.analysis_run_name else args.taxonomy_run_name

    # MAPPINGS
    mappings = {
        # FAIR eDNA TERMS
        'project_id': project_id,
        'assay_name': assay_name,
        'analysis_run_name': analysis_run_name,
        "sop_bioinformatics": tour['sop_bioinformatics'],
        "trim_method": trim_paramF(qaqc1,repseqs2,tour)[0],
        "trim_param": trim_paramF(qaqc1,repseqs2,tour)[1],
        "demux_tool": "",
        "demux_max_mismatch": "",
        "merge_tool": asv_tools(repseqs2,tour),
        "merge_min_overlap": 12 if repseqs2['asv_method'] == 'dada2pe' else "not applicable",
        "min_len_cutoff": qaqc1['minimum_length'], # CHECK!
        "min_len_tool": "Cutadapt "+str(tour['cutadapt_version']) if qaqc1['to_trim'] else "not applicable",
        "error_rate_tool": asv_tools(repseqs2,tour),
        "error_rate_cutoff": repseqs2['dada2_max_ee_f'],
        "error_rate_type": "expected error rate",
        "chimera_check_method": "denovo; "+asv_tools(repseqs2,tour),
        "chimera_check_param": "--chimera_method "+repseqs2['dada2_chimera_method']+" --min_parental_fold "+str(repseqs2['dada2_min_fold_parent_over_abundance']) if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "default",
        "otu_clust_tool": asv_tools(repseqs2,tour),
        "otu_clust_cutoff": 100,
        "min_reads_cutoff": min_reads(repseqs2)[0],
        "min_reads_cutoff_unit": min_reads(repseqs2)[1],
        "min_reads_tool": asv_tools(repseqs2,tour),
        "otu_db": "custom",
        "otu_db_custom": taxa3['database_name'],
        "tax_assign_cat": assign_tools(taxa3,tour)[1],
        "otu_seq_comp_appr": assign_tools(taxa3,tour)[0],
        "tax_class_id_cutoff": taxa3['perc_identity'] if taxa3['classify_method'] in ['consensus-blast','consensus-vsearch'] else "not applicable",
        "tax_class_query_cutoff": taxa3['query_cov'] if taxa3['classify_method'] in ['consensus-blast','consensus-vsearch'] else "not applicable",
        "tax_class_other": taxa3['classify_params'],
        "tax_class_collapse": assign_collapse(taxa3),
        "tax_class_other": "",
        "screen_contam_method": "",
        "screen_geograph_method": "",
        "screen_nontarget_method": "",
        "screen_other": "",
        "otu_raw_description": "",
        "otu_final_description": "",
        "bioinfo_method_additional": "",


        # CUSTOM TERMS
        #"analysis_run_name": " | ".join([qaqc1['run_name'],repseqs2['run_name'],taxa3['run_name']]),
        "discard_untrimmed": qaqc1['discard_untrimmed'],
        "qiime2_version": tour['qiime2_version'],
        "tourmaline_asv_method": repseqs2['asv_method'],
        "dada2_trunc_len_f": repseqs2['dada2_trunc_len_f'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "dada2pe_trunc_len_r": repseqs2['dada2pe_trunc_len_r'] if repseqs2['asv_method'] == 'dada2pe' else "not applicable",
        "dada2_trim_left_f": repseqs2['dada2_trim_left_f'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "dada2pe_trim_left_r": repseqs2['dada2pe_trim_left_r'] if repseqs2['asv_method'] == 'dada2pe' else "not applicable",
        "dada2_max_ee_f": repseqs2['dada2_max_ee_f'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "dada2pe_max_ee_r": repseqs2['dada2pe_max_ee_r'] if repseqs2['asv_method'] == 'dada2pe' else "not applicable",
        "dada2_trunc_q": repseqs2['dada2_trunc_q'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "dada2_pooling_method": repseqs2['dada2_pooling_method'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "dada2_chimera_method": repseqs2['dada2_chimera_method'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "dada2_min_fold_parent_over_abundance": repseqs2['dada2_min_fold_parent_over_abundance'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "dada2_n_reads_learn": repseqs2['dada2_n_reads_learn'] if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "not applicable",
        "deblur_trim_length": repseqs2['deblur_trim_length'] if repseqs2['asv_method'] == 'deblur' else "not applicable",
        ## do rest of deblur terms later
        "repseqs_min_abundance": repseqs2['repseq_min_abundance'] if repseqs2['to_filter'] else 0,
        "repseqs_min_length": repseqs2['repseq_min_length'] if repseqs2['to_filter'] else 0,
        "repseqs_max_length": repseqs2['repseq_max_length'] if repseqs2['to_filter'] else 0,
        "repseqs_min_prevalence": repseqs2['repseq_min_prevalence'] if repseqs2['to_filter'] else 0,
        "skl_confidence": taxa3['skl_confidence'] if taxa3['classify_method'] == 'naive-bayes' else "not applicable",
        "min_consensus": taxa3['min_consensus'] if taxa3['classify_method'] in ['consensus-blast','consensus-vsearch'] else "not applicable",
        "confidence_threshold": taxa3['confidence_threshold'] if taxa3['classify_method'] == 'bt2-blca' else "not applicable",

    }

    


    # Load the CSV file
    #metadata_terms = pd.read_csv(args.hecklist)

  

    # Save the combined data to the output TSV file
    try:
        # Generate output filename with analysis_run_name prefix
        metadata_filename = f"{analysis_run_name}_metadata.tsv"
        metadata_path = os.path.join(args.output_folder, metadata_filename)
        dict_to_tsv(mappings, metadata_path)
        print(f"Successfully generated metadata file: {metadata_path}")
        
        # Copy additional files
        print("\nCopying additional files...")
        
        # Copy asv_taxa_features.tsv from taxonomy step
        try:
            asv_taxa_file = find_and_copy_file(
                args.working_dir, 
                "taxonomy", 
                args.taxonomy_run_name, 
                "*asv_taxa_features.tsv", 
                args.output_folder, 
                analysis_run_name,
                "asv_taxa_features.tsv"
            )
        except FileNotFoundError as e:
            print(f"Warning: Could not find asv_taxa_features.tsv file: {e}")
        
        # Copy table.tsv from repseqs step
        try:
            table_file = find_and_copy_file(
                args.working_dir, 
                "repseqs", 
                args.repseqs_run_name, 
                "*-table.tsv", 
                args.output_folder, 
                analysis_run_name,
                "table.tsv"
            )
        except FileNotFoundError as e:
            print(f"Warning: Could not find table.tsv file: {e}")
        
        print(f"\nAll files saved to: {args.output_folder}")
        return 0
        
    except Exception as e:
        print(f"Error writing output file: {e}")
        return 1

if __name__ == "__main__":
    exit(main())