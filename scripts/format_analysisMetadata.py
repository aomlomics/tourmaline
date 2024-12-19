import argparse
import yaml
import pandas as pd

## ADD check for repeated run names, place to add project_id, assay_name, user provided terms

def load_yaml(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def dict_to_tsv(data, file_path):
    df = pd.DataFrame(list(data.items()), columns=['field_name', 'values'])
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

def trim_paramF(samples,repseqs,tour):
    output=""
    if samples['to_trim']:
        software = "; ".join([tour['qiime2_version'], "Cutadapt "+str(tour['cutadapt_version'])])
        if samples['paired_end']:
            output+=f"Trim forward reads of reverse complement of reverse primer, and reverse reads of reverse complement of forward primer. "
            #+=f"qiime cutadapt trim-paired --p-adapter-f {revcomp_primerR} --p-adapter-r {revcomp_primerF} --p-match-read-wildcards --p-match-adapter-wildcards --p-minimum-length {samples['minimum_length']} "
            if samples['discard_untrimmed']:
                output+=f"Then trim forward reads of forward primer, and reverse reads of reverse primer, discarding untrimmed reads. Minimum length of {str(samples['minimum_length'])}. "
                return (software,output)
            else:
                output+=f"Then trim forward reads of forward primer, and reverse reads of reverse primer. Minimum length of {str(samples['minimum_length'])} bp. "
                return (software,output)
        else:
            output+=f"Trim reverse complement of reverse primer. "
            if samples['discard_untrimmed']:
                output+=f"Then trim reads of forward primer, discarding untrimmed reads. Minimum length of {samples['minimum_length']} bp. "
                return (software,output)
            else:
                output+=f"Then trim reads of forward primer. Minimum length of {samples['minimum_length']}. "
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
        min_con = f"minimum consensus of {str(taxa['min_consensus'])}"
    elif taxa['classify_method'] == 'naive-bayes':
        min_con = f"confidence threshold of {str(taxa['skl_confidence'])}"
    level = taxa['collapse_taxalevel']
    output = f"collapse to {taxa['taxa_ranks'][level-1]} level with {min_con}"
    return output



def main():
    parser = argparse.ArgumentParser(description="Generate a single TSV file from multiple YAML files.")
    parser.add_argument('-s','--samples_config', required=True, help='Path to the first YAML file')
    parser.add_argument('-r','--repseqs_config', required=True, help='Path to the second YAML file')
    parser.add_argument('-t','--taxonomy_config', required=True, help='Path to the third YAML file')
    parser.add_argument('-T','--tourmaline_metadata',default="./00-data/tourmaline_metadata.yaml", help='Path to tourmaline metadata')
    #parser.add_argument('--checklist', required=True, help='Path to the CSV file with metadata terms')
    parser.add_argument('-o','--output', required=True, help='Path to the output YAML file')

    args = parser.parse_args()

    # Load the YAML files
    samples1 = load_yaml(args.samples_config)
    repseqs2 = load_yaml(args.repseqs_config)
    taxa3 = load_yaml(args.taxonomy_config)
    tour = load_yaml(args.tourmaline_metadata)

    # MAPPINGS
    mappings = {
        # FAIR eDNA TERMS
        'project_id': "",
        "sop_bioinformatics": tour['sop_bioinformatics'],
        "trim_method": trim_paramF(samples1,repseqs2,tour)[0],
        "trim_param": trim_paramF(samples1,repseqs2,tour)[1],
        "min_reads_cutoff": min_reads(repseqs2)[0],
        "min_reads_cutoff_unit": min_reads(repseqs2)[1],
        "min_reads_tool": asv_tools(repseqs2,tour),
        "error_rate_tool": asv_tools(repseqs2,tour),
        "error_rate_cutoff": repseqs2['dada2_max_ee_f'],
        "error_rate_type": "expected error rate",
        "merge_min_overlap": 12 if repseqs2['asv_method'] == 'dada2pe' else "not applicable",
        "merge_tool": asv_tools(repseqs2,tour),
        "otu_clust_cutoff": 100,
        "otu_clust_tool": asv_tools(repseqs2,tour),
        "chimera_check_method": "denovo; "+asv_tools(repseqs2,tour),
        "chimera_check_param": "--chimera_method "+repseqs2['dada2_chimera_method']+" --min_parental_fold "+str(repseqs2['dada2_min_fold_parent_over_abundance']) if repseqs2['asv_method'] in ['dada2pe','dada2se'] else "default",
        "tax_asign_cat": assign_tools(taxa3,tour)[1],
        "otu_seq_comp_appr": assign_tools(taxa3,tour)[0],
        "otu_db": taxa3['database_name'],
        "tax_class_id_cutoff": taxa3['perc_identity'] if taxa3['classify_method'] in ['consensus-blast','consensus-vsearch'] else "not applicable",
        "tax_class_query_cutoff": taxa3['query_cov'] if taxa3['classify_method'] in ['consensus-blast','consensus-vsearch'] else "not applicable",
        "tax_class_other": taxa3['classify_params'],
        "tax_class_collapse": assign_collapse(taxa3),

        # CUSTOM TERMS
        "analysis_run_name": " | ".join([samples1['run_name'],repseqs2['run_name'],taxa3['run_name']]),
        "discard_untrimmed": samples1['discard_untrimmed'],
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

    }

    


    # Load the CSV file
    #metadata_terms = pd.read_csv(args.hecklist)

  

    # Save the combined data to the output YAML file
    dict_to_tsv(mappings, args.output)

if __name__ == "__main__":
    main()