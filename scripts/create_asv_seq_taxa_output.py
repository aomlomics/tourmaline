import argparse
import pandas as pd
import numpy as np
import subprocess


def main():
    parser = argparse.ArgumentParser(description="Process taxonomy data.")
    parser.add_argument('--input_repseqs', required=True, help='Path to the input repseqs file')
    parser.add_argument('--input_taxonomy', required=True, help='Path to the input taxonomy file')
    parser.add_argument('--output', required=True, help='Path to the output file')
    parser.add_argument('--taxaranks', required=True, type=str, help='Comma-separated list of taxonomy ranks')

    args = parser.parse_args()

    if args.taxaranks is not None:
        args.taxaranks= [s.strip() for s in args.taxaranks.split(",")]

    # Run QIIME commands
    subprocess.run(
        f"qiime metadata tabulate "
        f"--m-input-file {args.input_repseqs} "
        f"--m-input-file {args.input_taxonomy} "
        f"--o-visualization merged-data.qzv",
        shell=True,
        check=True
    )

    subprocess.run(
        f"qiime tools export "
        f"--input-path merged-data.qzv "
        f"--output-path {args.output}",
        shell=True,
        check=True
    )

    # Move and process the exported metadata
    subprocess.run(f"mv {args.output}/metadata.tsv temp", shell=True, check=True)
    subprocess.run(f"rm -r {args.output}", shell=True, check=True)
    subprocess.run(
        f"sed -e '2d' temp | sed '1 s|id\\t|featureid\\t|' | sed '1 s|Taxon|taxonomy|' | sed '1 s|Sequence|dna_sequence|' > {args.output}",
        shell=True,
        check=True
    )
    subprocess.run(f"/bin/rm -r temp merged-data.qzv", shell=True, check=True)

    # Read the TSV file into pandas
    df = pd.read_csv(args.output, sep='\t')

    # Remove characters from the taxonomy column that match the pattern of an alpha character followed by __
    # and replace "; " with ";"
    df['verbatimIdentification'] = df['taxonomy']
    df['taxonomy'] = df['taxonomy'].str.replace(r'\b[a-zA-Z]__', '', regex=True).str.replace('; ', ';')


    df[args.taxaranks] = [""] * len(args.taxaranks) 
    #taxonomy_split = df['taxonomy'].str.split(';', expand=True)

    tax_ranks = args.taxaranks
    for index, row in df.iterrows():
        taxa = row['taxonomy'].split(";")
        for i in range(0,len(taxa)):
            if i < len(tax_ranks):
                df.loc[index,tax_ranks[i]] = taxa[i]
    
    # replace None with NA
    df = df.fillna(value=np.nan)
    
    # replace _,- with space, remove sp. 
    df[tax_ranks] = df[tax_ranks].replace('_',' ',regex=True)
    df[tax_ranks] = df[tax_ranks].replace(' sp\.','',regex=True)
    df[tax_ranks] = df[tax_ranks].replace(' spp\.','',regex=True)
    df[tax_ranks] = df[tax_ranks].replace('-',' ',regex=True)
    df[tax_ranks] = df[tax_ranks].replace('\/',' ',regex=True)
    

    # # Ensure the number of columns matches the length of params.taxaranks
    # if len(taxonomy_split.columns) == len(args.taxaranks):
    #     taxonomy_split.columns = args.taxaranks
    # else:
    #     raise ValueError("The number of taxonomy ranks does not match the number of columns in the taxonomy data.")

    # Concatenate the new columns with the original DataFrame
    # df = pd.concat([df, taxonomy_split], axis=1)

    # Move the 4th column (confidence or consensus) to the end
    cols = df.columns.tolist()
    cols.append(cols.pop(3))
    df = df[cols]

    # Save the modified DataFrame back to a TSV file
    df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()