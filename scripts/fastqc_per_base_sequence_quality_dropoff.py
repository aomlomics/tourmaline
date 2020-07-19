#!/usr/bin/env python

import click
import pandas as pd

@click.command()
@click.option('--input_file', '-f', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
              file_okay=True), 
              help="Output of `multiqc --export $DIR` containing per base sequence "
              "quality for each sequence, where $DIR is the location of the fastqc output for "
              "Read 1 or Read 2. File is usually in directory 'multiqc_data' and named "
              "'mqc_fastqc_per_base_sequence_quality_plot_1.txt'")
@click.option('--cutoff', '-c', required=False, type=float, default=0.9,
              help="Minimum fraction of sequences required to match "
                   "a diagnostic 5' tetramer [default: 0.5]")

def calculate_fastqc_dropoff(input_file, cutoff):
    """Determine the position where median per base sequence quality drops below
    some fraction (default: 0.90) of its maximum value. This is useful for 
    defining 3' truncation positions in DADA2 (trunc-len) and Deblur (trim-length).
    This script should be run separately for Read 1 and Read 2 fastqc/multiqc output.
    Prints to standard output.
    """
    mqc = pd.read_csv(input_file, sep='\t', index_col=0)
    mqc_max_to_end = mqc.loc[:, mqc.median().idxmax():]
    position = mqc_max_to_end.columns[mqc_max_to_end.median()/mqc_max_to_end.median().max() <= cutoff][0]
    print(position)

if __name__ == '__main__':
    calculate_fastqc_dropoff()
