#!/usr/bin/env python

import click
import pandas as pd

@click.command()
@click.option('--input', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
              file_okay=True), 
              help="File containing per base sequence quality for each sequence. "
              "It is usually named 'mqc_fastqc_per_base_sequence_quality_plot_1.txt' and "
              "is located in the directory 'multiqc_data'. These are the output of "
              "`multiqc --export $DIR`, where $DIR is the location of the fastqc output for "
              "Read 1 or Read 2 (but not both). The --export option creates these extra files.")
@click.option('--cutoff', '-c', required=False, type=float, default=0.9,
              help="Fraction of maximum median per base sequence quality "
                   "at which point sequence should be truncated [default: 0.9]")

def calculate_fastqc_dropoff(input, cutoff):
    """Determine the position where median per base sequence quality drops below
    some fraction (default: 0.90) of its maximum value and print to standard output.
    This is useful for defining 3' truncation positions in DADA2 (trunc-len) and
    Deblur (trim-length). This script MUST be run separately for Read 1 and Read 2.
    Do this by putting the R1 fastqc output in one directory and the R2 fastqc output
    in another directory, then running `multiqc --export $DIR` inside each directory.
    """
    firstrow = pd.read_csv(input, sep='\t', index_col=0, nrows=1)
    firstcol = pd.read_csv(input, sep='\t', index_col=0, usecols=[0])
    mqc = pd.read_csv(input, sep='\t', index_col=0, usecols=range(firstrow.shape[1]+1), skiprows=range(2,firstcol.shape[0],2))
    mqc_max_to_end = mqc.loc[:, mqc.median().idxmax():]
    positions_below_cutoff = mqc_max_to_end.columns[mqc_max_to_end.median()/mqc_max_to_end.median().max() <= cutoff]
    if len(positions_below_cutoff) > 0:
        print(positions_below_cutoff[0])
    else:
        print(mqc_max_to_end.columns[-1])

if __name__ == '__main__':
    calculate_fastqc_dropoff()
