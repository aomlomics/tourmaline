#!/usr/bin/env python

import click
import numpy as np
import pandas as pd

def count_starting_kmers(input_file, file_type, num_seqs, seed):
    """Generate value_counts dataframe of 5' tetramers for random subsample
    of a fasta file"""
    kmer_length = 4
    if file_type == 'fastq':
        interval = 4
    else:
        interval = 2
    if seed:
        np.random.seed(seed)
    starting_kmers = []
    with open(input_file) as handle:
        lines = pd.Series(handle.readlines())
        num_lines = len(lines)
        if num_lines/2 < num_seqs:
            rand_line_nos = np.random.choice(np.arange(1,num_lines,interval), 
                                             size=num_seqs, replace=True)
        else:
            rand_line_nos = np.random.choice(np.arange(1,num_lines,interval), 
                                             size=num_seqs, replace=False)
        rand_lines = lines[rand_line_nos]
    for sequence in rand_lines:
        starting_kmers.append(sequence[:kmer_length])
    starting_kmer_value_counts = pd.Series(starting_kmers).value_counts()
    if (starting_kmer_value_counts.shape[0] == 2):
        starting_kmer_value_counts.loc['NNNX'] = 0
    elif (starting_kmer_value_counts.shape[0] == 1):
        starting_kmer_value_counts.loc['NNNX'] = 0
        starting_kmer_value_counts.loc['NNNY'] = 0
    return(starting_kmer_value_counts)

@click.command()
@click.option('--input_file', '-i', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
              file_okay=True), 
              help="Input sequence file (.fa, .fna, .fasta, .fq, .fastq)")
@click.option('--file_type', '-t', required=False, type=str, default='fasta',
              help="Sequence file type (fasta, fastq) [default: fasta]")
@click.option('--num_seqs', '-n', required=False, type=int, default=10000,
              help="Number of sequences to randomly subsample [default: 10000]")
@click.option('--cutoff', '-c', required=False, type=float, default=0.5,
              help="Minimum fraction of sequences required to match "
                   "a diagnostic 5' tetramer [default: 0.5]")
@click.option('--seed', '-s', required=False, type=int, 
              help="Random number seed [default: None]")

def detect_amplicon_locus(input_file, file_type, num_seqs, cutoff, seed):
    """Determine the most likely amplicon locus of a fasta file based on the
    first four nucleotides.

    The most frequent 5' tetramer in a random subsample of sequences must 
    match, above a given cutoff fraction of sequences, one of the following  
    diagnostic tetramers:

      Tetramer\tAmplicon locus\tForward primer                                  
      TACG\t16S rRNA prok\t515f (EMP)                                           
      GTAG\tITS rRNA euk\tITS1f (EMP)                                           
      GCT[AC]\t18S rRNA euk\tEuk1391f (EMP)                                     
      GTCG\t12S rRNA mt\tMiFish-U-F                                             
    """
    starting_kmer_value_counts = count_starting_kmers(input_file, file_type, 
      num_seqs, seed)
    top_kmer = starting_kmer_value_counts.index[0]
    top_kmer_count = starting_kmer_value_counts[0]
    second_kmer = starting_kmer_value_counts.index[1]
    second_kmer_count = starting_kmer_value_counts[1]
    third_kmer = starting_kmer_value_counts.index[2]
    third_kmer_count = starting_kmer_value_counts[2]
    top_kmer_frac = top_kmer_count/num_seqs
    top2_kmer_frac = (top_kmer_count+second_kmer_count)/num_seqs
    top3_kmer_frac = (top_kmer_count+second_kmer_count+third_kmer_count)/num_seqs
    if (top_kmer == 'TACG') & (top_kmer_frac > cutoff):
        print('Amplicon locus: 16S/515f (%s%% of sequences start with %s)' %
              (round(top_kmer_frac*100, 1), top_kmer))
    elif (top_kmer == 'GTAG') & (top_kmer_frac > cutoff):
        print('Amplicon locus: ITS/ITS1f (%s%% of sequences start with %s)' %
              (round(top_kmer_frac*100, 1), top_kmer))
    elif (top_kmer in ['GCTA', 'GCTC', 'ACAC']) & (second_kmer in ['GCTA', 'GCTC', 
        'ACAC']) & (third_kmer in ['GCTA', 'GCTC', 'ACAC']) & (
        top3_kmer_frac > cutoff):
        print('Amplicon locus: 18S/Euk1391f (%s%% of sequences start with %s, %s, or %s)' %
              (round(top3_kmer_frac*100, 1), top_kmer, second_kmer, third_kmer))
    elif (top_kmer == 'GTCG') & (top_kmer_frac > cutoff):
        print('Amplicon locus: 12S/MiFish-U-F (%s%% of sequences start with %s)' %
              (round(top_kmer_frac*100, 1), top_kmer))
    else:
        print('Could not determine amplicon locus.'),
        print('Most frequent starting tetramer was %s (%s%% of sequences).' %
              (top_kmer, round(top_kmer_frac*100, 1)))

if __name__ == '__main__':
    detect_amplicon_locus()
