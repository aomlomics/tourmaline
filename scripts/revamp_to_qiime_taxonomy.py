"""Convert a REVAMP ASV taxonomy table into a QIIME 2 taxonomy TSV.

REVAMP's `asv_taxonomy_processing_figureOuts.pl` writes
`{run_name}_asvTaxonomyTable.txt`: one row per ASV with seven fixed rank columns
(Kingdom..Species). This script rewrites it as the three-column taxonomy TSV
Tourmaline's other classify methods produce, so the downstream QIIME 2 rules
(import, collapse, barplot, ODE export) work unchanged.

Conversions applied:
  * `NA`, empty and REVAMP's gap-filler ranks (e.g. `Embiotocidae__o`, meaning an
    unnamed order containing family Embiotocidae) all become `NA`.
  * Trailing `NA` ranks are dropped; internal `NA` ranks keep their position, so a
    rank never shifts up into the slot above it.
  * `Unknown` / `Environmental Unknown` assignments, and ASVs absent from the REVAMP
    table, become `Unassigned`.
  * The third column is `percent_id`: the best-hit BLAST percent identity for that
    ASV on a 0-100 scale (REVAMP has no confidence score), or 0 when unassigned.

Run by Tourmaline during the taxonomy step.

Usage:
  python scripts/revamp_to_qiime_taxonomy.py \
    --asv-taxonomy-table run_asvTaxonomyTable.txt \
    --formatted-blast ASV_blastn_nt_formatted.txt \
    --repseqs-fasta ASVs.fa \
    --output run-taxonomy.tsv \
    --taxaranks kingdom,phylum,class,order,family,genus,species
"""

import argparse
import re
import sys

UNASSIGNED = 'Unassigned'
# REVAMP fills a missing intermediate rank with the name below it plus a rank suffix,
# e.g. "Embiotocidae__o" for an unnamed order holding the family Embiotocidae.
GAP_FILLER = re.compile(r'__[kpcofgs]$')
EMPTY_VALUES = {'', 'NA', 'na', 'N/A'}
UNKNOWN_VALUES = {'Unknown', 'Environmental Unknown'}


def read_fasta_ids(path):
    """Return the sequence IDs of a FASTA file, in file order."""
    ids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_revamp_table(path):
    """Return {ASV id: [rank values]} from a REVAMP asvTaxonomyTable.txt."""
    assignments = {}
    with open(path) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        n_ranks = len(header) - 1
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            # Pad short rows so every ASV has the same number of ranks.
            ranks = fields[1:] + [''] * (n_ranks - len(fields) + 1)
            assignments[fields[0]] = ranks
    return assignments, n_ranks


def read_percent_ids(path):
    """Return {ASV id: best-hit percent identity} from ASV_blastn_nt_formatted.txt."""
    percent_ids = {}
    with open(path) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        try:
            asv_col = header.index('ASV')
            percent_col = header.index('percent')
        except ValueError:
            sys.exit(f'ERROR: {path} is missing the ASV and/or percent column. '
                     f'Found: {header}')
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= percent_col:
                continue
            # Several best hits are joined with ", "; they share the same percent
            # identity by construction, so the first value is the best hit.
            value = fields[percent_col].split(',')[0].strip()
            try:
                percent_ids[fields[asv_col]] = float(value)
            except ValueError:
                continue
    return percent_ids


def build_taxon(ranks):
    """Collapse REVAMP rank values into a QIIME taxonomy string."""
    cleaned = []
    for rank in ranks:
        rank = rank.strip()
        if rank in EMPTY_VALUES or rank in UNKNOWN_VALUES or GAP_FILLER.search(rank):
            cleaned.append('NA')
        else:
            cleaned.append(rank)
    # Drop trailing NA ranks only; internal ones keep their position.
    while cleaned and cleaned[-1] == 'NA':
        cleaned.pop()
    if not cleaned:
        return UNASSIGNED
    return ';'.join(cleaned)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--asv-taxonomy-table', required=True,
                        help='REVAMP {run_name}_asvTaxonomyTable.txt')
    parser.add_argument('--formatted-blast', required=True,
                        help='REVAMP ASV_blastn_nt_formatted.txt (for percent identity)')
    parser.add_argument('--repseqs-fasta', required=True,
                        help='Representative sequences exported to FASTA; sets output order')
    parser.add_argument('--output', required=True, help='Taxonomy TSV to write')
    parser.add_argument('--taxaranks', required=True,
                        help='Comma-separated rank names; REVAMP always assigns seven')
    args = parser.parse_args()

    taxa_ranks = [rank.strip() for rank in args.taxaranks.split(',') if rank.strip()]

    assignments, n_ranks = read_revamp_table(args.asv_taxonomy_table)
    if len(taxa_ranks) != n_ranks:
        sys.exit(f'ERROR: taxa_ranks has {len(taxa_ranks)} ranks ({args.taxaranks}) but '
                 f'{args.asv_taxonomy_table} has {n_ranks}. REVAMP assigns seven ranks '
                 '(kingdom,phylum,class,order,family,genus,species); set taxa_ranks to match.')

    percent_ids = read_percent_ids(args.formatted_blast)
    repseqs_ids = read_fasta_ids(args.repseqs_fasta)

    n_unassigned = 0
    with open(args.output, 'w') as out:
        out.write('Feature ID\tTaxon\tpercent_id\n')
        for asv in repseqs_ids:
            if asv in assignments:
                taxon = build_taxon(assignments[asv])
            else:
                taxon = UNASSIGNED
            if taxon == UNASSIGNED:
                n_unassigned += 1
                percent_id = 0
            else:
                percent_id = percent_ids.get(asv, 0)
            out.write(f'{asv}\t{taxon}\t{percent_id}\n')

    print(f'Wrote {args.output}: {len(repseqs_ids)} ASVs, {n_unassigned} unassigned.')


if __name__ == '__main__':
    main()
