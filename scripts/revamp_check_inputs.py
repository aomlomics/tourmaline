"""Check that REVAMP taxonomy inputs all describe the same set of ASVs.

The REVAMP classify method can be handed a BLASTn result file produced on another
machine (`revamp_blast_results`), because the NCBI nt database usually lives on a
different computer than Tourmaline. A stale BLAST file is easy to supply by mistake
and fails silently: every ASV that is missing from it simply comes back Unassigned.
This script fails loudly instead.

Checks performed:
  1. feature IDs in the feature table == sequence IDs in the exported repseqs FASTA
  2. every query ID in the BLAST btab is a known repseqs ID (when a btab is given)
  3. at least one repseqs ID has a BLAST hit (when a btab is given)

ASVs with no BLAST hit at all are legitimate (nothing in nt matched above the
cutoffs), so repseqs IDs missing from the btab are reported, not treated as errors.

Run by Tourmaline during the taxonomy step, in the QIIME 2 environment.

Usage:
  python scripts/revamp_check_inputs.py \
    --repseqs-fasta ASVs.fa \
    --table table.qza \
    --blast-results ASV_blastn_nt.btab \
    --output input_check.txt
"""

import argparse
import sys


def read_fasta_ids(path):
    """Return the sequence IDs of a FASTA file, in file order."""
    ids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_table_ids(path):
    """Return the feature (observation) IDs of a QIIME 2 feature table artifact."""
    import biom
    import qiime2

    table = qiime2.Artifact.load(path).view(biom.Table)
    return list(table.ids('observation'))


def read_btab_query_ids(path):
    """Return the unique query IDs in a BLASTn tabular (outfmt 6) file."""
    ids = set()
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            ids.add(line.split('\t')[0])
    return ids


def preview(items, limit=5):
    """Format a few IDs for an error message."""
    items = sorted(items)
    shown = ', '.join(items[:limit])
    if len(items) > limit:
        shown += f', ... ({len(items)} total)'
    return shown


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--repseqs-fasta', required=True,
                        help='Representative sequences exported to FASTA')
    parser.add_argument('--table', required=True,
                        help='Feature table artifact (.qza)')
    parser.add_argument('--blast-results',
                        help='BLASTn btab file to validate against the repseqs (optional)')
    parser.add_argument('--output', required=True,
                        help='Summary file written when all checks pass')
    args = parser.parse_args()

    errors = []
    notes = []

    repseqs_ids = read_fasta_ids(args.repseqs_fasta)
    repseqs_set = set(repseqs_ids)
    if not repseqs_ids:
        errors.append(f'No sequences found in {args.repseqs_fasta}')
    if len(repseqs_ids) != len(repseqs_set):
        errors.append(f'Duplicate sequence IDs in {args.repseqs_fasta}')
    notes.append(f'repseqs sequences: {len(repseqs_set)}')

    table_set = set(read_table_ids(args.table))
    notes.append(f'feature table features: {len(table_set)}')
    table_only = table_set - repseqs_set
    repseqs_only = repseqs_set - table_set
    if table_only or repseqs_only:
        errors.append(
            'The feature table and the representative sequences describe different ASVs. '
            'They must come from the same repseqs run.\n'
            f'  in table but not repseqs ({len(table_only)}): {preview(table_only)}\n'
            f'  in repseqs but not table ({len(repseqs_only)}): {preview(repseqs_only)}')

    if args.blast_results:
        btab_ids = read_btab_query_ids(args.blast_results)
        notes.append(f'ASVs with BLAST hits: {len(btab_ids & repseqs_set)}')
        unknown = btab_ids - repseqs_set
        if unknown:
            errors.append(
                f'{args.blast_results} contains query IDs that are not in the representative '
                'sequences. The BLAST results were produced from a different set of ASVs.\n'
                f'  unknown query IDs ({len(unknown)}): {preview(unknown)}')
        elif not btab_ids & repseqs_set:
            errors.append(
                f'{args.blast_results} shares no query IDs with the representative sequences.')
        else:
            no_hits = repseqs_set - btab_ids
            if no_hits:
                notes.append(f'ASVs with no BLAST hit (will be Unassigned): {len(no_hits)}')

    if errors:
        sys.exit('ERROR: REVAMP input check failed.\n\n' + '\n\n'.join(errors) + '\n')

    with open(args.output, 'w') as handle:
        for note in notes:
            handle.write(note + '\n')
    print('REVAMP input check passed: ' + '; '.join(notes))


if __name__ == '__main__':
    main()
