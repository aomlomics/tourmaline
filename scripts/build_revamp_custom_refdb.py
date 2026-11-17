"""Prepare a user-supplied reference database for the `revamp-custom` classify method.

REVAMP normally takes taxonomy from NCBI taxids: BLAST reports a taxid per hit and
taxonkit turns it into a lineage. `revamp-custom` uses the reference database's own
taxonomy instead. Every unique lineage in the reference taxonomy gets an integer
"pseudo-taxid"; after BLAST, each hit's subject ID is replaced by the pseudo-taxid of
its lineage, and `scripts/revamp_lca.py` looks lineages up from the table written here.
No NCBI taxonomy is involved, and the lineage keeps every rank in `taxa_ranks`.

Outputs (in --outdir):
  ref_blast.fasta         reference sequences renamed r1..rN, so makeblastdb
                          -parse_seqids never rejects or rewrites an original ID
  seqid_map.tsv           blast_id, original_id, pseudo_taxid
  pseudo_taxonomy.tsv     pseudo_taxid, lineage (one field per rank, ';'-joined,
                          empty string for a missing rank)
  refdb_report.txt        counts, rank-depth histogram and warnings

Normalization applied to each lineage:
  * Placeholder values ('', NA, Unknown, unidentified, ...; and an empty rank prefix
    such as 'g__') become empty ranks. They would otherwise be compared as real taxon
    names, so two unrelated 'NA' genera would count as the same genus.
  * A lineage with fewer fields than taxa_ranks is padded with empty trailing ranks.
    A lineage with more fields is an error.
  * Taxon names are otherwise kept exactly as given (rank prefixes included).

Run by Tourmaline during the taxonomy step (classify_method: revamp-custom).

Usage:
  python scripts/build_revamp_custom_refdb.py \
    --ref-seqs refseqs.fasta --ref-taxa reftaxa.tsv \
    --taxaranks domain,kingdom,phylum,class,order,family,genus,species \
    --outdir results/run-taxonomy/revamp/refdb
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

# Whole-field values treated as "no name at this rank" (compared case-insensitively).
PLACEHOLDERS = {'', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', 'unidentified',
                'unassigned', 'unclassified'}
# An empty rank prefix such as "g__" or "s__" (Greengenes/SILVA style).
EMPTY_PREFIX = re.compile(r'^[A-Za-z]{1,2}__$')
HEADER_IDS = {'feature id', 'featureid', 'feature-id', '#otu id', 'id', 'seqid',
              'sequence id', '#feature id'}
MAX_EXAMPLES = 20


def is_placeholder(value):
    return value.lower() in PLACEHOLDERS or bool(EMPTY_PREFIX.match(value))


def read_fasta(path):
    """Return an ordered list of (id, sequence); the ID is the first word of the header."""
    records = []
    seen = set()
    current_id = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    records.append((current_id, ''.join(chunks)))
                fields = line[1:].split()
                if not fields:
                    sys.exit(f'ERROR: {path} has a FASTA header with no sequence ID.')
                current_id = fields[0]
                if current_id in seen:
                    sys.exit(f'ERROR: duplicate sequence ID in {path}: {current_id}')
                seen.add(current_id)
                chunks = []
            else:
                if current_id is None:
                    sys.exit(f'ERROR: {path} does not start with a FASTA header.')
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, ''.join(chunks)))
    return records


def read_taxonomy(path, n_ranks):
    """Return ({seq id: lineage tuple}, stats) from a two-column taxonomy TSV.

    A header line is detected and skipped. Columns after the second (for example a
    QIIME 2 Confidence column) are ignored.
    """
    lineages = {}
    stats = Counter()
    placeholder_values = Counter()
    with open(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.rstrip('\n').rstrip('\r')
            if not line.strip():
                continue
            fields = line.split('\t')
            if line_no == 1 and fields[0].strip().lower() in HEADER_IDS:
                continue
            if len(fields) < 2:
                sys.exit(f'ERROR: {path} line {line_no} has no tab-separated taxonomy '
                         f'column: {line[:100]}')
            seq_id = fields[0].strip()
            ranks = [rank.strip() for rank in fields[1].split(';')]
            # A trailing ';' is common and means nothing.
            while len(ranks) > n_ranks and ranks[-1] == '':
                ranks.pop()
            if len(ranks) > n_ranks:
                sys.exit(f'ERROR: {path} line {line_no} ({seq_id}) has {len(ranks)} ranks '
                         f'but taxa_ranks lists {n_ranks}:\n  {fields[1]}\n'
                         'Set taxa_ranks to every rank in the reference taxonomy.')
            if len(ranks) < n_ranks:
                stats['padded'] += 1
                ranks += [''] * (n_ranks - len(ranks))
            cleaned = []
            for rank in ranks:
                if is_placeholder(rank):
                    if rank:
                        placeholder_values[rank] += 1
                    cleaned.append('')
                else:
                    cleaned.append(rank)
            lineage = tuple(cleaned)
            if seq_id in lineages and lineages[seq_id] != lineage:
                sys.exit(f'ERROR: {path} lists {seq_id} more than once with different '
                         'lineages.')
            if seq_id in lineages:
                stats['duplicate_identical'] += 1
            lineages[seq_id] = lineage
    return lineages, stats, placeholder_values


def deepest_named(lineage):
    """Index of the deepest non-empty rank, or -1."""
    for i in range(len(lineage) - 1, -1, -1):
        if lineage[i]:
            return i
    return -1


def unfillable_gaps(lineage):
    """Rank indices of internal gaps that revamp_lca.py cannot fill.

    revamp_lca.py fills an empty rank from the nearest named rank below it, as
    REVAMP's fillIn_taxonkit.pl does, but only from ranks above the deepest one.
    A gap with only the deepest rank named below it (e.g. genus missing, species
    present) stays a gap and truncates the lineage there.
    """
    n = len(lineage)
    last = deepest_named(lineage)
    anchor = max((i for i in range(1, n - 1) if lineage[i]), default=0)
    return [i for i in range(1, last) if not lineage[i] and i > anchor]


def find_conflicts(lineages, taxa_ranks):
    """Names that sit under more than one parent lineage at the same rank."""
    parents = defaultdict(set)
    for lineage in set(lineages):
        for i in range(1, len(lineage)):
            if lineage[i]:
                parents[(i, lineage[i])].add(lineage[:i])
    conflicts = []
    for (i, name), parent_set in sorted(parents.items()):
        if len(parent_set) > 1:
            shown = [';'.join(p) for p in sorted(parent_set)[:3]]
            conflicts.append((taxa_ranks[i], name, len(parent_set), shown))
    return conflicts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ref-seqs', required=True, help='Reference sequences (FASTA)')
    parser.add_argument('--ref-taxa', required=True,
                        help='Reference taxonomy TSV: sequence ID, ;-separated lineage')
    parser.add_argument('--taxaranks', required=True,
                        help='Comma-separated rank names, one per lineage field')
    parser.add_argument('--outdir', required=True, help='Directory for the outputs')
    args = parser.parse_args()

    taxa_ranks = [rank.strip() for rank in args.taxaranks.split(',') if rank.strip()]
    n_ranks = len(taxa_ranks)
    if n_ranks < 2:
        sys.exit('ERROR: revamp-custom needs at least two ranks in taxa_ranks.')

    records = read_fasta(args.ref_seqs)
    lineages, tax_stats, placeholder_values = read_taxonomy(args.ref_taxa, n_ranks)
    if not records:
        sys.exit(f'ERROR: no sequences in {args.ref_seqs}')

    os.makedirs(args.outdir, exist_ok=True)
    pseudo_ids = {}
    kept_lineages = []
    no_taxonomy = []
    with open(os.path.join(args.outdir, 'ref_blast.fasta'), 'w') as fasta_out, \
            open(os.path.join(args.outdir, 'seqid_map.tsv'), 'w') as map_out:
        map_out.write('blast_id\toriginal_id\tpseudo_taxid\n')
        n_written = 0
        for seq_id, seq in records:
            if seq_id not in lineages:
                no_taxonomy.append(seq_id)
                continue
            lineage = lineages[seq_id]
            if lineage not in pseudo_ids:
                pseudo_ids[lineage] = len(pseudo_ids) + 1
            n_written += 1
            blast_id = f'r{n_written}'
            fasta_out.write(f'>{blast_id}\n{seq}\n')
            map_out.write(f'{blast_id}\t{seq_id}\t{pseudo_ids[lineage]}\n')
            kept_lineages.append(lineage)

    if n_written == 0:
        sys.exit(f'ERROR: none of the {len(records)} sequence IDs in {args.ref_seqs} '
                 f'appear in {args.ref_taxa}.')

    with open(os.path.join(args.outdir, 'pseudo_taxonomy.tsv'), 'w') as tax_out:
        tax_out.write('pseudo_taxid\tlineage\n')
        for lineage, pseudo_id in pseudo_ids.items():
            tax_out.write(f'{pseudo_id}\t{";".join(lineage)}\n')

    seq_ids = {seq_id for seq_id, _ in records}
    taxonomy_only = [seq_id for seq_id in lineages if seq_id not in seq_ids]
    depth = Counter(deepest_named(lineage) for lineage in kept_lineages)
    empty_first = sum(1 for lineage in kept_lineages if not lineage[0])
    gap_ranks = Counter()
    for lineage in kept_lineages:
        for i in unfillable_gaps(lineage):
            gap_ranks[taxa_ranks[i]] += 1
    conflicts = find_conflicts(kept_lineages, taxa_ranks)

    lines = [
        'revamp-custom reference database',
        f'  reference sequences:   {args.ref_seqs}',
        f'  reference taxonomy:    {args.ref_taxa}',
        f'  taxa_ranks ({n_ranks}):        {",".join(taxa_ranks)}',
        '',
        f'Sequences used:                 {n_written}',
        f'Unique lineages (pseudo-taxids): {len(pseudo_ids)}',
        f'Sequences with no taxonomy (skipped): {len(no_taxonomy)}',
        f'Taxonomy entries with no sequence (ignored): {len(taxonomy_only)}',
        f'Lineages padded to {n_ranks} ranks: {tax_stats["padded"]}',
        '',
        'Deepest named rank per sequence:',
    ]
    for i in range(n_ranks - 1, -2, -1):
        if depth[i]:
            label = taxa_ranks[i] if i >= 0 else '(no rank named)'
            lines.append(f'  {label:<20} {depth[i]}')
    if placeholder_values:
        lines += ['', 'Placeholder values treated as empty ranks:']
        for value, count in placeholder_values.most_common():
            lines.append(f'  {value!r:<20} {count}')

    warnings = []
    if no_taxonomy:
        warnings.append(f'{len(no_taxonomy)} reference sequences have no taxonomy entry and '
                        f'were left out of the database, e.g. {", ".join(no_taxonomy[:5])}')
    if empty_first:
        warnings.append(f'{empty_first} sequences have no {taxa_ranks[0]}; as in REVAMP, hits '
                        'to them are treated as "Environmental Unknown" and never outweigh a '
                        'hit with a full lineage.')
    if gap_ranks:
        detail = ', '.join(f'{rank}: {count}' for rank, count in gap_ranks.items())
        warnings.append('Some lineages have a missing rank that cannot be filled from the '
                        f'ranks below it ({detail}); hits to them are truncated above the '
                        'gap.')
    if conflicts:
        warnings.append(f'{len(conflicts)} taxon names occur under more than one parent '
                        'lineage. ASVs hitting more than one of them are assigned above the '
                        'conflicting rank. Examples:')
        for rank, name, n_parents, shown in conflicts[:MAX_EXAMPLES]:
            warnings.append(f'    {rank} {name!r} under {n_parents} parents, e.g. '
                            + ' | '.join(shown))
    if warnings:
        lines += ['', 'WARNINGS:'] + [f'  {w}' if w.startswith('    ') else f'  - {w}'
                                      for w in warnings]

    report = '\n'.join(lines) + '\n'
    with open(os.path.join(args.outdir, 'refdb_report.txt'), 'w') as report_out:
        report_out.write(report)
    print(report)


if __name__ == '__main__':
    main()
