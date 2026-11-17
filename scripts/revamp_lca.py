"""REVAMP's taxonomy assignment, for a reference taxonomy with any number of ranks.

Used by the `revamp-custom` classify method. REVAMP's own assignment script
(assets/asv_taxonomy_processing_figureOuts.pl) is fixed at seven ranks
(Kingdom..Species) and six cutoffs, so this script reimplements its algorithm for N
ranks and N-1 cutoffs. Hit selection is unchanged: REVAMP's BLAST loop and its
assets/reformat_blast.R still choose the best hits (query-coverage filter, then all
hits at the best percent identity), and this script reads that R script's output.
The `revamp` method (NCBI nt) still uses REVAMP's Perl script.

The algorithm, per ASV, as in REVAMP:
  1. Depth. Compare the best hit's percent identity with the cutoffs (deepest rank
     first). At or above cutoff k (0-based), the ASV may be assigned N-k ranks; below
     the last cutoff it is "Unknown". So the first rank is never assigned on its own,
     exactly as REVAMP never assigns a kingdom alone.
  2. Lineage preparation (REVAMP's fillIn_taxonkit.pl and its ';;' handling):
       - an empty first rank makes the lineage "Environmental Unknown";
       - an empty rank above the deepest named rank among ranks 2..N-1 is filled from
         the rank below it (REVAMP writes e.g. "Embiotocidae__o"; here it is a
         placeholder that compares the same way and is written as NA);
       - any other empty rank with a named rank below it is a gap: the lineage is
         cut at the gap.
  3. Each best hit's lineage is cut to the depth from step 1. A hit is "complete"
     if it has a name (or filled placeholder) at every one of those ranks; otherwise
     it is cut further at its first missing rank.
  4. If any hit is complete, the assignment is the longest common prefix of the
     complete hits; otherwise of all hits. No common first rank gives "Unknown".

Deliberate differences from REVAMP's Perl, all on lineages that lack trailing ranks
(the Perl loses track of rank positions when Perl's split() drops trailing empty
fields):
  * A lineage missing only its last rank (e.g. no species) counts as incomplete at
    full depth, like any other shallow lineage; the Perl treats it as complete.
  * At a reduced depth such a lineage is cut positionally; the Perl removes one
    rank too many (a genus-only reference at genus depth gives family).
  * REVAMP's NCBI-specific rule that "Eukaryota" with nothing below it is
    "Environmental Unknown" is not applied.
  * Taxon names are compared exactly; REVAMP first replaces every character outside
    [A-Za-z0-9 ;] with '_'.
On lineages with every rank named, or with internal gaps only, the output matches
REVAMP's; tests/test_revamp_lca.py checks this against the Perl script.

Outputs:
  --output    a REVAMP-style asvTaxonomyTable: ASV plus one column per rank, NA for
              missing or filled ranks, "Unknown" / "Environmental Unknown" in the first
              column when unassigned. scripts/revamp_to_qiime_taxonomy.py converts it.
  --details   optional per-ASV table: best percent identity, depth, hit counts.

Usage:
  python scripts/revamp_lca.py \
    --formatted-blast revamp/blast_results/ASV_blastn_nt_formatted.txt \
    --pseudo-taxonomy revamp/refdb/pseudo_taxonomy.tsv \
    --asv-fasta revamp/dada2/ASVs.fa \
    --taxaranks domain,kingdom,phylum,class,order,family,genus,species \
    --cutoffs 98,96,92,85,80,75,70 \
    --output revamp/ASV2Taxonomy/run_asvTaxonomyTable.txt
"""

import argparse
import sys

GAP = object()       # a missing rank with a named rank below it
UNKNOWN = 'Unknown'
ENV_UNKNOWN = 'Environmental Unknown'


class Filled:
    """An empty rank filled from the rank below it (REVAMP's 'Name__o').

    Two filled ranks are equal when they were filled from equal values, which is
    what comparing REVAMP's 'Name__o' strings amounts to.
    """
    __slots__ = ('below',)

    def __init__(self, below):
        self.below = below

    def __eq__(self, other):
        return isinstance(other, Filled) and self.below == other.below

    def __hash__(self):
        return hash(('filled', self.below))

    def __repr__(self):
        return f'Filled({self.below!r})'


def parse_cutoffs(text, n_ranks):
    """Return the cutoffs as floats, checking there are N-1 of them, non-increasing."""
    try:
        cutoffs = [float(value) for value in str(text).split(',') if value.strip()]
    except ValueError:
        raise ValueError(f'revamp_taxonomy_cutoffs must be comma-separated numbers: {text}')
    if len(cutoffs) != n_ranks - 1:
        raise ValueError(
            f'revamp_taxonomy_cutoffs has {len(cutoffs)} values ({text}) but taxa_ranks has '
            f'{n_ranks} ranks; revamp-custom needs {n_ranks - 1} cutoffs, one per rank from '
            'the deepest rank up to the second rank (the first rank has none).')
    if any(b > a for a, b in zip(cutoffs, cutoffs[1:])):
        raise ValueError(f'revamp_taxonomy_cutoffs must not increase from left (deepest '
                         f'rank) to right: {text}')
    return cutoffs


def assigned_depth(percent, cutoffs):
    """Number of ranks an ASV may be assigned at this percent identity, or 0."""
    n_ranks = len(cutoffs) + 1
    for k, cutoff in enumerate(cutoffs):
        if percent >= cutoff:
            return n_ranks - k
    return 0


def prepare_lineage(ranks):
    """Turn a list of rank names ('' = missing) into REVAMP's working form.

    Returns a list of length N whose entries are a name, a Filled placeholder, GAP, or
    None for trailing ranks with nothing named below them.
    """
    n = len(ranks)
    if not ranks[0]:
        return [ENV_UNKNOWN] + [GAP] * (n - 1)
    values = list(ranks)
    # Fill upwards from the deepest named rank among ranks 2..N-1 (REVAMP: phylum..genus).
    anchor = max((i for i in range(1, n - 1) if values[i]), default=0)
    for i in range(anchor - 1, 0, -1):
        if not values[i]:
            values[i] = Filled(values[i + 1])
    last = max(i for i in range(n) if values[i])
    prepared = []
    in_gap = False
    for i, value in enumerate(values):
        if i > last:
            prepared.append(None)
        elif in_gap or not value:
            in_gap = True
            prepared.append(GAP)
        else:
            prepared.append(value)
    return prepared


def cut_to_depth(prepared, depth):
    """Return (ranks, complete) for a prepared lineage cut to *depth* ranks."""
    head = prepared[:depth]
    kept = []
    for value in head:
        if value is GAP or value is None:
            return kept, False
        kept.append(value)
    return kept, True


def common_prefix(lineages):
    prefix = list(lineages[0])
    for lineage in lineages[1:]:
        n = 0
        while n < len(prefix) and n < len(lineage) and prefix[n] == lineage[n]:
            n += 1
        del prefix[n:]
    return prefix


def assign(percent, hit_lineages, cutoffs):
    """Assign one ASV. Returns (ranks or None for Unknown, depth, n_complete)."""
    depth = assigned_depth(percent, cutoffs)
    if depth == 0 or not hit_lineages:
        return None, depth, 0
    complete, incomplete = [], []
    for prepared in hit_lineages:
        ranks, is_complete = cut_to_depth(prepared, depth)
        (complete if is_complete else incomplete).append(ranks)
    chosen = complete if complete else incomplete
    prefix = common_prefix(chosen)
    return (prefix or None), depth, len(complete)


def read_pseudo_taxonomy(path, n_ranks):
    """Return {pseudo-taxid: prepared lineage}."""
    lineages = {}
    with open(path) as handle:
        handle.readline()
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            taxid, _, lineage = line.partition('\t')
            ranks = lineage.split(';')
            if len(ranks) != n_ranks:
                sys.exit(f'ERROR: {path} has a {len(ranks)}-rank lineage for pseudo-taxid '
                         f'{taxid} but taxa_ranks has {n_ranks}. Rebuild the reference '
                         'database with the same taxa_ranks.')
            lineages[taxid] = prepare_lineage(ranks)
    return lineages


def read_formatted_blast(path):
    """Return {ASV: (best percent identity, [subject taxids])} from reformat_blast.R output."""
    hits = {}
    with open(path) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        try:
            asv_col, pct_col, taxid_col = (header.index(c) for c in ('ASV', 'percent', 'taxid'))
        except ValueError:
            sys.exit(f'ERROR: {path} is missing an ASV, percent or taxid column: {header}')
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= max(asv_col, pct_col, taxid_col):
                continue
            # reformat_blast.R keeps only hits at the best percent identity and joins
            # unique values with ", ", so the first value is the best hit's.
            percent = float(fields[pct_col].split(',')[0])
            taxids = [t.strip() for t in fields[taxid_col].split(',') if t.strip()]
            hits[fields[asv_col]] = (percent, taxids)
    return hits


def read_fasta_ids(path):
    ids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def render(ranks, n_ranks):
    """Rank values for the output table: names, NA for filled or missing ranks."""
    cells = [value if isinstance(value, str) else 'NA' for value in ranks]
    return cells + ['NA'] * (n_ranks - len(cells))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--formatted-blast', required=True,
                        help='ASV_blastn_nt_formatted.txt from REVAMP reformat_blast.R, with '
                             'pseudo-taxids in the taxid column')
    parser.add_argument('--pseudo-taxonomy', required=True,
                        help='pseudo_taxonomy.tsv from build_revamp_custom_refdb.py')
    parser.add_argument('--asv-fasta', required=True, help='ASV sequences; sets output rows')
    parser.add_argument('--taxaranks', required=True, help='Comma-separated rank names')
    parser.add_argument('--cutoffs', required=True,
                        help='N-1 percent identity cutoffs, deepest rank first')
    parser.add_argument('--output', required=True, help='asvTaxonomyTable to write')
    parser.add_argument('--details', help='Optional per-ASV details table to write')
    args = parser.parse_args()

    taxa_ranks = [rank.strip() for rank in args.taxaranks.split(',') if rank.strip()]
    n_ranks = len(taxa_ranks)
    if n_ranks < 2:
        sys.exit('ERROR: revamp-custom needs at least two ranks in taxa_ranks.')
    try:
        cutoffs = parse_cutoffs(args.cutoffs, n_ranks)
    except ValueError as err:
        sys.exit(f'ERROR: {err}')

    lineages = read_pseudo_taxonomy(args.pseudo_taxonomy, n_ranks)
    hits = read_formatted_blast(args.formatted_blast)
    asv_ids = read_fasta_ids(args.asv_fasta)

    unknown_taxids = set()
    counts = {'assigned': 0, 'no_hit': 0, 'below_cutoffs': 0, 'unknown': 0}
    detail_rows = []
    with open(args.output, 'w') as out:
        out.write('ASV\t' + '\t'.join(rank.capitalize() for rank in taxa_ranks) + '\n')
        for asv in asv_ids:
            percent, depth, n_complete, n_hits = '', 0, 0, 0
            if asv not in hits:
                result, status = None, 'no_hit'
            else:
                percent, taxids = hits[asv]
                hit_lineages = []
                for taxid in taxids:
                    if taxid in lineages:
                        hit_lineages.append(lineages[taxid])
                    else:
                        unknown_taxids.add(taxid)
                n_hits = len(hit_lineages)
                result, depth, n_complete = assign(percent, hit_lineages, cutoffs)
                if result is not None and result[0] != ENV_UNKNOWN:
                    status = 'assigned'
                elif depth == 0:
                    status = 'below_cutoffs'
                else:
                    status = 'unknown'
            counts[status] += 1
            cells = render(result, n_ranks) if result else [UNKNOWN] + ['NA'] * (n_ranks - 1)
            out.write(asv + '\t' + '\t'.join(cells) + '\n')
            depth_label = taxa_ranks[depth - 1] if depth else ''
            detail_rows.append((asv, percent, depth_label, n_hits, n_complete, status,
                                ';'.join(cells)))

    if args.details:
        with open(args.details, 'w') as det:
            det.write('ASV\tbest_percent_id\tmax_depth_rank\tn_best_hits\t'
                      'n_complete_hits\tstatus\ttaxonomy\n')
            for row in detail_rows:
                det.write('\t'.join(str(v) for v in row) + '\n')

    if unknown_taxids:
        example = ', '.join(sorted(unknown_taxids)[:5])
        sys.exit(f'ERROR: {len(unknown_taxids)} taxids in {args.formatted_blast} are not in '
                 f'{args.pseudo_taxonomy} (e.g. {example}). The BLAST results do not '
                 'come from this reference database.')

    print(f'revamp-custom: {len(asv_ids)} ASVs; {counts["assigned"]} assigned, '
          f'{counts["no_hit"]} with no BLAST hit passing query coverage, '
          f'{counts["below_cutoffs"]} below the last cutoff, '
          f'{counts["unknown"]} with no common rank.')


if __name__ == '__main__':
    main()
