"""Build Krona text inputs from a Tourmaline taxonomy table and feature table.

Krona's `ktImportText` takes a tab-separated file per dataset, one line per lineage:
`count<TAB>rank1<TAB>rank2<TAB>...`. This script writes one such file per sample, plus
one summing all samples, and a manifest the krona_plot rule turns into `file,label`
arguments.

Works with every Tourmaline classify method, because it reads only the `Feature ID` and
`Taxon` columns that all of them produce. Rank prefixes (`k__`, `d__`) are stripped,
trailing empty/`NA` ranks are dropped, internal ones are kept so ranks don't shift, and
features with no assignment (or absent from the taxonomy) collapse into `Unassigned`.

Run by Tourmaline during the taxonomy step, in the QIIME 2 environment.

Usage:
  python scripts/taxonomy_to_krona.py \
    --table table.qza \
    --taxonomy run-taxonomy.qza \
    --outdir figures/krona_inputs \
    --manifest figures/krona_inputs/krona_datasets.tsv
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# Matches QIIME-style rank prefixes such as k__, d__, p__ at the start of a rank.
RANK_PREFIX = re.compile(r'^[a-zA-Z]__')
EMPTY_VALUES = {'', 'NA', 'na', 'N/A', 'Unassigned', 'Unclassified', 'unassigned'}
UNASSIGNED = 'Unassigned'
# Krona labels become shell arguments, so keep them free of spaces and commas.
UNSAFE_CHARS = re.compile(r'[^A-Za-z0-9._-]')


def load_table(path):
    """Return a biom.Table from a .qza artifact or a .biom file."""
    import biom

    if path.endswith('.qza'):
        import qiime2
        return qiime2.Artifact.load(path).view(biom.Table)
    with open(path, 'rb') as handle:
        return biom.parse.parse_biom_table(handle)


def load_taxonomy(path):
    """Return {feature id: taxon string} from a .qza artifact or a taxonomy TSV."""
    import pandas as pd

    if path.endswith('.qza'):
        import qiime2
        frame = qiime2.Artifact.load(path).view(pd.DataFrame)
    else:
        frame = pd.read_csv(path, sep='\t', index_col=0)
    column = 'Taxon' if 'Taxon' in frame.columns else frame.columns[0]
    return {str(idx): str(value) for idx, value in frame[column].items()}


def split_taxon(taxon):
    """Split a taxonomy string into the rank list Krona should nest."""
    if taxon is None:
        return [UNASSIGNED]
    ranks = []
    for rank in str(taxon).split(';'):
        rank = RANK_PREFIX.sub('', rank.strip())
        ranks.append('NA' if rank in EMPTY_VALUES else rank)
    while ranks and ranks[-1] == 'NA':
        ranks.pop()
    if not ranks:
        return [UNASSIGNED]
    return ranks


def write_dataset(path, lineage_counts):
    """Write one Krona text file: count then one column per rank."""
    with open(path, 'w') as handle:
        for lineage, count in sorted(lineage_counts.items(), key=lambda kv: -kv[1]):
            if count <= 0:
                continue
            printable = int(count) if float(count).is_integer() else count
            handle.write('\t'.join([str(printable)] + list(lineage)) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--table', required=True, help='Feature table (.qza or .biom)')
    parser.add_argument('--taxonomy', required=True, help='Taxonomy (.qza or TSV)')
    parser.add_argument('--outdir', required=True, help='Directory for Krona text files')
    parser.add_argument('--manifest', required=True,
                        help='TSV of file path and Krona dataset label, in plot order')
    parser.add_argument('--per-sample', choices=['yes', 'no'], default='yes',
                        help='Also write one dataset per sample (default: yes)')
    args = parser.parse_args()

    table = load_table(args.table)
    taxonomy = load_taxonomy(args.taxonomy)

    os.makedirs(args.outdir, exist_ok=True)

    lineages = {feature: tuple(split_taxon(taxonomy.get(feature)))
                for feature in map(str, table.ids('observation'))}
    n_missing = sum(1 for feature in lineages if feature not in taxonomy)

    summed = defaultdict(float)
    datasets = []

    # One pass over the table: accumulate per-sample and summed counts together.
    per_sample = defaultdict(lambda: defaultdict(float))
    for values, feature, _ in table.iter(axis='observation', dense=True):
        lineage = lineages[str(feature)]
        for sample, value in zip(table.ids('sample'), values):
            if value:
                per_sample[sample][lineage] += float(value)
                summed[lineage] += float(value)

    # The summed dataset is written first so Krona opens on it.
    summed_path = os.path.join(args.outdir, '00_all_samples.txt')
    write_dataset(summed_path, summed)
    datasets.append((summed_path, 'all_samples'))

    if args.per_sample == 'yes':
        for sample in table.ids('sample'):
            label = UNSAFE_CHARS.sub('_', str(sample))
            path = os.path.join(args.outdir, f'{label}.txt')
            write_dataset(path, per_sample[sample])
            datasets.append((path, label))

    with open(args.manifest, 'w') as handle:
        for path, label in datasets:
            handle.write(f'{path}\t{label}\n')

    print(f'Wrote {len(datasets)} Krona dataset(s) to {args.outdir}; '
          f'{len(lineages)} features, {n_missing} missing from the taxonomy.')
    if not summed:
        sys.exit('ERROR: no counts found; the feature table and taxonomy may not match.')


if __name__ == '__main__':
    main()
