#!/usr/bin/env python
"""Build a small, structure-preserving test subsample of a reference database.

Produces drop-in replacements for a Tourmaline tax-credit ``refseqs_file`` /
``taxa_file`` pair that are small enough for the tax-credit step to finish in
minutes, while keeping enough taxonomic structure for cross-validated and
novel-taxa evaluation to remain meaningful.

Why not a plain random sample
-----------------------------
``tax_credit.framework_functions.generate_cross_validated_sequences`` does not
split randomly: it builds a taxonomy tree and calls ``get_strata(tree, k)``,
which keeps any lineage node with at least ``k`` (= ``iterations``) sequences as
its own stratum and aggregates smaller siblings into a "misc" pool before
running ``StratifiedKFold``.  A uniform random sample of a reference database is
almost all singleton species, so every stratum collapses toward the root and the
train/test split stops testing anything.  Novel-taxa evaluation has a second
requirement: holding out a genus/family/order only works when *sibling* taxa
remain in the reference set.

This script therefore samples in tiers:

1. **Multi-sequence tier** - species carrying >= ``--min-multi-seqs`` sequences,
   contributing up to ``--max-seqs-per-species`` each.  These become the real
   species-level strata.
2. **Sibling backfill** - for every retained family, guarantee >= 2 genera; for
   every retained order, guarantee >= 2 families, so novel-taxa levels 6/5/4
   have both query and reference material.
3. **Singleton tier** - one sequence per species, spread round-robin across the
   retained families, reproducing the long tail of a real database.

Determinism
-----------
All sampling goes through a single ``random.Random`` seeded with ``--seed``
(default ``20260910``).  Candidate lists are sorted before shuffling so results
never depend on set/dict iteration order.  Re-running with the same inputs,
seed and size arguments reproduces the output byte for byte.

Format fidelity
---------------
Sequence records are emitted with their original FASTA header and sequence
string verbatim (both source databases are unwrapped, one line per sequence,
with bare accession headers such as ``AB021901`` or
``OP537910.1_representative_of_4_identical_accessions``).  Taxonomy rows keep
the ``Feature ID``/``Taxon`` header and all seven semicolon-separated ranks,
including literal ``NA`` placeholders - ``clean_database``'s junk filter does
not treat ``NA`` as junk, so those lineages are deliberately preserved.  Output
records follow input file order.

Inputs may be QIIME 2 ``.qza`` artifacts or plain FASTA/TSV; ``.qza`` payloads
are read directly from the artifact.  With ``--write-qza`` (requires the
``qiime2-amplicon-2024.10`` environment) the subsample is also imported back to
``.qza`` so it is a literal drop-in for a config that points at artifacts.

Example
-------
    conda run -n qiime2-amplicon-2024.10 python scripts/subsample_reference_db.py \
        --seqs  /path/to/db-seqs.qza \
        --taxa  /path/to/db-taxa.qza \
        --n-target 400 \
        --out-prefix 00-data/tax-credit-test/mifishGom-test \
        --write-qza
"""

from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict, OrderedDict
from os.path import abspath, basename, dirname, exists, join, splitext

DEFAULT_SEED = 20260910
N_RANKS = 7
TAXONOMY_HEADER = ("Feature ID", "Taxon")


# --------------------------------------------------------------------------
# Reading (.qza payloads are read straight out of the artifact zip)
# --------------------------------------------------------------------------

def _read_qza_member(path: str, suffix: str) -> list:
    """Return the lines of the single ``data/`` member ending in ``suffix``."""
    with zipfile.ZipFile(path) as zf:
        members = [
            n for n in zf.namelist()
            if "/data/" in n and n.endswith(suffix) and not n.endswith("/")
        ]
        if len(members) != 1:
            raise SystemExit(
                f"{path}: expected exactly one data/*{suffix} member, "
                f"found {len(members)}"
            )
        with zf.open(members[0]) as fh:
            return fh.read().decode("utf-8").splitlines()


def _read_lines(path: str, qza_suffix: str) -> list:
    if path.lower().endswith(".qza"):
        return _read_qza_member(path, qza_suffix)
    with open(path) as fh:
        return fh.read().splitlines()


def read_sequences(path: str) -> "OrderedDict[str, str]":
    """Parse FASTA into ``{header_id: sequence}`` preserving file order."""
    lines = _read_lines(path, "dna-sequences.fasta")
    seqs: "OrderedDict[str, str]" = OrderedDict()
    current = None
    chunks: list = []
    for line in lines:
        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(chunks)
            current = line[1:].strip()
            chunks = []
        elif current is not None:
            chunks.append(line.strip())
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


def read_taxonomy(path: str) -> "OrderedDict[str, str]":
    """Parse a taxonomy TSV into ``{feature_id: taxon}`` preserving file order.

    Tolerates both the headed ``TSVTaxonomyFormat`` and the headerless variant.
    """
    lines = _read_lines(path, "taxonomy.tsv")
    taxa: "OrderedDict[str, str]" = OrderedDict()
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            raise SystemExit(f"{path}: line {i + 1} is not tab-separated")
        if i == 0 and parts[0].strip().lower() in ("feature id", "featureid"):
            continue
        taxa[parts[0].strip()] = parts[1].strip()
    return taxa


# --------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------

class Record:
    __slots__ = ("fid", "seq", "taxon", "ranks")

    def __init__(self, fid: str, seq: str, taxon: str):
        self.fid = fid
        self.seq = seq
        self.taxon = taxon
        self.ranks = [r.strip() for r in taxon.split(";")]

    def key(self, depth: int) -> str:
        return ";".join(self.ranks[:depth])


def build_records(seqs: dict, taxa: dict, min_len: int) -> list:
    """Pair sequences with taxonomy, dropping unpaired or short records."""
    records = []
    missing_taxa = 0
    too_short = 0
    wrong_ranks = 0
    for fid, seq in seqs.items():
        taxon = taxa.get(fid)
        if taxon is None:
            missing_taxa += 1
            continue
        if len(seq) < min_len:
            too_short += 1
            continue
        rec = Record(fid, seq, taxon)
        if len(rec.ranks) != N_RANKS:
            wrong_ranks += 1
            continue
        records.append(rec)
    if missing_taxa:
        print(f"  note: {missing_taxa} sequences had no taxonomy row (skipped)")
    if too_short:
        print(f"  note: {too_short} sequences shorter than {min_len} bp (skipped)")
    if wrong_ranks:
        print(f"  note: {wrong_ranks} lineages were not {N_RANKS} ranks (skipped)")
    return records


def subsample(records: list, n_target: int, rng: random.Random,
              n_multi_species: int, min_multi_seqs: int,
              max_seqs_per_species: int) -> list:
    """Select records in tiers; returns the chosen records in input order."""
    by_species = OrderedDict()
    for rec in records:
        by_species.setdefault(rec.key(N_RANKS), []).append(rec)

    # Family (5 ranks) and order (4 ranks) membership, for sibling backfill.
    genera_in_family = defaultdict(set)
    families_in_order = defaultdict(set)
    species_in_genus = defaultdict(list)
    for sp_key, recs in by_species.items():
        ranks = recs[0].ranks
        genera_in_family[";".join(ranks[:5])].add(";".join(ranks[:6]))
        families_in_order[";".join(ranks[:4])].add(";".join(ranks[:5]))
        species_in_genus[";".join(ranks[:6])].append(sp_key)

    chosen: "OrderedDict[str, list]" = OrderedDict()   # species key -> records
    budget = n_target

    def take(sp_key: str, n_seqs: int) -> int:
        """Add ``n_seqs`` sequences of ``sp_key``; returns how many were added."""
        nonlocal budget
        if budget <= 0 or sp_key in chosen:
            return 0
        pool = sorted(by_species[sp_key], key=lambda r: r.fid)
        n = min(n_seqs, len(pool), budget)
        if n <= 0:
            return 0
        chosen[sp_key] = rng.sample(pool, n) if n < len(pool) else list(pool)
        budget -= n
        return n

    # --- Tier 1: multi-sequence species, spread round-robin across families ---
    multi_by_family = defaultdict(list)
    for sp_key, recs in by_species.items():
        if len(recs) >= min_multi_seqs:
            multi_by_family[";".join(recs[0].ranks[:5])].append(sp_key)

    # Prefer families that have >= 2 genera: those keep novel-taxa level 6 honest.
    fams = sorted(multi_by_family)
    rng.shuffle(fams)
    fams.sort(key=lambda f: 0 if len(genera_in_family[f]) >= 2 else 1)
    for fam in fams:
        rng.shuffle(multi_by_family[fam])

    picked_multi = 0
    round_idx = 0
    while picked_multi < n_multi_species and budget > 0:
        progressed = False
        for fam in fams:
            if picked_multi >= n_multi_species or budget <= 0:
                break
            candidates = multi_by_family[fam]
            if round_idx >= len(candidates):
                continue
            sp_key = candidates[round_idx]
            n_avail = len(by_species[sp_key])
            n_take = min(max_seqs_per_species, max(min_multi_seqs, n_avail))
            if take(sp_key, n_take):
                picked_multi += 1
                progressed = True
        if not progressed:
            break
        round_idx += 1

    # --- Tier 2: sibling backfill so novel-taxa has something to hold out ---
    def singleton_from(sp_keys: list) -> bool:
        for sp_key in sp_keys:
            if sp_key not in chosen and take(sp_key, 1):
                return True
        return False

    # Every retained family should carry >= 2 genera.
    for fam in sorted({";".join(r.ranks[:5]) for recs in chosen.values() for r in recs}):
        present = {";".join(r.ranks[:6]) for recs in chosen.values()
                   for r in recs if ";".join(r.ranks[:5]) == fam}
        if len(present) >= 2:
            continue
        others = sorted(genera_in_family[fam] - present)
        rng.shuffle(others)
        for gen in others:
            if singleton_from(sorted(species_in_genus[gen])):
                break

    # Every retained order should carry >= 2 families.
    for order in sorted({";".join(r.ranks[:4]) for recs in chosen.values() for r in recs}):
        present = {";".join(r.ranks[:5]) for recs in chosen.values()
                   for r in recs if ";".join(r.ranks[:4]) == order}
        if len(present) >= 2:
            continue
        others = sorted(families_in_order[order] - present)
        rng.shuffle(others)
        for fam in others:
            gens = sorted(genera_in_family[fam])
            rng.shuffle(gens)
            if any(singleton_from(sorted(species_in_genus[g])) for g in gens):
                break

    # --- Tier 3: singleton tail, round-robin across families for breadth ---
    singles_by_family = defaultdict(list)
    for sp_key, recs in by_species.items():
        if sp_key not in chosen:
            singles_by_family[";".join(recs[0].ranks[:5])].append(sp_key)
    single_fams = sorted(singles_by_family)
    rng.shuffle(single_fams)
    # Families already represented come first, so the tail thickens real lineages
    # before reaching for unrelated ones.
    retained_fams = {";".join(r.ranks[:5]) for recs in chosen.values() for r in recs}
    single_fams.sort(key=lambda f: 0 if f in retained_fams else 1)
    for fam in single_fams:
        rng.shuffle(singles_by_family[fam])

    round_idx = 0
    while budget > 0:
        progressed = False
        for fam in single_fams:
            if budget <= 0:
                break
            candidates = singles_by_family[fam]
            if round_idx >= len(candidates):
                continue
            if take(candidates[round_idx], 1):
                progressed = True
        if not progressed:
            break
        round_idx += 1

    # --- Keep at least one NA-bearing lineage, which clean_database preserves ---
    has_na = any("NA" in r.ranks for recs in chosen.values() for r in recs)
    if not has_na:
        na_species = sorted(k for k, recs in by_species.items()
                            if "NA" in recs[0].ranks and k not in chosen)
        if na_species:
            rng.shuffle(na_species)
            budget += 1          # allow a single overshoot rather than drop it
            take(na_species[0], 1)

    selected_ids = {r.fid for recs in chosen.values() for r in recs}
    return [r for r in records if r.fid in selected_ids]


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------

def write_fasta(records: list, path: str) -> None:
    """Write unwrapped FASTA, matching the source databases' layout."""
    with open(path, "w") as fh:
        for rec in records:
            fh.write(f">{rec.fid}\n{rec.seq}\n")


def write_taxonomy(records: list, path: str) -> None:
    """Write a headed TSVTaxonomyFormat file, ranks preserved verbatim."""
    with open(path, "w") as fh:
        fh.write("\t".join(TAXONOMY_HEADER) + "\n")
        for rec in records:
            fh.write(f"{rec.fid}\t{rec.taxon}\n")


def qiime_import(input_path: str, out_path: str, semantic_type: str,
                 input_format: str = None) -> None:
    cmd = [
        "qiime", "tools", "import",
        "--type", semantic_type,
        "--input-path", input_path,
        "--output-path", out_path,
    ]
    if input_format:
        cmd += ["--input-format", input_format]
    print("  " + " ".join(cmd))
    if exists(out_path):
        os.remove(out_path)
    subprocess.run(cmd, check=True)


def summarize(records: list, label: str) -> None:
    species = Counter(r.key(N_RANKS) for r in records)
    genera = {r.key(6) for r in records}
    families = {r.key(5) for r in records}
    orders = {r.key(4) for r in records}
    multi = sum(1 for v in species.values() if v >= 3)
    fam_2gen = sum(
        1 for f in families
        if len({r.key(6) for r in records if r.key(5) == f}) >= 2
    )
    ord_2fam = sum(
        1 for o in orders
        if len({r.key(5) for r in records if r.key(4) == o}) >= 2
    )
    lengths = sorted(len(r.seq) for r in records)
    print(f"  {label}: {len(records)} seqs | {len(species)} species "
          f"({multi} with >=3 seqs) | {len(genera)} genera | "
          f"{len(families)} families ({fam_2gen} with >=2 genera) | "
          f"{len(orders)} orders ({ord_2fam} with >=2 families)")
    print(f"  length min/median/max: {lengths[0]}/{lengths[len(lengths) // 2]}"
          f"/{lengths[-1]} | NA-bearing lineages: "
          f"{sum(1 for r in records if 'NA' in r.ranks)}")


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--seqs", required=True,
                   help="Reference sequences (.qza or FASTA)")
    p.add_argument("--taxa", required=True,
                   help="Reference taxonomy (.qza or TSV)")
    p.add_argument("--out-prefix", required=True,
                   help="Output path prefix, e.g. 00-data/tax-credit-test/db-test")
    p.add_argument("--n-target", type=int, default=400,
                   help="Approximate number of sequences to keep (default: 400)")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED,
                   help=f"Random seed (default: {DEFAULT_SEED})")
    p.add_argument("--n-multi-species", type=int, default=50,
                   help="Species contributing multiple sequences (default: 50)")
    p.add_argument("--min-multi-seqs", type=int, default=3,
                   help="Minimum sequences for the multi-sequence tier "
                        "(default: 3; keep >= the tax-credit `iterations` value)")
    p.add_argument("--max-seqs-per-species", type=int, default=5,
                   help="Cap on sequences kept per species (default: 5)")
    p.add_argument("--min-length", type=int, default=100,
                   help="Drop sequences shorter than this (default: 100, "
                        "matching the config's min_read_length)")
    p.add_argument("--write-qza", action="store_true",
                   help="Also import outputs to .qza (needs the QIIME 2 env)")
    args = p.parse_args()

    seqs_path = abspath(args.seqs)
    taxa_path = abspath(args.taxa)
    for path in (seqs_path, taxa_path):
        if not exists(path):
            raise SystemExit(f"input not found: {path}")

    out_prefix = abspath(args.out_prefix)
    os.makedirs(dirname(out_prefix), exist_ok=True)

    print(f"reading {basename(seqs_path)} / {basename(taxa_path)}")
    seqs = read_sequences(seqs_path)
    taxa = read_taxonomy(taxa_path)
    records = build_records(seqs, taxa, args.min_length)
    summarize(records, "source")

    rng = random.Random(args.seed)
    selected = subsample(
        records,
        n_target=args.n_target,
        rng=rng,
        n_multi_species=args.n_multi_species,
        min_multi_seqs=args.min_multi_seqs,
        max_seqs_per_species=args.max_seqs_per_species,
    )
    summarize(selected, f"subsample (seed={args.seed})")

    fasta_out = out_prefix + "-seqs.fasta"
    tsv_out = out_prefix + "-taxa.tsv"
    write_fasta(selected, fasta_out)
    write_taxonomy(selected, tsv_out)
    print(f"  wrote {fasta_out}")
    print(f"  wrote {tsv_out}")

    if args.write_qza:
        qiime_import(fasta_out, out_prefix + "-seqs.qza",
                     "FeatureData[Sequence]")
        qiime_import(tsv_out, out_prefix + "-taxa.qza",
                     "FeatureData[Taxonomy]", "TSVTaxonomyFormat")
        print(f"  wrote {out_prefix}-seqs.qza")
        print(f"  wrote {out_prefix}-taxa.qza")

    return 0


if __name__ == "__main__":
    sys.exit(main())
