#!/usr/bin/env python
"""Build the small mock-community test fixture in 00-data/tax-credit-test/mock/.

Deterministic and standard-library only. ASV sequences are taken from the
mifishGom-test database, so these are smoke-test fixtures, not a realistic
benchmark (every ASV has an exact match in that database).

Writes:
  mock-rep-seqs.fasta           ASV sequences (md5 feature ids, like DADA2 in QIIME 2)
  mock-feature-table.tsv        read counts: mock-even, mock-staggered, blank-1
  mock-composition-gom.tsv      expected composition, mifishGom-test backbone
  mock-asv-taxonomy-addJ.tsv    known ASV taxonomy, addJ-test backbone

Run from the tourmaline repository root:
  python scripts/make_mock_test_fixture.py
"""

import csv
import hashlib
from pathlib import Path

FIXTURE_DIR = Path("00-data/tax-credit-test")
OUT_DIR = FIXTURE_DIR / "mock"

# Species in both test databases. Antigonia combatia and Scatophagus argus
# have different lineages in the two databases (the backbone gotcha).
SHARED_SPECIES = [
    "Antigonia combatia",
    "Scatophagus argus",
    "Rachycentron canadum",
    "Lophius americanus",
    "Lampris guttatus",
    "Acipenser oxyrinchus",
    "Hexanchus griseus",
    "Platax orbicularis",
    "Brama dussumieri",
    "Anoplogaster cornuta",
]
# Species only in mifishGom-test: 'not_in_database' for addJ-test.
GOM_ONLY_SPECIES = 2
STAGGERED = [0.30, 0.20, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04, 0.03, 0.01, 0.005, 0.005]
READS_PER_SAMPLE = 20000


def read_taxonomy(fp):
    with open(fp) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return {r[0]: r[1] for r in rows[1:]}


def read_fasta(fp):
    seqs, current = {}, None
    with open(fp) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    gom_tax = read_taxonomy(FIXTURE_DIR / "mifishGom-test-taxa.tsv")
    gom_seqs = read_fasta(FIXTURE_DIR / "mifishGom-test-seqs.fasta")
    addj_tax = read_taxonomy(FIXTURE_DIR / "addJ-test-taxa.tsv")
    addj_by_species = {t.split(";")[-1]: t for t in addj_tax.values()}

    gom_by_species = {}
    for seq_id in sorted(gom_seqs):
        species = gom_tax[seq_id].split(";")[-1]
        gom_by_species.setdefault(species, seq_id)
    gom_only = [s for s in sorted(gom_by_species)
                if s not in addj_by_species and ";NA" not in gom_tax[gom_by_species[s]]
                and s.count(" ") == 1][:GOM_ONLY_SPECIES]
    species = SHARED_SPECIES + gom_only

    # one contaminant ASV with no known taxonomy, from addJ-test
    addj_seqs = read_fasta(FIXTURE_DIR / "addJ-test-seqs.fasta")
    contaminant = sorted(addj_seqs)[0]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    asvs = []
    for name in species:
        seq = gom_seqs[gom_by_species[name]]
        asvs.append((hashlib.md5(seq.encode()).hexdigest(), seq, name))
    contaminant_seq = addj_seqs[contaminant]
    contaminant_id = hashlib.md5(contaminant_seq.encode()).hexdigest()

    with open(OUT_DIR / "mock-rep-seqs.fasta", "w") as fh:
        for asv_id, seq, _ in asvs:
            fh.write(f">{asv_id}\n{seq}\n")
        fh.write(f">{contaminant_id}\n{contaminant_seq}\n")

    even = 1 / len(species)
    # read counts deviate from the designed proportions, as PCR bias would
    bias = [1.0, 0.6, 1.3, 1.1, 0.9, 1.5, 0.7, 1.0, 1.2, 0.8, 1.0, 0.9]
    with open(OUT_DIR / "mock-feature-table.tsv", "w") as fh:
        fh.write("Feature ID\tmock-even\tmock-staggered\tblank-1\n")
        for (asv_id, _, _), b, stag in zip(asvs, bias, STAGGERED):
            fh.write(f"{asv_id}\t{round(READS_PER_SAMPLE * even * b)}\t"
                     f"{round(READS_PER_SAMPLE * stag * b)}\t0\n")
        fh.write(f"{contaminant_id}\t150\t400\t35\n")

    with open(OUT_DIR / "mock-composition-gom.tsv", "w") as fh:
        fh.write("Taxonomy\tmock-even\tmock-staggered\n")
        for name, stag in zip(species, STAGGERED):
            fh.write(f"{gom_tax[gom_by_species[name]]}\t{even:.6f}\t{stag}\n")

    with open(OUT_DIR / "mock-asv-taxonomy-addJ.tsv", "w") as fh:
        fh.write("Feature ID\tTaxon\n")
        for asv_id, _, name in asvs:
            lineage = addj_by_species.get(name, gom_tax[gom_by_species[name]])
            fh.write(f"{asv_id}\t{lineage}\n")

    print(f"wrote {len(asvs) + 1} ASVs to {OUT_DIR}")


if __name__ == "__main__":
    main()
