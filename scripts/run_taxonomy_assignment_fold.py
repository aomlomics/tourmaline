#!/usr/bin/env python
"""Run one tax-credit fold assignment using Tourmaline taxonomy QIIME commands."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BT2_CONDA_ENV = "bt2-blca"
BT2_INDEX_MARKER = "bowtie2_index.1.bt2"


def run_cmd(cmd: str, check: bool = True) -> None:
    print(cmd, flush=True)
    subprocess.run(cmd, shell=True, check=check)


def run_bt2_cmd(cmd: str, check: bool = True) -> None:
    run_cmd(f"conda run -n {BT2_CONDA_ENV} --no-capture-output {cmd}", check=check)


def _qiime_export(qza_path: Path, dest_path: Path, output_format: str) -> None:
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if dest_path.exists():
        dest_path.unlink()
    run_cmd(
        "qiime tools export "
        f"--input-path '{qza_path}' "
        f"--output-path '{dest_path}' "
        f"--output-format {output_format}"
    )


def _write_headerless_taxonomy_tsv(src_path: Path, dest_path: Path) -> None:
    with src_path.open(encoding="utf-8", errors="replace") as src, dest_path.open(
        "w", encoding="utf-8"
    ) as dest:
        for line in src:
            stripped = line.strip()
            if not stripped or stripped.startswith("Feature ID"):
                continue
            parts = stripped.split("\t")
            if len(parts) >= 2:
                dest.write(f"{parts[0]}\t{parts[1]}\n")


def ensure_fasta(path: Path, staging_dir: Path) -> Path:
    if path.suffix.lower() != ".qza":
        return path
    dest = staging_dir / f"{path.stem}.fasta"
    if not dest.exists() or dest.stat().st_mtime < path.stat().st_mtime:
        _qiime_export(path, dest, "DNAFASTAFormat")
    return dest


def ensure_taxonomy_tsv(path: Path, staging_dir: Path) -> Path:
    if path.suffix.lower() != ".qza":
        return path
    dest = staging_dir / f"{path.stem}.tsv"
    if dest.exists() and dest.stat().st_mtime >= path.stat().st_mtime:
        return dest
    try:
        _qiime_export(path, dest, "HeaderlessTSVTaxonomyFormat")
    except subprocess.CalledProcessError:
        tmp_tsv = dest.with_suffix(".qiime_export.tsv")
        if tmp_tsv.exists():
            tmp_tsv.unlink()
        _qiime_export(path, tmp_tsv, "TSVTaxonomyFormat")
        _write_headerless_taxonomy_tsv(tmp_tsv, dest)
        tmp_tsv.unlink()
    return dest


def build_bowtie_index(ref_fasta: Path, index_dir: Path, threads: int) -> Path:
    index_dir.mkdir(parents=True, exist_ok=True)
    prefix = index_dir / "bowtie2_index"
    marker = index_dir / BT2_INDEX_MARKER
    if marker.exists():
        return prefix
    run_bt2_cmd(
        "bowtie2-build "
        f"--threads {threads} "
        f"-f '{ref_fasta}' "
        f"'{prefix}'"
    )
    return prefix


def param_id_from_config(cfg: dict, confidence: float) -> str:
    method = cfg.get("classify_method", "naive-bayes")
    fit_params = (cfg.get("fit_params") or "").strip()
    if method == "naive-bayes" and fit_params:
        digest = hashlib.md5(fit_params.encode()).hexdigest()[:8]
        base = f"nb-{digest}"
    else:
        base = method.replace("-", "_")
    return f"{base}-conf{confidence}"


def assign_naive_bayes(
    query_qza: Path,
    ref_seqs: Path,
    ref_taxa: Path,
    out_dir: Path,
    cfg: dict,
    confidence: float,
    skip_fit: bool = False,
    classifier_qza: Path | None = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    classifier = classifier_qza or out_dir / "classifier.qza"
    taxonomy_qza = out_dir / "taxonomy.qza"
    taxonomy_tsv = out_dir / "taxonomy.tsv"
    assignments = out_dir / "query_tax_assignments.txt"
    fit_params = (cfg.get("fit_params") or "").strip()
    classify_params = (cfg.get("classify_params") or "").strip()
    threads = int(cfg.get("classify_threads", 1))

    if not skip_fit:
        run_cmd(
            "qiime feature-classifier fit-classifier-naive-bayes "
            f"--i-reference-reads {ref_seqs} "
            f"--i-reference-taxonomy {ref_taxa} "
            f"--o-classifier {classifier} "
            f"{fit_params}"
        )

    run_cmd(
        "qiime feature-classifier classify-sklearn "
        f"--i-classifier {classifier} "
        f"--i-reads {query_qza} "
        f"--p-confidence {confidence} "
        f"--o-classification {taxonomy_qza} "
        f"--p-n-jobs {threads} "
        f"{classify_params}"
    )
    run_cmd(
        "qiime tools export "
        f"--input-path {taxonomy_qza} "
        f"--output-path {taxonomy_tsv} "
        "--output-format TSVTaxonomyFormat"
    )
    shutil.copyfile(taxonomy_tsv, assignments)
    return assignments


def assign_consensus_blast(
    query_qza: Path,
    ref_seqs: Path,
    ref_taxa: Path,
    out_dir: Path,
    cfg: dict,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    taxonomy_qza = out_dir / "taxonomy.qza"
    taxonomy_tsv = out_dir / "taxonomy.tsv"
    assignments = out_dir / "query_tax_assignments.txt"
    classify_params = (cfg.get("classify_params") or "").strip()
    run_cmd(
        "qiime feature-classifier classify-consensus-blast "
        f"--i-reference-reads {ref_seqs} "
        f"--i-reference-taxonomy {ref_taxa} "
        f"--i-query {query_qza} "
        f"--p-perc-identity {cfg['perc_identity']} "
        f"--p-query-cov {cfg['query_cov']} "
        f"--p-min-consensus {cfg['min_consensus']} "
        f"--o-classification {taxonomy_qza} "
        f"--o-search-results {out_dir / 'search_results.qza'} "
        f"{classify_params}"
    )
    run_cmd(
        "qiime tools export "
        f"--input-path {taxonomy_qza} "
        f"--output-path {taxonomy_tsv} "
        "--output-format TSVTaxonomyFormat"
    )
    shutil.copyfile(taxonomy_tsv, assignments)
    return assignments


def assign_consensus_vsearch(
    query_qza: Path,
    ref_seqs: Path,
    ref_taxa: Path,
    out_dir: Path,
    cfg: dict,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    taxonomy_qza = out_dir / "taxonomy.qza"
    taxonomy_tsv = out_dir / "taxonomy.tsv"
    assignments = out_dir / "query_tax_assignments.txt"
    classify_params = (cfg.get("classify_params") or "").strip()
    threads = int(cfg.get("classify_threads", 1))
    run_cmd(
        "qiime feature-classifier classify-consensus-vsearch "
        f"--i-reference-reads {ref_seqs} "
        f"--i-reference-taxonomy {ref_taxa} "
        f"--i-query {query_qza} "
        f"--p-perc-identity {cfg['perc_identity']} "
        f"--p-query-cov {cfg['query_cov']} "
        f"--p-min-consensus {cfg['min_consensus']} "
        f"--o-classification {taxonomy_qza} "
        f"--o-search-results {out_dir / 'search_results.qza'} "
        f"--p-threads {threads} "
        f"{classify_params}"
    )
    run_cmd(
        "qiime tools export "
        f"--input-path {taxonomy_qza} "
        f"--output-path {taxonomy_tsv} "
        "--output-format TSVTaxonomyFormat"
    )
    shutil.copyfile(taxonomy_tsv, assignments)
    return assignments


def assign_bt2_blca(
    query_qza: Path,
    ref_seqs: Path,
    ref_taxa: Path,
    out_dir: Path,
    cfg: dict,
    confidence: float,
    skip_fit: bool = False,
    bowtie_index_dir: Path | None = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    staging_dir = out_dir / "staging"
    query_fasta = ensure_fasta(query_qza, staging_dir)
    ref_fasta = ensure_fasta(ref_seqs, staging_dir)
    ref_taxonomy = ensure_taxonomy_tsv(ref_taxa, staging_dir)

    index_dir = bowtie_index_dir or out_dir / "bowtie2_index"
    if skip_fit and bowtie_index_dir is None:
        raise ValueError("bt2-blca --skip-fit requires --bowtie-index-dir")
    if not skip_fit:
        build_bowtie_index(ref_fasta, index_dir, int(cfg.get("classify_threads", 1)))

    index_prefix = index_dir / "bowtie2_index"
    temp_dir = out_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    sam_path = out_dir / "bowtie2_all.sam"
    raw_taxonomy = out_dir / "raw-taxonomy.tsv"
    taxonomy_tsv = out_dir / "taxonomy.tsv"
    assignments = out_dir / "query_tax_assignments.txt"
    taxa_ranks = cfg.get("taxa_ranks", "kingdom,phylum,class,order,family,genus,species")
    threads = int(cfg.get("classify_threads", 1))

    run_bt2_cmd(
        "bowtie2 "
        f"-x '{index_prefix}' -f -U '{query_fasta}' "
        f"-S '{temp_dir}/end_to_end.sam' "
        f"--no-hd --no-sq --very-sensitive --end-to-end --no-unal "
        f"-p {threads} -k 100 --un '{temp_dir}/end_to_end_reject.fasta'"
    )
    run_bt2_cmd(
        "bowtie2 "
        f"-x '{index_prefix}' -f -U '{temp_dir}/end_to_end_reject.fasta' "
        f"-S '{temp_dir}/local.sam' "
        f"--no-hd --no-sq --very-sensitive --local --no-unal "
        f"-p {threads} -k 100 --un '{temp_dir}/end_to_end_and_local_reject.fasta'"
    )
    with sam_path.open("w", encoding="utf-8") as out_sam:
        for part in ("end_to_end.sam", "local.sam"):
            part_path = temp_dir / part
            if part_path.exists():
                out_sam.write(part_path.read_text(encoding="utf-8", errors="replace"))

    blca_script = SCRIPT_DIR / "blca_from_bowtie.py"
    reformat_script = SCRIPT_DIR / "reformat_summary_for_r.py"
    run_cmd(
        "python "
        f"'{blca_script}' "
        f"-i '{sam_path}' -r '{ref_taxonomy}' -q '{ref_fasta}' "
        f"-b {cfg['perc_identity']} -l {cfg['query_cov']} "
        f"-p muscle -n 100 -m 1.0 -f 2.5 -g -2 "
        f"-tr '{taxa_ranks}' -o '{raw_taxonomy}'"
    )
    run_cmd(
        "python "
        f"'{reformat_script}' '{raw_taxonomy}' '{taxonomy_tsv}' "
        f"{confidence} '{taxa_ranks}'"
    )

    reject_fasta = temp_dir / "end_to_end_and_local_reject.fasta"
    if reject_fasta.exists() and reject_fasta.stat().st_size > 0:
        with taxonomy_tsv.open("a", encoding="utf-8") as out_tsv, reject_fasta.open(
            encoding="utf-8", errors="replace"
        ) as reject:
            for line in reject:
                if line.startswith(">"):
                    feature_id = line[1:].strip()
                    out_tsv.write(f"{feature_id}\tUnassigned\t0\n")

    shutil.copyfile(taxonomy_tsv, assignments)
    return assignments


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-qza", required=True)
    parser.add_argument("--ref-seqs", required=True)
    parser.add_argument("--ref-taxa", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--classify-method", default="naive-bayes")
    parser.add_argument("--confidence", type=float, default=0.7)
    parser.add_argument("--fit-params", default="")
    parser.add_argument("--classify-params", default="")
    parser.add_argument("--classify-threads", type=int, default=5)
    parser.add_argument("--perc-identity", type=float, default=0.8)
    parser.add_argument("--query-cov", type=float, default=0.8)
    parser.add_argument("--min-consensus", type=float, default=0.51)
    parser.add_argument("--skip-fit", action="store_true")
    parser.add_argument("--fit-only", action="store_true")
    parser.add_argument("--classifier-qza", default=None)
    parser.add_argument("--bowtie-index-dir", default=None)
    parser.add_argument(
        "--taxa-ranks",
        default="kingdom,phylum,class,order,family,genus,species",
    )
    args = parser.parse_args()

    cfg = {
        "classify_method": args.classify_method,
        "fit_params": args.fit_params,
        "classify_params": args.classify_params,
        "classify_threads": args.classify_threads,
        "perc_identity": args.perc_identity,
        "query_cov": args.query_cov,
        "min_consensus": args.min_consensus,
        "taxa_ranks": args.taxa_ranks,
    }
    query = Path(args.query_qza)
    ref_seqs = Path(args.ref_seqs)
    ref_taxa = Path(args.ref_taxa)
    out_dir = Path(args.output_dir)
    classifier = Path(args.classifier_qza) if args.classifier_qza else None
    bowtie_index_dir = Path(args.bowtie_index_dir) if args.bowtie_index_dir else None
    method = args.classify_method

    if args.fit_only:
        if method == "naive-bayes":
            fit_params = (cfg.get("fit_params") or "").strip()
            classifier = classifier or out_dir / "classifier.qza"
            run_cmd(
                "qiime feature-classifier fit-classifier-naive-bayes "
                f"--i-reference-reads {ref_seqs} "
                f"--i-reference-taxonomy {ref_taxa} "
                f"--o-classifier {classifier} "
                f"{fit_params}"
            )
            return 0
        if method == "bt2-blca":
            ref_fasta = ensure_fasta(ref_seqs, out_dir / "staging")
            build_bowtie_index(
                ref_fasta,
                bowtie_index_dir or out_dir / "bowtie2_index",
                int(cfg.get("classify_threads", 1)),
            )
            return 0
        print(f"--fit-only supports naive-bayes and bt2-blca only, not {method}", file=sys.stderr)
        return 1

    if method == "naive-bayes":
        assign_naive_bayes(
            query, ref_seqs, ref_taxa, out_dir, cfg, args.confidence,
            skip_fit=args.skip_fit, classifier_qza=classifier,
        )
    elif method == "bt2-blca":
        assign_bt2_blca(
            query, ref_seqs, ref_taxa, out_dir, cfg, args.confidence,
            skip_fit=args.skip_fit, bowtie_index_dir=bowtie_index_dir,
        )
    elif method == "consensus-blast":
        assign_consensus_blast(query, ref_seqs, ref_taxa, out_dir, cfg)
    elif method == "consensus-vsearch":
        assign_consensus_vsearch(query, ref_seqs, ref_taxa, out_dir, cfg)
    else:
        print(f"Unsupported classify_method for tax-credit fold: {method}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
