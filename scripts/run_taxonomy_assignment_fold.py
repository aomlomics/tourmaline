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


def run_cmd(cmd: str, check: bool = True) -> None:
    print(cmd, flush=True)
    subprocess.run(cmd, shell=True, check=check)


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
    args = parser.parse_args()

    cfg = {
        "classify_method": args.classify_method,
        "fit_params": args.fit_params,
        "classify_params": args.classify_params,
        "classify_threads": args.classify_threads,
        "perc_identity": args.perc_identity,
        "query_cov": args.query_cov,
        "min_consensus": args.min_consensus,
    }
    query = Path(args.query_qza)
    ref_seqs = Path(args.ref_seqs)
    ref_taxa = Path(args.ref_taxa)
    out_dir = Path(args.output_dir)
    classifier = Path(args.classifier_qza) if args.classifier_qza else None
    method = args.classify_method

    if args.fit_only:
        if method != "naive-bayes":
            print("--fit-only supports naive-bayes only", file=sys.stderr)
            return 1
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

    if method == "naive-bayes":
        assign_naive_bayes(
            query, ref_seqs, ref_taxa, out_dir, cfg, args.confidence,
            skip_fit=args.skip_fit, classifier_qza=classifier,
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
