#!/usr/bin/env bash
# Run REVAMP's taxonomy assignment steps on ASVs from a Tourmaline repseqs run.
#
# REVAMP (https://github.com/McAllister-NOAA/REVAMP) assigns taxonomy by BLASTing
# ASVs against a local NCBI nt database and merging the best hits to their lowest
# common ancestor. `revamp.sh` also runs cutadapt, DADA2, tables and figures, and
# prompts for input, so this script calls only the taxonomy scripts in REVAMP's
# assets/ directory, in the same order and with the same arguments revamp.sh uses.
#
# Nothing inside the REVAMP clone is modified: scripts are read from --revamp-dir
# and all output is written under --workdir, which mirrors a REVAMP output folder
# (dada2/, blast_results/, ASV2Taxonomy/) because REVAMP's scripts use relative
# paths between those directories.
#
# TAXONKIT_DB is pointed straight at <blastdb>/taxdump for the duration of this
# script. revamp.sh instead deletes and re-copies the dmp files in the user's
# TAXONKIT_DB directory on every run; that is skipped here.
#
# Modes:
#   blast   BLASTn ASVs against nt, raising -max_target_seqs until all best hits fit
#   assign  reformat BLAST hits, look up lineages with taxonkit, assign taxonomy
#
# Run by Tourmaline during the taxonomy step, in the `revamp` conda environment.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_revamp_taxonomy.sh --mode blast|assign --revamp-dir DIR --workdir DIR
                              --blastdb DIR [options]

  --mode         blast: run BLASTn against nt; assign: taxonomy from BLAST results
  --revamp-dir   REVAMP clone (its assets/ scripts are called, never modified)
  --workdir      REVAMP-style working directory (dada2/, blast_results/, ASV2Taxonomy/)
  --blastdb      Database directory: nt volumes (blast mode) and taxdump/ (assign mode)
  --run-name     Basename for REVAMP outputs (assign mode; default: revamp)
  --blast-mode   allIN | allEnvOUT | mostEnvOUT (blast mode; default: mostEnvOUT)
  --query-cov    Percent of ASV length a hit must cover (assign mode; default: 90)
  --cutoffs      Percent ID cutoffs S,G,F,O,C,P (assign mode; default: 97,95,90,80,70,60)
  --threads      Threads for BLASTn (blast mode; default: 1)
EOF
}

mode=""
revamp_dir=""
workdir=""
blastdb=""
run_name="revamp"
blast_mode="mostEnvOUT"
query_cov="90"
cutoffs="97,95,90,80,70,60"
threads="1"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) mode="$2"; shift 2 ;;
        --revamp-dir) revamp_dir="$2"; shift 2 ;;
        --workdir) workdir="$2"; shift 2 ;;
        --blastdb) blastdb="$2"; shift 2 ;;
        --run-name) run_name="$2"; shift 2 ;;
        --blast-mode) blast_mode="$2"; shift 2 ;;
        --query-cov) query_cov="$2"; shift 2 ;;
        --cutoffs) cutoffs="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

for required in mode revamp_dir workdir blastdb; do
    if [[ -z "${!required}" ]]; then
        echo "ERROR: --${required//_/-} is required" >&2
        usage >&2
        exit 1
    fi
done

if [[ ! -d "$revamp_dir/assets" ]]; then
    echo "ERROR: no assets/ directory in revamp_dir: $revamp_dir" >&2
    echo "       revamp_dir must point at a clone of https://github.com/McAllister-NOAA/REVAMP" >&2
    exit 1
fi
if [[ ! -d "$blastdb" ]]; then
    echo "ERROR: revamp_blastdb does not exist: $blastdb" >&2
    exit 1
fi

revamp_dir="$(cd "$revamp_dir" && pwd)"
workdir="$(cd "$workdir" && pwd)"
blastdb="$(cd "$blastdb" && pwd)"

# A host install earlier in PATH (e.g. /usr/local/bin/Rscript) can shadow the conda
# environment's, and the host R usually lacks dplyr/Biostrings. Prefer the env's own.
prefer_env_bin() {
    local name="$1"
    if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/$name" ]]; then
        echo "$CONDA_PREFIX/bin/$name"
    else
        echo "$name"
    fi
}
RSCRIPT="$(prefer_env_bin Rscript)"
PERL="$(prefer_env_bin perl)"

# Largest per-ASV hit count in a blast_assessment.pl report. Computed with awk rather
# than `sort | head`, which exits 141 (SIGPIPE) under `set -o pipefail` once the file is
# long enough that sort is still writing when head closes the pipe.
max_hit_count() {
    awk -F'\t' 'NF>1 && $2+0>m {m=$2+0} END {print m+0}' "$1"
}

# BLAST looks for its taxonomy files (taxdb.*, taxonomy4blast.sqlite3) along BLASTDB, not
# only beside -db. Put the database directory on that path so taxid filtering can work.
export BLASTDB="${blastdb}${BLASTDB:+:$BLASTDB}"

case "$mode" in
blast)
    if [[ ! -f "$workdir/dada2/ASVs.fa" ]]; then
        echo "ERROR: missing $workdir/dada2/ASVs.fa" >&2
        exit 1
    fi
    mkdir -p "$workdir/blast_results"

    # REVAMP excludes environmental/uncultured subjects using taxid lists built by
    # its ncbi_db_cleanup.sh. Resolve the list before BLASTing so a missing file is
    # reported here rather than as a BLAST error.
    negative_list=""
    case "$blast_mode" in
        allIN) ;;
        mostEnvOUT) negative_list="$blastdb/taxdump/taxid_exclusion_list_leavesinUnclassified.txt" ;;
        allEnvOUT) negative_list="$blastdb/taxdump/taxid_exclusion_list_removesUnclassified.txt" ;;
        *) echo "ERROR: revamp_blast_mode must be allIN, allEnvOUT or mostEnvOUT (got: $blast_mode)" >&2; exit 1 ;;
    esac
    if [[ -n "$negative_list" && ! -f "$negative_list" ]]; then
        echo "ERROR: missing taxid exclusion list for blast mode '$blast_mode':" >&2
        echo "       $negative_list" >&2
        echo "       Prepare the database with REVAMP's ncbi_db_cleanup.sh, or use revamp_blast_mode: allIN." >&2
        exit 1
    fi

    # -negative_taxidlist additionally needs BLAST's own taxonomy files. Without them
    # BLAST reports "The -taxids command line option requires additional data files"
    # and the exclusion list is not applied, silently giving allIN results.
    if [[ -n "$negative_list" ]]; then
        IFS=':' read -ra blastdb_paths <<< "$BLASTDB"
        # taxdb.btd/.bti carry the names; taxonomy4blast.sqlite3 is what BLAST 2.12+
        # uses to expand a taxid list to its subtree. A missing sqlite3 file is not
        # always reported cleanly -- it can crash the search minutes in.
        for required in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
            found=FALSE
            for dir in "${blastdb_paths[@]}"; do
                if [[ -f "$dir/$required" ]]; then
                    found=TRUE
                    break
                fi
            done
            if [[ "$found" = FALSE ]]; then
                echo "ERROR: blast mode '$blast_mode' filters by taxid, which needs BLAST's" >&2
                echo "       taxonomy files; $required was not found on BLASTDB ($BLASTDB)." >&2
                echo "       Install them into the database directory:" >&2
                echo "         cd $blastdb && update_blastdb.pl taxdb && tar -xzf taxdb.tar.gz" >&2
                echo "       (or: wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz)" >&2
                echo "       All three must come from the same taxdb download." >&2
                echo "       If they cannot be installed, set revamp_blast_mode: allIN, which" >&2
                echo "       does no taxid filtering." >&2
                exit 1
            fi
        done
        echo "Taxid filtering: $(wc -l < "$negative_list") taxids excluded ($blast_mode)"
    fi
    echo "BLAST version: $(blastn -version 2>&1 | head -1)"

    n_asvs=$(grep -c ">" "$workdir/dada2/ASVs.fa")
    max_target_seqs=4000
    runthroughcount=3
    pass_blast_scrutiny=FALSE

    # Same escalation loop as revamp.sh: if any ASV's hit list was truncated at
    # max_target_seqs, the best-hit set may be incomplete, so BLAST again with a
    # higher ceiling (up to three extra attempts).
    while [[ "$pass_blast_scrutiny" = FALSE ]]; do
        runthroughcount=$((runthroughcount + 1))
        echo "Running BLASTn (max_target_seqs=${max_target_seqs}): $(date)"
        # shellcheck disable=SC2086
        blastn -db "$blastdb/nt" \
            -query "$workdir/dada2/ASVs.fa" \
            -outfmt '6 qseqid pident length staxids sacc' \
            -subject_besthit \
            -max_target_seqs "$max_target_seqs" \
            -num_threads "$threads" \
            -out "$workdir/blast_results/ASV_blastn_nt.btab" \
            ${negative_list:+-negative_taxidlist "$negative_list"}

        if [[ ! -s "$workdir/blast_results/ASV_blastn_nt.btab" ]]; then
            echo "ERROR: BLASTn produced no hits ($workdir/blast_results/ASV_blastn_nt.btab" >&2
            echo "       is empty). Check the messages above and that -db $blastdb/nt exists." >&2
            exit 1
        fi

        "$PERL" "$revamp_dir/assets/blast_assessment.pl" \
            -i "$workdir/blast_results/ASV_blastn_nt.btab" \
            -c "$n_asvs" > "$workdir/blast_results/checkmaxtargetseqs.txt"
        highest=$(max_hit_count "$workdir/blast_results/checkmaxtargetseqs.txt")

        if [[ "$highest" -lt "$max_target_seqs" ]]; then
            pass_blast_scrutiny=TRUE
        elif [[ "$runthroughcount" -le 3 ]]; then
            max_target_seqs=$((highest * (runthroughcount + 2)))
            echo "Rerun BLAST number ${runthroughcount}: new max_target_seqs = ${max_target_seqs}"
        else
            echo "WARNING: max_target_seqs inclusive of all top hits could not be reached"
            echo "         (used max_target_seqs = ${max_target_seqs}); some hit lists are truncated."
            pass_blast_scrutiny=TRUE
        fi
    done
    echo "Finished BLASTn: $(date)"
    ;;

assign)
    btab="$workdir/blast_results/ASV_blastn_nt.btab"
    if [[ ! -f "$btab" ]]; then
        echo "ERROR: missing $btab" >&2
        exit 1
    fi
    common_names="$blastdb/taxdump/common_names.dmp"
    if [[ ! -f "$common_names" ]]; then
        echo "ERROR: missing $common_names" >&2
        echo "       REVAMP's taxonomy script requires it. Create it from the NCBI taxonomy dump:" >&2
        echo "       grep \"genbank common name\" $blastdb/taxdump/names.dmp > $common_names" >&2
        exit 1
    fi
    for dmp in names.dmp nodes.dmp merged.dmp delnodes.dmp; do
        if [[ ! -f "$blastdb/taxdump/$dmp" ]]; then
            echo "ERROR: missing $blastdb/taxdump/$dmp (needed by taxonkit)" >&2
            exit 1
        fi
    done
    # `reformat2` was added in taxonkit 0.20.0; older versions fail obscurely here.
    if ! taxonkit reformat2 --help >/dev/null 2>&1; then
        echo "ERROR: this taxonkit has no 'reformat2' command; REVAMP needs taxonkit >= 0.20.0." >&2
        echo "       Installed: $(taxonkit version 2>&1 | head -1)" >&2
        exit 1
    fi

    export TAXONKIT_DB="$blastdb/taxdump"

    # Report truncated hit lists. With externally supplied BLAST results we cannot
    # re-run BLAST with a higher ceiling, so this is a warning, not an error.
    n_asvs=$(grep -c ">" "$workdir/dada2/ASVs.fa")
    "$PERL" "$revamp_dir/assets/blast_assessment.pl" -i "$btab" -c "$n_asvs" \
        > "$workdir/blast_results/checkmaxtargetseqs.txt"
    highest=$(max_hit_count "$workdir/blast_results/checkmaxtargetseqs.txt")
    if [[ "$highest" -ge 4000 ]]; then
        echo "WARNING: an ASV has ${highest} BLAST hits at its best percent identity."
        echo "         If BLASTn was run with -max_target_seqs ${highest} or lower, hit lists are"
        echo "         truncated and assignments may be deeper than the evidence supports."
        echo "         Re-run BLASTn with a higher -max_target_seqs."
    fi

    echo "Reformatting BLAST output: $(date)"
    "$RSCRIPT" --vanilla "$revamp_dir/assets/reformat_blast.R" "$workdir/blast_results" "$query_cov"

    mkdir -p "$workdir/ASV2Taxonomy"
    cd "$workdir/ASV2Taxonomy"

    echo "Looking up lineages with taxonkit: $(date)"
    cut -f4 ../blast_results/ASV_blastn_nt_formatted.txt \
        | sed -E "s/ //g" | tr ',' '\n' | tr ';' '\n' | sort | uniq | grep -v "taxid" > taxids.txt
    taxonkit lineage taxids.txt | awk '$2!=""' > taxonkit_out.txt
    taxonkit reformat2 taxonkit_out.txt | cut -f1,3 > reformatted_taxonkit_out.txt

    # Same as revamp.sh under -y: fill gaps in the K/P/C/O/F/G/S string, then replace
    # characters that would break the downstream perl/R scripts.
    "$PERL" "$revamp_dir/assets/fillIn_taxonkit.pl" -i reformatted_taxonkit_out.txt \
        > reformatted_taxonkit_out.txt_temp
    mv reformatted_taxonkit_out.txt reformatted_taxonkit_out_ORIGINAL.txt
    sed -E 's/[^A-Za-z0-9;[:blank:]]/_/g' reformatted_taxonkit_out.txt_temp \
        > reformatted_taxonkit_out.txt
    rm reformatted_taxonkit_out.txt_temp

    echo "Assigning taxonomy: $(date)"
    "$PERL" "$revamp_dir/assets/asv_taxonomy_processing_figureOuts.pl" \
        -a ../dada2/ASVs_counts.tsv \
        -s ../blast_results/ASV_blastn_nt_formatted.txt \
        -t reformatted_taxonkit_out.txt \
        -f "$cutoffs" \
        -n "$run_name" \
        -c "$common_names"

    if [[ ! -f "${run_name}_asvTaxonomyTable.txt" ]]; then
        echo "ERROR: REVAMP did not write ${run_name}_asvTaxonomyTable.txt" >&2
        exit 1
    fi
    echo "Finished REVAMP taxonomy: $(date)"
    ;;

*)
    echo "ERROR: --mode must be blast or assign (got: $mode)" >&2
    exit 1
    ;;
esac
