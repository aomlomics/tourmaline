# Adapted from https://github.com/qunfengdong/BLCA/ by Jesse Gomer
# for the University of California Conservation Consortium's CALeDNA Program

# Example usage:
# python blca_from_bowtie.py -i take_3_local.sam -r CO1_labeled_taxonomy.txt -q CO1_.fasta -b 0.8 -p path/to/muscle
#
# Parallelism (--threads N): each query's muscle alignment, scoring and
# bootstrap run in a pool of N worker processes. The parent reads the SAM,
# reference and taxonomy once and sends each worker a small, self-contained
# task (the query, its hits' sequences and taxonomies), so workers never need a
# copy of the reference and any multiprocessing start method works. Output is
# written in the same order as a serial run.
#
# Reproducibility: each query's bootstrap draws come from its own RNG, seeded
# from --seed and the query ID, so results do not depend on --threads or on
# scheduling. (The original seeded nothing, so runs differed between calls.)

import sys
import os

if sys.version_info < (3, 0):
    from StringIO import StringIO
else:
    from io import StringIO

try:
    from Bio import AlignIO, SeqIO
except ImportError:
    sys.stderr.write("Error! BioPython is not detected!\n")
    sys.exit(1)

import random
import subprocess
import re
import zlib
import multiprocessing
from collections import deque, namedtuple, defaultdict
import argparse

'''
BLCA Core annotation tool
'''


class SamEntry(object):
    def __init__(self, raw_row):
        self.qname = raw_row[0]
        self.flag = raw_row[1]
        self.rname = raw_row[2]
        self.pos = raw_row[3]
        self.mapq = raw_row[4]
        self.cigar = raw_row[5]
        self.rnext = raw_row[6]
        self.pnext = raw_row[7]
        self.tlen = raw_row[8]
        self.seq = raw_row[9]
        self.qual = raw_row[10]
        self.alignment_scores = [int(score.split(':')[-1]) for score in raw_row[11:13]]
        # sometimes it is a column early  because there is no second best match
        if raw_row[17][:2] == 'MD':
            self.md_z_flags = raw_row[17].split(':')[-1]
        else:
            self.md_z_flags = raw_row[18].split(':')[-1]

        total_match_count, total_count = self.calulate_match_count(self.md_z_flags)
        self.identity_ratio_old = total_match_count / total_count
        self.identity_ratio = total_match_count / float(len(self.seq))

        self.soft_clipped_ratio = self.cigar_total_soft_clipping() / total_count

        self.total_match = total_match_count

    def calulate_match_count(self, md_z_flags):
        tokenizer = re.compile(r'(\d+)|(\^[A-Z])|([A-Z])')
        total_match_count = 0.0
        total_mismatch_count = 0.0
        for item in tokenizer.finditer(md_z_flags):
            match_count, deletion, mismatch = item.groups()
            if match_count:
                total_match_count += int(match_count)
            if deletion:
                total_mismatch_count += 1
            if mismatch:
                total_mismatch_count += 1

        return total_match_count, (total_match_count + total_mismatch_count)

    # I am not completely sure if this is the logic that you want to use for checking unmapped at the ends
    # This function returns the max S at either end of the CIGAR score, 0  if there is no S at both ends
    def cigar_max_s(self):
        beginning_s, ending_s = self.get_soft_clipping()
        return max(beginning_s, ending_s)

    def cigar_total_soft_clipping(self):
        beginning_s, ending_s = self.get_soft_clipping()
        return beginning_s + ending_s

    def get_soft_clipping(self):
        tokenizer = re.compile(r'(?:\d+)|(?:[A-Z=])')
        cigar_elements = tokenizer.findall(self.cigar)

        if cigar_elements[1] == 'S':
            beginning_s = int(cigar_elements[0])
        else:
            beginning_s = 0

        if cigar_elements[-1] == 'S':
            ending_s = int(cigar_elements[-2])
        else:
            ending_s = 0

        return (beginning_s, ending_s)


class NotAvailableHandler(object):
    def __init__(self):
        self.count = 0
        self.na_forms = {'na', 'NA', 'Not Available', 'not available', 'Not available'}
        self.encoded_na_form = 'NA;;'

    def encode_if_na(self, taxon):
        if taxon not in self.na_forms:
            return taxon
        self.count += 1
        return 'NA;;{}'.format(self.count)

    def decode_if_na(self, taxon):
        if not taxon.startswith(self.encoded_na_form):
            return taxon

        return 'NA'


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description='Bayesian-based LCA taxonomic classification method')
    ##### Required arguments #####
    required = parser.add_argument_group('required arguments')
    required.add_argument("-i", "--sam", help="Input SAM file", type=str, required=True)
    required.add_argument('-q', '--reference', help="Reference fasta file", type=str, required=True)
    required.add_argument("-r", "--tax", help="reference taxonomy file for the Database", type=str, required=True)

    ##### Taxonomy arguments #####
    taxoptions = parser.add_argument_group('taxonomy profiling options [filtering of hits]')
    taxoptions.add_argument("-tr", "--ranks",type=str, help="Taxa ranks in database.")
    taxoptions.add_argument("-n", "--nper", help="number of times to bootstrap. Default: 100", type=int, default=100)
    taxoptions.add_argument("-b", "--iset", help="minimum identity score to include", type=float, default=0.8)
    taxoptions.add_argument('-l', '--length', help="minimum length of hit to include relative to query", type=float,
                            default=0.5)
    taxoptions.add_argument('-s', '--softclipping', help='maximum soft clipped ratio to include', type=float, default=0.2)
    ##### Alignment control arguments #####
    alignoptions = parser.add_argument_group('alignment control arguments')
    alignoptions.add_argument("-m", "--match", default=1.0, help="alignment match score. Default: 1", type=float)
    alignoptions.add_argument("-f", "--mismatch", default=-2.5, help="alignment mismatch penalty. Default: -2.5",
                              type=float)
    alignoptions.add_argument("-g", "--ngap", default=-2.0, help="alignment gap penalty. Default: -2", type=float)
    ##### Other arguments #####
    optional = parser.add_argument_group('other arguments')
    optional.add_argument('-p', '--muscle', help='Path to call muscle default: muscle', default='muscle')

    optional.add_argument("-o","--outfile",help="output file name. Default: <fasta>.blca.out",type=str)
    optional.add_argument('--continue_mode', help="accepted for compatibility; has no effect (output is always rewritten)", action='store_true')
    optional.add_argument("--muscle_use_diags", help="pass the diag argument to muscle", action='store_true')
    optional.add_argument("--muscle_max_iterations", help="set the max number of iterations for muscle", type=int, default=16)
    optional.add_argument("-t", "--threads", type=int, default=1,
                          help="worker processes for the per-query muscle + bootstrap work. Default: 1")
    optional.add_argument("--seed", type=int, default=0,
                          help="base seed for the per-query bootstrap RNGs. Default: 0")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    return args


def check_program(prgname):
    '''Check whether a program has been installed and put in the PATH'''
    path = os.popen("which " + prgname).read().rstrip()
    if len(path) > 0 and os.path.exists(path):
        print(prgname + " is located in your PATH!")
    else:
        print("ERROR: " + prgname + " is NOT in your PATH, please set up " + prgname + "!")
        sys.exit(1)


def get_dic_from_aln(aln):
    '''Read in alignment and convert it into a dictionary'''
    alignment = AlignIO.read(aln, "clustal")
    alndic = {}
    for r in alignment:
        alndic[r.id] = list(r.seq)
    return alndic


def build_muscle_fasta(query_name, query_seq, hit_names, reference_sequences):
    """Build MUSCLE input with short IDs; CLUSTAL output truncates long names.

    Aliases (``r0``..``rN`` for hits, ``q`` for the query) are what make the
    query distinguishable from a reference hit that carries the *same* name --
    which is the normal case for self-validation, where the query set is the
    reference set. Scoring therefore stays in alias space and only resolves
    names for the taxonomy lookup; mapping back to names any earlier would
    collapse the query and its own reference record into one entry.
    """
    alias_to_name = {}
    lines = []
    for i, hit in enumerate(hit_names):
        alias = f"r{i}"
        alias_to_name[alias] = hit
        lines.append(f">{alias}\n{reference_sequences[hit]}\n")
    query_alias = "q"
    alias_to_name[query_alias] = query_name
    lines.append(f">{query_alias}\n{query_seq}\n")
    return "".join(lines), alias_to_name, query_alias


def pairwise_score(alndic, query, match, mismatch, ngap):
    '''Calculate pairwise alignment score given a query'''
    nt = ["A", "C", "T", "G", "g", "a", "c", "t"]
    hitscore = {}
    for k, v in alndic.items():
        if k != query:
            hitscore[k] = 0
            for i in range(len(v)):
                if (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(match)
                elif (alndic[query][i] not in nt) and (v[i] not in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(0)
                elif ((alndic[query][i] not in nt) or (v[i] not in nt)) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(mismatch)
                elif (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(ngap)
    total = float(sum(hitscore.values()))
    if total <= 0:
        total = 1
    for k, v in hitscore.items():
        hitscore[k] = v / total
    return hitscore


def random_aln_score(alndic, query, match, mismatch, ngap, rng=random):
    '''Randomize the alignment, and calculate the score'''
    nt = ["A", "C", "T", "G", "g", "a", "c", "t"]
    idx = []
    for i in range(len(list(alndic.values())[0])):
        idx.append(rng.choice(range(len(list(alndic.values())[0]))))

    hitscore = {}
    for k, v in alndic.items():
        if k != query:
            hitscore[k] = 0
            for i in idx:
                if (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(match)
                elif (alndic[query][i] not in nt) and (v[i] not in nt) and (alndic[query][i] == v[i]):
                    hitscore[k] += float(0)
                elif ((alndic[query][i] not in nt) or (v[i] not in nt)) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(mismatch)
                elif (alndic[query][i] in nt) and (v[i] in nt) and (alndic[query][i] != v[i]):
                    hitscore[k] += float(ngap)
    return hitscore

def get_gap_pos(query, alndic):
    '''Get the gap position in the alignment'''
    for i in range(len(alndic[query])):
        if alndic[query][i] != "-":
            start = i
            break
    for i in range(len(alndic[query]) - 1, 0, -1):
        if alndic[query][i] != "-":
            end = i
            break
    return start, end


def cut_gap(alndic, start, end):
    '''Given a start and end gap position, truncate the alignmnet'''
    trunc_alndic = {}
    for k_truc, v_truc in alndic.items():
        trunc_alndic[k_truc] = v_truc[start:end]
    return trunc_alndic


def read_tax_acc(taxfile, not_available_handler, levels):
    tx = open(taxfile)
    acctax = {}
    for l in tx:
        lne = l.rstrip().strip(";").split("\t")
        if len(lne) != 2:
            continue
        if (levels[0] + ':') not in l:
            taxons = [not_available_handler.encode_if_na(taxon) for taxon in lne[1].split(';')]
            acctax[lne[0].split('.')[0]] = dict(zip(levels, taxons))
        else:
            pairs = [x.split(":", 1) for x in lne[1].split(";")]
            encoded = [(level, not_available_handler.encode_if_na(taxon)) for level, taxon in pairs]
            acctax[lne[0].split(".")[0]] = dict(encoded)
    tx.close()
    return acctax


def query_seed(base_seed, query_name):
    '''Seed for one query's bootstrap RNG: stable across runs and processes.

    Derived from the query ID with crc32 rather than hash(), which Python
    randomises per process for strings.
    '''
    return (base_seed * 1000003) ^ zlib.crc32(query_name.encode("utf-8"))


# Per-worker settings, set once by _init_worker (and directly for --threads 1),
# so each task only carries what differs between queries.
_SETTINGS = {}


def _init_worker(settings):
    _SETTINGS.clear()
    _SETTINGS.update(settings)


def classify_query(task):
    '''Align, score and bootstrap one query; return its output line.

    *task* is ``(seqn, fifsa, alias_to_name, query_alias, hit_taxonomies,
    hits)``: the MUSCLE input already built from the query and its valid hits,
    and the taxonomy of every hit alias (``None`` where the accession has no
    taxonomy). Everything reference-sized stays in the parent; NA ranks were
    already encoded there, in query order, exactly as a serial run does.
    '''
    seqn, fifsa, alias_to_name, query_alias, hit_taxonomies, hits = task
    s = _SETTINGS
    levels = s["levels"]
    decode = NotAvailableHandler().decode_if_na

    muscle_call = [s["muscle"], '-quiet', '-clw', '-maxiters', str(s["maxiters"])]
    if s["use_diags"]:
        muscle_call.append('-diags')
    proc = subprocess.Popen(muscle_call, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    outs, errs = proc.communicate(fifsa.encode('utf-8'))

    # Keep MUSCLE's aliases: scoring must be able to tell the query apart from a
    # reference hit with an identical name (see build_muscle_fasta).
    alndic = get_dic_from_aln(StringIO(outs.decode('utf-8')))
    if query_alias not in alndic:
        return (f"[WARNING] Query {seqn} missing from MUSCLE alignment; skipping.",
                seqn + "\tUnclassified\n")
    ### get gap position and truncate the alignment###
    start, end = get_gap_pos(query_alias, alndic)
    trunc_alndic = cut_gap(alndic, start, end)
    match, mismatch, ngap = s["match"], s["mismatch"], s["ngap"]
    orgscore = pairwise_score(trunc_alndic, query_alias, match, mismatch, ngap)
    if not orgscore:
        # Only reachable if MUSCLE dropped every reference row; valid_hits is
        # non-empty by this point, so the alignment normally has >= 2 rows.
        return (f"[WARNING] No scoreable hits for {seqn} after alignment; skipping.",
                seqn + "\tUnclassified\n")
    ### start bootstrap ###
    rng = random.Random(query_seed(s["seed"], seqn))
    pervote = {}  # record vote after nper bootstrap
    for j in range(s["nper"]):
        random_scores = random_aln_score(trunc_alndic, query_alias, match, mismatch, ngap, rng)
        max_score = max(random_scores.values())
        hits_with_max_score = [k3 for k3, v3 in random_scores.items() if v3 == max_score]
        vote_share = 1.0 / len(hits_with_max_score)
        for hit in hits_with_max_score:
            if hit in pervote:
                pervote[hit] += vote_share
            else:
                pervote[hit] = vote_share

    ### normalize vote by total votes ###
    ttlvote = sum(pervote.values())
    for k4, v4 in pervote.items():
        pervote[k4] = v4 / ttlvote

    votes_by_level = {}
    for level in levels:
        votes_by_level[level] = defaultdict(int)

    for hit_alias in orgscore.keys():
        # orgscore/pervote are keyed by MUSCLE alias; the taxonomy was resolved
        # from the real accession in the parent.
        hit_taxonomy = hit_taxonomies.get(hit_alias)
        if hit_taxonomy is None:
            continue
        for level in levels:
            votes_by_level[level][hit_taxonomy[level]] += pervote.get(hit_alias, 0)

    out = [seqn + "\t"]
    for level in levels:
        levels_votes = votes_by_level[level]
        out.append(level + ":" + decode(max(levels_votes, key=levels_votes.get)) + ";")
    out.append("\t")
    for level in levels:
        levels_votes = votes_by_level[level]
        out.append(level + ":" + str(max(levels_votes.values())) + ";")
    out.append("\t" + ";".join(hits))
    out.append("\n")
    return None, "".join(out)


def build_tasks(input_sequences, reference_sequences, acc2tax, na_handler, levels, counters):
    '''Yield one task per query (or a finished line), in input order.

    Everything that reads or mutates the large shared structures happens here,
    in the parent and in query order, so the encoded NA ranks a query sees are
    exactly those of a serial run.
    '''
    for seqn, info in input_sequences.items():
        counters["count"] += 1

        # A query whose ID is also a reference ID is expected for self-validation
        # (the query set IS the reference set), so it is counted and reported once
        # at the end rather than skipped. The query is still distinguishable from
        # its own reference record during scoring because MUSCLE aliases are kept.
        # NB: acc2tax keys are truncated at the first dot by read_tax_acc, so the
        # query ID must be truncated the same way or this never matches an
        # accession like "MT901437.1".
        if seqn.split(".")[0] in acc2tax:
            counters["self_id_queries"] += 1

        ### Get all the hits list belong to the same query ###
        valid_hits = []
        seen_hits = set()
        for hit in info.hits:
            if hit not in reference_sequences:
                print("Missing reference sequence for " + hit)
                continue
            if hit in seen_hits:
                continue
            seen_hits.add(hit)
            valid_hits.append(hit)
        if not valid_hits:
            yield ("line", seqn + "\tUnclassified\n")
            continue

        fifsa, alias_to_name, query_alias = build_muscle_fasta(
            seqn, info.seq, valid_hits, reference_sequences
        )
        # Resolve each hit's taxonomy now, filling missing ranks with encoded NAs
        # (mutating acc2tax, as the serial loop always did, so a hit reused by a
        # later query keeps the same encoded NA).
        hit_taxonomies = {}
        for hit_alias, hit_name in alias_to_name.items():
            if hit_alias == query_alias:
                continue
            short_hit_name = hit_name.split(".")[0]
            if short_hit_name not in acc2tax:
                print("Missing taxonomy info for ", short_hit_name)
                hit_taxonomies[hit_alias] = None
                continue
            hit_taxonomy = acc2tax[short_hit_name]
            for level in levels:
                # deal with missing values in the taxonomy
                if level not in hit_taxonomy:
                    hit_taxonomy[level] = na_handler.encode_if_na("NA")
            hit_taxonomies[hit_alias] = dict(hit_taxonomy)
        yield ("task", (seqn, fifsa, alias_to_name, query_alias, hit_taxonomies, info.hits))


def main(argv=None):
    args = parse_args(argv)

    ### taxa ranks
    levels = [s.strip() for s in args.ranks.split(",")]
    sam_file_name = args.sam
    outfile_name = args.outfile or (sam_file_name + '.blca.out')

    ### read in pre-formatted lineage information ###
    na_handler = NotAvailableHandler()
    acc2tax = read_tax_acc(args.tax, na_handler, levels)
    print("> 1 > Read in taxonomy information!")

    reference_sequences = {}
    with open(args.reference) as f:
        for r in SeqIO.parse(f, "fasta"):
            reference_sequences[r.id] = str(r.seq)

    print("> 2 > Read in reference db")

    SequenceInfo = namedtuple('SequenceInfo', ['seq', 'hits'])
    ### read in input fasta file ###
    input_sequences = {}
    possible_rejects = set()
    with open(sam_file_name) as sam_file:
        for line in sam_file:
            pieces = line.strip().split('\t')
            entry = SamEntry(pieces)

            if entry.identity_ratio_old < args.iset:
                possible_rejects.add(entry.qname)
            if entry.soft_clipped_ratio > args.softclipping:
                possible_rejects.add(entry.qname)
            elif entry.rname not in reference_sequences:
                possible_rejects.add(entry.qname)
            elif len(reference_sequences[entry.rname]) / float(len(entry.seq)) < args.length:
                possible_rejects.add(entry.qname)
            elif entry.qname not in input_sequences:
                input_sequences[entry.qname] = SequenceInfo(seq=entry.seq, hits=[entry.rname])
            else:
                input_sequences[entry.qname].hits.append(entry.rname)

    rejects = possible_rejects.difference(set(input_sequences))

    print("> 3 > Read in bowtie2 output!")

    settings = {
        "levels": levels,
        "muscle": args.muscle,
        "maxiters": args.muscle_max_iterations,
        "use_diags": args.muscle_use_diags,
        "match": args.match,
        "mismatch": args.mismatch,
        "ngap": args.ngap,
        "nper": args.nper,
        "seed": args.seed,
    }
    counters = {"count": 0, "self_id_queries": 0}

    # Write to a per-process temp file and rename on success, so a killed run never
    # leaves a partial file behind and two concurrent runs can't interleave writes.
    final_outfile_name = outfile_name
    tmp_outfile_name = f"{final_outfile_name}.tmp{os.getpid()}"
    with open(tmp_outfile_name, 'w') as outfile:
        outfile.write("featureid\ttaxonomy\ttaxonomy_confidence\taccessions\n")
        items = build_tasks(input_sequences, reference_sequences, acc2tax, na_handler,
                            levels, counters)

        def write_result(result):
            warning, out_line = result
            if warning:
                print(warning)
            outfile.write(out_line)

        if args.threads == 1:
            _init_worker(settings)
            for kind, payload in items:
                if kind == "line":
                    outfile.write(payload)
                else:
                    write_result(classify_query(payload))
        else:
            # Sliding window: at most max_pending queries in flight, written in
            # input order as the oldest finishes. Bounds the parent's memory
            # (each task carries its hits' sequences) instead of queueing a whole
            # fold's tasks at once. Lines finished in the parent (no valid hit)
            # wait in the window too, so they keep their place in the output.
            max_pending = args.threads * 16
            window = deque()

            def flush_oldest():
                kind, payload = window.popleft()
                if kind == "line":
                    outfile.write(payload)
                else:
                    write_result(payload.get())

            with multiprocessing.Pool(args.threads, initializer=_init_worker,
                                      initargs=(settings,)) as pool:
                for kind, payload in items:
                    if kind == "line":
                        window.append(("line", payload))
                    else:
                        window.append(("task", pool.apply_async(classify_query, (payload,))))
                    while len(window) > max_pending:
                        flush_oldest()
                while window:
                    flush_oldest()

        # sorted: set order changes per process (string hashing is randomised)
        for seqn in sorted(rejects):
            outfile.write(seqn + "\tUnclassified\n")

    if counters["self_id_queries"]:
        print(
            f"> Note: {counters['self_id_queries']} of {counters['count']} queries share an ID with the "
            "reference database. This is expected for self-validation; each query "
            "was still scored against its own reference record. If you did not "
            "intend to classify reference sequences against themselves, check that "
            "your query IDs are distinct from the reference accessions."
        )

    os.replace(tmp_outfile_name, final_outfile_name)


if __name__ == "__main__":
    main()
