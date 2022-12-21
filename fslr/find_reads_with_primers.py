import os
import pysam
from collections import defaultdict
from skbio.alignment import StripedSmithWaterman
import subprocess
import sys

"""
Try and pull out reads with the expected primer sequences
"""


def rev_comp(s):
    d = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join([d[i] for i in s])[::-1]


def check_primer(primer_pairs, seq):
    p1 = False
    p2 = False
    m = {}  # memoize alignments
    for primer1, primer2, p1name, p2name, strand1, strand2 in primer_pairs:
        # check ends of input for primer pair
        # 0.4 = original thresh, 0.7 for high stringency
        target_score1 = int((len(primer1) * 2) * 0.6)
        target_score2 = int((len(primer2) * 2) * 0.6)
        if p1name + strand1 in m:
            p1 = m[p1name + strand1]
        else:
            # check alignment
            p1 = StripedSmithWaterman(primer1)(seq[:len(primer1)*4]).optimal_alignment_score >= target_score1
            if p1:
                p1 = p1name + strand1
            m[p1name + strand1] = p1
        if not p1:
            continue
        if p2name + strand2 in m:
            p2 = m[p2name + strand2]
        else:
            p2 = StripedSmithWaterman(primer2)(seq[-len(primer2)*4:]).optimal_alignment_score >= target_score2
            if p2:
                p2 = p2name + strand2
            m[p2name + strand2] = p2
        if p1 and p2:
            break

    return tuple(sorted([str(p1), str(p2)]))


def run_porechop(samp, out):
    subprocess.run(
        f"/home/kez/Tools/Porechop/porechop-runner.py --discard_middle -t1 -i {samp} -o {out}",
                   shell=True, stdout=subprocess.PIPE)


def filter_fastq(pth, key_to_p, key_to_p_rev, basename):
    fq = pysam.FastqFile(pth)
    primer_pairs = []
    for k in key_to_p:
        forward = key_to_p[k]
        for k2 in key_to_p:
            rev = key_to_p_rev[k2]
            primer_pairs.append((forward, rev, k, k2, "F", "R"))

    total = 0
    counts = defaultdict(int)
    with open(f"{basename}.primers_labelled.fq", "w") as out, open(f"{basename}.no_primers.fq", "w") as out2:
    # out = open(f"{basename}.primers_labelled.fq", "w")
    # out2 = open(f"{basename}.no_primers.fq", "w")
        recs = []
        for c, aln in enumerate(fq):
            seq = aln.sequence
            key = check_primer(primer_pairs, seq)
            recs.append({"p": f"{key[0]}_{key[1]}", "l": len(aln.sequence), "name": aln.name})

            # aln.name = f"{aln.name}.{key[0]}_{key[1]}"
            # out.write(str(aln) + "\n")

            if key != ('False', 'False'):
                aln.name = f"{aln.name}.{key[0]}_{key[1]}"
                out.write(str(aln) + "\n")
            else:
                aln.name = f"{aln.name}.{key[0]}_{key[1]}"
                out2.write(str(aln) + "\n")
            counts[key] += 1
            total += 1

    return counts, recs


def func(args):
    path, key_to_p, lock, filter_counts, keep_temp = args

    basename = path.replace('.filtered_junk.fq', '')
    key_to_p_rev = {k: rev_comp(v) for k, v in key_to_p.items()}

    run_porechop(path, f'{basename}.porechop.fq')
    counts, recs = filter_fastq(f'{basename}.porechop.fq', key_to_p, key_to_p_rev, basename)

    with lock:
        print('find_reads_with_primers counts:', path, dict(counts), file=sys.stderr)
        for k, v in counts.items():
            filter_counts['_'.join(k)] += v

    if not keep_temp:
        os.remove(f'{basename}.filtered_junk.fq')
        os.remove(f'{basename}.porechop.fq')
