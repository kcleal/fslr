import os
import pysam
from collections import defaultdict
from skbio.alignment import StripedSmithWaterman
import sys

"""
Try and pull out reads with the expected primer sequences
"""


def rev_comp(s):
    d = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join([d[i] for i in s])[::-1]


def check_primer(primer_pairs, seq):  # legacy
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


def check_primer2(primer_pairs, read, trim_thresh):
    seq = read.sequence
    res = []
    ss = 500
    for primer1, primer2, p1name, p2name, strand1, strand2 in primer_pairs:
        max_score1 = len(primer1) * 2
        max_score2 = len(primer2) * 2
        p1_space = min(int(len(seq)/2), ss)
        p2_space = min(int(len(seq)/2), ss)
        aln1 = StripedSmithWaterman(primer1, suppress_sequences=True)(seq[:p1_space])
        aln2 = StripedSmithWaterman(primer2, suppress_sequences=True)(seq[-p2_space:])
        score1 = aln1.optimal_alignment_score / max_score1
        score2 = aln2.optimal_alignment_score / max_score2
        name1 = 'False' if score1 < trim_thresh else p1name + strand1
        name2 = 'False' if score2 < trim_thresh else p2name + strand2
        res.append((round(score1, 2), round(score2, 2), name1, name2, aln1, aln2, p1_space, p2_space))

    best = sorted(res, key=lambda x: (x[0] + x[1]))[-1]
    if best[2] == 'False' and best[3] == 'False':
        return best[0], best[1], best[2], best[3], 0

    target_begin = best[4].target_begin
    target_end = len(seq) - best[7] + best[5].target_end_optimal
    trimmed = target_begin + (len(seq) - target_end)
    read.sequence = read.sequence[target_begin:target_end]
    read.quality = read.quality[target_begin:target_end]
    return best[0], best[1], best[2], best[3], trimmed


def label_and_chop_primers(pth, key_to_p, basename, trim_thresh):
    fq = pysam.FastqFile(pth)
    primer_pairs = set([])
    for k1 in key_to_p:
        k1_forward = key_to_p[k1]
        k1_reverse = rev_comp(k1_forward)
        for k2 in key_to_p:
            k2_forward = key_to_p[k2]
            k2_reverse = rev_comp(k2_forward)
            primer_pairs.add((k1_forward, k2_reverse, k1, k2, "F", "R"))
            if k1 != k2:
                primer_pairs.add((k1_reverse, k2_forward, k1, k2, "R", "F"))

    total = 0
    counts = defaultdict(int)
    counts['starting_bases'] = 0
    counts['trimmed_bases'] = 0
    counts['trimmed_reads'] = 0

    with open(f"{basename}.primers_labelled.fq", "w") as out, open(f"{basename}.no_primers.fq", "w") as out2:
        recs = []
        for c, aln in enumerate(fq):
            counts['starting_bases'] += len(aln.sequence)
            key = check_primer2(primer_pairs, aln, trim_thresh)
            counts['trimmed_bases'] += key[4]
            if key[4] > 0:
                counts['trimmed_reads'] += 1
            recs.append({"p": f"{key[0]}_{key[1]}", "l": len(aln.sequence), "name": aln.name})
            if key[2] != 'False' or key[3] != 'False':
                aln.name = f"{aln.name}.{key[0]}_{key[1]}.{key[2]}_{key[3]}"
                out.write(str(aln) + "\n")
            else:
                aln.name = f"{aln.name}.{key[0]}_{key[1]}.{key[2]}_{key[3]}"
                out2.write(str(aln) + "\n")
            counts[f'{key[2]}_{key[3]}'] += 1
            total += 1
    # print('Total processed: ', total, 'Trimmed bases: ', counts['trimmed_bases'], file=sys.stderr)
    return counts, recs


def func(args):

    path, key_to_p, lock, filter_counts, keep_temp, trim_thresh = args
    basename = path.replace('.filtered_junk.fq', '')
    counts, recs = label_and_chop_primers(f'{basename}.filtered_junk.fq', key_to_p, basename, trim_thresh)
    for k, v in counts.items():
        if k not in filter_counts:
            filter_counts[k] = 0
        filter_counts[k] += v
    # print(f'find_reads_with_primers counts', filter_counts, file=sys.stderr)
    if not keep_temp:
        os.remove(f'{basename}.filtered_junk.fq')
