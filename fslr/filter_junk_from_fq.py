# goal is to remove repetitive sequences from nanopore data - 'junk' sequences and concatemers
# concatemers are identified by primer sites within the main sequence
import glob
import subprocess
import pysam
from math import exp
from skbio.alignment import StripedSmithWaterman
import sys
import uuid
import os
from sys import stderr
from collections import deque


def find_lower_case(s):
    i = 0
    while i < len(s):
        if s[i].islower():
            end = i + 1
            for j in range(end, len(s)):
                if not s[j].islower():
                    break
                end += 1
            yield i, end
            i = end
        else:
            i += 1


def compute_rep(seq):

    last_visited = {}
    tot_amount = 0
    total_seen = 0

    for k in (2, 3, 4, 5, 6):

        decay = 0.25 * 1/k
        max_amount = exp(-decay) * k  # If last kmer was the same as current kmer

        sub_ptr = 0
        for i in range(len(seq) - k):

            a = seq[sub_ptr:sub_ptr + k]
            if a in last_visited:
                diff = i - last_visited[a]
                x = exp(-decay * diff)
                amount = (k * x) / max_amount

            else:
                amount = 0
            if i > k:
                tot_amount += amount
                total_seen += 1
            last_visited[a] = i
            sub_ptr += 1

    if total_seen == 0:
        return 0

    return tot_amount / total_seen


def check_for_concatemer(seq, target_primers, primers, primers_r):

    for k in target_primers:
        for s in [primers[k], primers_r[k]]:
            ls = len(s) * 4
            if len(seq) < 4 * ls:
                return True  # drop anyway
            trim = seq[ls:len(seq)-ls]
            # check alignment
            if len(trim) > 10000:
                seq_len = len(trim)
                start = 0
                end = 10_000
                while True:
                    # include small overlap
                    sub = trim[max(0, start-len(s) - 1):min(end, seq_len)]
                    aln = StripedSmithWaterman(s)(sub)
                    if aln.optimal_alignment_score >= 28:
                        return True
                    start += 10_000
                    end += 10_000
                    if start >= seq_len:
                        break
            else:
                aln = StripedSmithWaterman(s)(trim)
                if aln.optimal_alignment_score >= 28:
                    return True

    return False


def telmer_pct(rot, s):
    telmer_count = 0
    tot = 0
    for kmer in (s[ii:ii + 6] for ii in range(len(s) - 6 + 1)):
        if kmer in rot:
            telmer_count += 1
        tot += 1
    return telmer_count / tot


def get_seqs_to_drop(fq_input, primer_list, primers, primers_r, outfile, filter_counts_all, tantanfile, lock, rot):

    length = 150

    f = pysam.FastxFile(tantanfile)

    bad = set([])
    name = fq_input.split('/')[-1].split('.')[0]

    filter_counts = {'total_kept': 0, 'concatemers_dropped': 0, 'total_dropped': 0, 'junk_seqs_dropped': 0}

    for l in f:
        seq = l.sequence

        # find blocks of lowercase
        drop = False
        for start, end in find_lower_case(seq):
            if end - start > length:
                s = seq[start:end].upper()

                pct_tel = telmer_pct(rot, s)

                if pct_tel > 0.3:
                    continue
                rep = compute_rep(s)
                if rep < 0.3:
                    continue

                bad.add(l.name)
                filter_counts['junk_seqs_dropped'] += 1
                drop = True
                break
        else:
            concat = check_for_concatemer(seq, primer_list, primers, primers_r)
            if concat:
                filter_counts['conactemers_dropped'] += 1
                drop = True

        if not drop:
            filter_counts['total_kept'] += 1
            outfile.write(str(l) + '\n')
        else:
            filter_counts['total_dropped'] += 1

    with lock:
        for k, v in filter_counts.items():
            if k in filter_counts_all:
                filter_counts_all[k] += v
    return filter_counts


def rev_comp(s):
    d = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join([d[i] for i in s])[::-1]


def tel_tokens(tel,):
    d = deque(tel)
    f_rotations = []
    for i in range(len(tel)-1):
        d.rotate()
        f_rotations.append("".join(d))
    return f_rotations


def make_rotations(tta):
    variant_rotations = set([])
    for key in tta:
        variant_rotations.update(tel_tokens(key))
    return variant_rotations


def func(args):
    targets = ["CCCTAA", "CCCTGA", "CCCGAA", "CCCTAC", "CCCTCA", "CCCCAA", "CCCTTA", "CCCTAT", "CCCTAG", "CCCAAA", "CCCACT", "CCCCAT", "CCCGCA", "CCCGCT", "CCCTCT"]
    targets += [rev_comp(t) for t in targets]

    rot = make_rotations(targets)

    temp_name = str(uuid.uuid4())
    fq_input, primers, outfolder, outname, filter_counts_all, lock, keep_temp = args
    with open(f'{outfolder}/{outname}.{temp_name}.filtered_junk.fq', 'w') as outfile:
        primers_r = {k: rev_comp(v) for k, v in primers.items()}
        if fq_input[-3:] == ".gz":
            subprocess.run(f'gzip -dc {fq_input} | tantan - > {outfolder}/tmp.tantan.{temp_name}.fasta', shell=True)
        else:
            subprocess.run(f'tantan {fq_input} > {outfolder}/tmp.tantan.{temp_name}.fasta', shell=True)
        tantanfile = glob.glob(f'{outfolder}/tmp.tantan.{temp_name}.fasta')[0]
        primer_list = primers.keys()
        fc = get_seqs_to_drop(fq_input, primer_list, primers, primers_r, outfile, filter_counts_all, tantanfile, lock, rot)
        if not keep_temp:
            os.remove(f'{outfolder}/tmp.tantan.{temp_name}.fasta')

        with lock:
            print(fq_input, file=sys.stderr)
            print(fc, file=stderr)
