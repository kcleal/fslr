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
    # Minimum sequence length threshold - drop sequences shorter than 200bp
    if len(seq) < 200:
        return '_short'

    # Define the end trim size - either 100bp or less if the sequence is short
    trim_size = 100
    for k in target_primers:
        for s in (primers[k], primers_r[k]):

            trim = seq[trim_size:len(seq) - trim_size]
            if not trim:
                return '_short'

            # Check alignment with chunking for long sequences
            if len(trim) > 10000:
                seq_len = len(trim)
                start = 0
                end = 10000
                chunk_overlap = len(
                    s) + 10  # Add a bit more overlap to ensure we don't miss alignments at chunk boundaries

                while start < seq_len:  # Changed condition to be more explicit
                    # Calculate chunk boundaries with overlap
                    chunk_start = max(0, start - chunk_overlap if start > 0 else 0)
                    chunk_end = min(end + chunk_overlap if end < seq_len else seq_len, seq_len)

                    sub = trim[chunk_start:chunk_end]
                    aln = StripedSmithWaterman(s)(sub)
                    if aln.optimal_alignment_score >= 28:
                        return '_concatemer'

                    # Move to next chunk
                    if end >= seq_len:
                        break
                    start += 10000
                    end += 10000
            else:
                # For shorter sequences, analyze the whole trimmed sequence at once
                aln = StripedSmithWaterman(s)(trim)
                if aln.optimal_alignment_score >= 28:
                    return '_concatemer'

    return ''


def telmer_pct(rot, s):
    telmer_count = 0
    tot = 0
    for kmer in (s[ii:ii + 6] for ii in range(len(s) - 6 + 1)):
        if kmer in rot:
            telmer_count += 1
        tot += 1
    return telmer_count / tot


def get_seqs_to_drop(primer_list, primers, primers_r, outfile, filter_counts_all, tantanfile, lock, rot,
                     junkfile):

    length = 150

    f = pysam.FastxFile(tantanfile)

    bad = set([])

    filter_counts = {'total_kept': 0, 'concatemers_dropped': 0, 'total_dropped': 0, 'junk_seqs_dropped': 0,
                     'short_seqs_dropped': 0}

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
                if junkfile:
                    l.name += '_junk'
                break
        else:
            reason = check_for_concatemer(seq, primer_list, primers, primers_r)
            if reason:
                if reason == '_short':
                    filter_counts['short_seqs_dropped'] += 1
                elif reason == '_concatemer':
                    filter_counts['concatemers_dropped'] += 1
                else:
                    raise RuntimeError(reason)

                # filter_counts['concatemers_dropped'] += 1
                drop = True
                if junkfile:
                    l.name += reason

        if not drop:
            filter_counts['total_kept'] += 1
            outfile.write(str(l) + '\n')
        else:
            filter_counts['total_dropped'] += 1
            if junkfile:
                junkfile.write(str(l) + '\n')


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
    with open(f'{outfolder}/{outname}.{temp_name}.filtered_junk.fq', 'w') as outfile, \
         open(f'{outfolder}/{outname}.{temp_name}.junk.fq', 'w') as junkfile:
        primers_r = {k: rev_comp(v) for k, v in primers.items()}
        if fq_input[-3:] == ".gz":
            subprocess.run(f'gzip -dc {fq_input} | tantan - > {outfolder}/tmp.tantan.{temp_name}.fasta', shell=True)
        else:
            subprocess.run(f'tantan {fq_input} > {outfolder}/tmp.tantan.{temp_name}.fasta', shell=True)
        tantanfile = glob.glob(f'{outfolder}/tmp.tantan.{temp_name}.fasta')[0]
        primer_list = primers.keys()
        fc = get_seqs_to_drop(primer_list, primers, primers_r, outfile, filter_counts_all, tantanfile, lock, rot,
                              junkfile if keep_temp else None)
        if not keep_temp:
            os.remove(f'{outfolder}/tmp.tantan.{temp_name}.fasta')
            os.remove(f'{outfolder}/{outname}.{temp_name}.junk.fq')
        with lock:
            print(fq_input, file=sys.stderr)
            print(fc, file=stderr)
