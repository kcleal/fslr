import sys
import pysam
from collections import defaultdict
import kmerhash
from kmerhash import kmerhasher, hashes2seq
import numpy as np
from scipy.spatial.distance import hamming

def partners(bam_path):

    bam = pysam.AlignmentFile(bam_path, 'rb')

    names = {}
    names_r = {}
    c = 0

    kmer_array = {}

    kmers = defaultdict(set)
    klength = 16
    for a in bam.fetch(until_eof=True):
        if a.qname not in names:
            names[a.qname] = c
            names[c] = a.qname
            c += 1
        k_id = names[a.qname]
        kh = kmerhasher(a.seq, klength)
        kmer_array[a.qname] = kh
        for k in kh:
            kmers[k].add(k_id)

    # now find mappings
    bam = pysam.AlignmentFile(bam_path, 'rb')
    querys = np.zeros(c)
    n = 0
    for a in bam.fetch(until_eof=True):
        if a.flag & 2304:
            continue
        kh = kmer_array[a.qname]
        r_id = names[a.qname]
        print(a.qname)
        for k in kh:
            targets = kmers[k]
            for c_id in targets:
                if c_id != r_id:
                    querys[c_id] += 1
        print(querys.max() / len(kh))
        srt = sorted(enumerate(querys), key=lambda x: x[1])[-10:]

        window_length = 100
        for t_id, score in srt:
            akh = kmer_array[t_id]
            yes = 0
            no = 0
            for w in range(0, len(kh) - window_length + 1):
                q = kh[w]
                for wj in range(w, w+window_length):
                    if wj >= len(akh):
                        no += 1
                        break
                    if akh[wj] == q:
                        yes += 1
                        break
                else:
                    no += 1


            print(names[t_id], score / len(kh), yes/ (yes + no))

        querys[:] = 0
        n += 1

        if n > 10:
            quit()
        print()



    print(c)

