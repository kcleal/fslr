import numpy as np
import pandas as pd
import pysam
import pandas
import sys
import networkx as nx
from networkx import community, connected_components
import pafpy
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import subprocess
import glob
import shutil
import re
from skbio.alignment import StripedSmithWaterman
import multiprocessing
from sklearn.ensemble import IsolationForest
from pywfa.align import WavefrontAligner

from networkx.algorithms import approximation


def find_targets(filter_counts, path, cluster_fraction='auto'):
    bam = pysam.AlignmentFile(path, 'rb')
    basename = path.replace('.bwa_dodi.bam', '')

    allowed_qnames = set(idx for idx, grp in pd.read_csv(f"{basename}.mappings.bed", delimiter='\t').groupby('qname') if len(grp) > 1)

    lengths = []
    qnames = []
    for r in bam.fetch(until_eof=True):
        if r.flag & 2304 and r.qname in allowed_qnames:  # not (not primary, supplementary) --> is primary
            qnames.append(r.qname)
            lengths.append(r.infer_query_length())

    bins = np.arange(0, 15_000, 50)
    lens = np.array(lengths)[:, np.newaxis]
    clf = IsolationForest(random_state=0, contamination=cluster_fraction).fit(lens)
    pred = clf.predict(lens)
    outliers = [l for l, v in zip(lengths, pred) if v != 1]
    inliers = [l for l, v in zip(lengths, pred) if v == 1]
    print('N outlier reads', len(outliers))
    filter_counts['n_cluster_targets'] = len(outliers)

    plt.figure(figsize=(6, 4))
    plt.ylabel('Cound')
    plt.xlabel('Read length (bp)')
    plt.title(basename.split("/")[-1])
    plt.yscale('log')
    plt.hist(inliers, bins=bins, label='Skip-cluster')
    plt.hist(outliers, bins=bins, alpha=0.5, color='red', label='Cluster targets')
    plt.savefig(basename + 'cluster_targets_log.pdf')
    plt.close()
    plt.figure(figsize=(6, 4))
    plt.ylabel('Count')
    plt.xlabel('Read length (bp)')
    plt.title(basename.split("/")[-1])
    plt.hist(inliers, bins=bins, label='Skip-cluster')
    plt.hist(outliers, bins=bins, alpha=0.5, color='red', label='Cluster targets')
    plt.savefig(basename + 'cluster_targets.pdf')
    plt.close()

    outlier_names = set(q for q, v in zip(qnames, pred) if v != 1)
    bam = pysam.AlignmentFile(path, 'rb')
    with open(basename + '.cluster_targets.fasta', 'w') as fq:
        for r in bam:
            if not r.flag & 2304: # and r.qname in outlier_names:
                fq.write(f'>{r.qname}\n{r.seq}\n')

    return


def merge_intervals(interval):
    temp_tuple = sorted(interval)
    merged = [temp_tuple[0]]
    for current in temp_tuple:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


def intervals_sum(intervals):
    return sum(j - i for i, j in intervals)


def gaps_in_cigar(cigar):
    c = re.split(r'(\d+)', cigar)[1:]  # Drop leading empty string
    aligned_bases = 0
    gaps = False
    for i in range(0, len(c) - 1, 2):
        if c[i + 1] not in "DI":
            if c[i + 1] == 'M':
                aligned_bases += int(c[i])
            continue  # Don't count deletions, or soft/hard clips at right-hand side, or X mismatches
        if int(c[i]) >= 30:
            gaps = True
    return aligned_bases, gaps


def proc_job(args):
    v_qseq, tseq, u, v = args
    #
    a = WavefrontAligner(v_qseq)(tseq)
    wfa_score = abs(a.score / max(len(v_qseq), len(tseq)))
    aligned_bases, gaps = gaps_in_cigar(a.cigarstring)
    w = aligned_bases / len(v_qseq)
    # if gaps or w < 0.9:
    if gaps or wfa_score > 0.4:
        return False

    return u, v, (1 - wfa_score)**2 #w

    # old method:
    query = StripedSmithWaterman(v_qseq, suppress_sequences=True)
    res = query(tseq)
    aligned_bases, gaps = gaps_in_cigar(res.cigar)
    w = aligned_bases / len(v_qseq)
    if gaps or w < 0.98:
        return False
    return u, v, w


def cluster_paf(basename, procs):

    old_chars = "ACGTN"
    replace_chars = "TGCAN"
    tab = str.maketrans(old_chars, replace_chars)
    name = basename.split('/')[-1]
    print('Clustering ', basename, file=sys.stderr)

    # get lists of candidate matching reads
    # candidates = defaultdict(set)
    G = nx.Graph()
    c = 0
    with pafpy.PafFile(basename + '.minimap2_cluster.paf') as paf:
        for r in paf:
            if r.tname != r.qname:
                #candidates[r.qname].add(r.tname)
                G.add_edge(r.qname, r.tname)
                c += 1
    print('Candidate alignments', c)

    pool = multiprocessing.Pool(procs)
    name_seq = {}
    jobs = []
    af = pysam.AlignmentFile(basename + '.bwa_dodi.bam', 'rb')
    for a in af:
        if a.flag & 2304:
            continue
        s = a.get_forward_sequence()
        if not s:
            continue
        name_seq[a.qname] = s


    # par keys:
    # 'blast_identity', 'blen', 'count', 'from_str', 'get_tag', 'index', 'is_inversion', 'is_primary', 'is_secondary',
    # 'is_unmapped', 'mapq', 'mlen', 'qend', 'qlen', 'qname', 'qstart', 'query_aligned_length', 'query_coverage',
    # 'relative_length', 'strand', 'tags', 'target_aligned_length', 'target_coverage', 'tend', 'tlen', 'tname', 'tstart'


    # G2 = nx.Graph()
    # paf_alignments = {}
    # for k, v in candidates.items():
    #     tnames = defaultdict(list)
    #     for r in v:  # enumerate tnames
    #         tnames[r.tname].append(r)
    #     for kk, vv in tnames.items():
    #         for i in vv:
    #             # qc = i.query_coverage
    #             # tc = i.target_coverage
    #             # if tc > min_cov and qc > min_cov:
    #                 paf_alignments[(i.qname, i.tname)] = i
    #                 # weight = min(tc, qc) / max(tc, qc)
    #                 # G.add_edge(k, kk, weight=i.blast_identity())
    #                 G.add_edge(k, kk)
    #                 # G2.add_edge(k, kk, weight=i.blast_identity())
    #                 # break

    print('Finding connected components', file=sys.stderr)
    # do second round of clustering within connected components. all-vs-all using global alignment

    G2 = nx.Graph()
    com = sorted(connected_components(G), key=lambda x: -len(x))
    for clst_id, dta in enumerate(com):
        seen = set([])
        print(clst_id, len(dta))
        jobs = []
        for u in dta:
            for v in dta:
                if (u, v) in seen or (v, u) in seen:
                    continue
                # if (u, v) in paf_alignments:
                #     i = paf_alignments[(u, v)]
                # elif (v, u) in paf_alignments:
                #     i = paf_alignments[(v, u)]
                # else:
                #     continue
                seen.add((u, v))
                # if i.qname in name_seq:
                #     qseq = name_seq[i.qname]
                # else:
                #     continue
                # if i.tname in name_seq:
                #     tseq = name_seq[i.tname]
                # else:
                #     continue
                # if not qseq or not tseq:
                #     continue

                qseq = name_seq[u]
                tseq = name_seq[v]

                v_qseq = qseq.translate(tab)[::-1]

                jobs.append((qseq, tseq, u, v))
                jobs.append((v_qseq, tseq, u, v))
        print(len(jobs))
        with multiprocessing.Pool(procs) as p:
            res = p.map(proc_job, jobs)
            for item in res:
                if item:
                    u, v, w = item
                    G2.add_edge(u, v, weight=w)


    if G2.size() == 0:
        print('Size of community graph was 0. No clusters found sorry.', file=sys.stderr)
        return
    print('Edges in cluster graph before community detection', len(G2.edges()), file=sys.stderr)
    com = community.greedy_modularity_communities(G2, weight='weight')
    # com = sorted(connected_components(G), key=lambda x: -len(x))
    print('Cluster_ids, N_reads', file=sys.stderr)
    names_2_id = {}
    cluster_counts = defaultdict(int)
    for clst_id, dta in enumerate(com):

        sub = G2.subgraph(dta)
        cluster_counts[len(dta)] += 1

        print(clst_id, len(dta), approximation.average_clustering(sub))
        print()

        if len(dta) == 1:
            continue

        for k in dta:
            names_2_id[k] = clst_id

    print(cluster_counts)

    bam = pysam.AlignmentFile(basename + '.bwa_dodi.bam', 'rb')

    if os.path.exists(basename + '_clusters'):
        shutil.rmtree(basename + '_clusters')
    os.mkdir(basename + '_clusters')

    outfs = [pysam.AlignmentFile(basename + f'_clusters/{name}.{cid}.{len(dta)}.bam', 'wb', template=bam) for cid, dta in enumerate(com)]
    singles = 0
    cids = []
    for a in bam.fetch(until_eof=True):
        if a.qname in names_2_id:
            cid = names_2_id[a.qname]
            outfs[cid].write(a)
            cids.append((cid, a.qname))
        else:
            singles += 1

    with open(f"{basename}.clusterIDs.txt", "w") as cl_out:
        for i in sorted(cids):
            cl_out.write(f'{i[0]}\t{i[1]}\n')

    print('n single read clusters', singles, file=sys.stderr)

    for item in outfs:
        item.close()

    bam.close()
    for f in glob.glob(f'./{name}_clusters/*.bam'):
        subprocess.run(f'samtools index {f}', shell=True)

    for cid, dta in enumerate(com):
        b = basename + f'_clusters/{name}.{cid}.{len(dta)}.bam'
        c = f"samtools fasta {b} > {basename}_clusters/{name}.{cid}.cluster_reads.fa"
        subprocess.run(c, shell=True)
        c = f"abpoa {basename}_clusters/{name}.{cid}.cluster_reads.fa | sed 's/Consensus_sequence/{name}.{cid}/g' > {basename}_clusters/{name}.{cid}.cons.fa"
        subprocess.run(c, shell=True)
        subprocess.run(f'rm {b}', shell=True)

    print('Clustering done', file=sys.stderr)
