import pandas as pd
import networkx as nx
import numpy as np
import pysam
from sortedintersect import IntervalSet
from collections import defaultdict, namedtuple

IntervalItem = namedtuple('interval_item',
                          ['chrom', 'start', 'end', 'aln_size', 'qname', 'n_alignments', 'qlen2', 'index'])


# mask high coverage regions - mosdepth algorithm
# def mosdepth_coverage():

def keep_fillings(bed_file):
    first = {}
    last = {}
    for idx, qname in zip(bed_file.index, bed_file['qname']):
        if qname not in first:
            first[qname] = idx
            last[qname] = idx
        else:
            last[qname] = idx
    to_drop = list(first.values()) + list(last.values())
    bed_file = bed_file[~bed_file.index.isin(to_drop)]

    qlen2 = {}
    for qname, grp in bed_file.groupby('qname'):
        qlen2[qname] = grp['qend'].max() - grp['qstart'].min()
    bed_file['qlen2'] = [qlen2[q] for q in bed_file['qname']]

    return bed_file


def delete_false(bed_file):
    # Define the string to search for
    string_to_search = 'False'
    to_keep = ~bed_file['qname'].str.contains(string_to_search)
    filtered_df = bed_file[to_keep]

    return filtered_df



def mask_sequences2(read_alignments, mask, chromosome_lengths, threshold=500_000):
    if not mask:
        return read_alignments
    new_alignments = []
    before = len(read_alignments)
    for a in read_alignments:
        if 'subtelomere' in mask:
            if a.start < threshold or (
                    a.chrom in chromosome_lengths and chromosome_lengths[a.chrom] - a.end < threshold):
                continue
        if a.chrom in mask:
            continue
        new_alignments.append(a)
    if len(read_alignments) == 1 and before >= 4:
        return []
    return new_alignments


def prepare_data(bed_df, cluster_mask, chromosome_lengths, threshold=500_000):
    # need to make sure rend > rstart for sortedintersect, and intervals are sorted
    bed_df['start'] = np.minimum(bed_df['rstart'], bed_df['rend'])
    bed_df['end'] = np.maximum(bed_df['rstart'], bed_df['rend'])
    bed_df = bed_df.sort_values('start')
    # calculate interval_sizes
    columns = ['chrom', 'start', 'end', 'aln_size', 'qname', 'n_alignments', 'qlen2']
    data = []
    for i in zip(*(bed_df[col] for col in columns), bed_df.index):
        data.append(IntervalItem(*i))
    if cluster_mask:
        data = mask_sequences2(data, cluster_mask, chromosome_lengths, threshold)
    return data


def build_interval_trees(data):
    interval_tree = defaultdict(lambda: IntervalSet(with_data=True))
    for itv in data:
        interval_tree[itv.chrom].add(itv.start, itv.end, itv)
    return interval_tree


def calculate_overlap(interval1, interval2):
    overlap = max(0, (min(interval1.end, interval2.end) - max(interval1.start, interval2.start)))
    reciprocal_overlap = min(overlap / interval1.aln_size, overlap / interval2.aln_size)
    return reciprocal_overlap


# see the sizes of comparisons
def overall_jaccard_similarity_optimized(l1, l2, l1_comparisons, l2_comparisons, percentage, min_threshold):
    if not l1 or not l2:
        return 0, 0
    len1 = len(l1)
    len2 = len(l2)
    len_product = len1 * len2
    l1_comparisons[:len1] = 0
    l2_comparisons[:len2] = 0
    zeros = len1 + len2
    intersection = 0
    count = 0
    for i, interval1 in enumerate(l1):
        for j, interval2 in enumerate(l2):
            count += 1
            if l2_comparisons[j]:
                continue
            if interval1.chrom == interval2.chrom and calculate_overlap(interval1, interval2) >= percentage:
                l1_comparisons[i] = 1
                l2_comparisons[j] = 1
                intersection += 1
                zeros -= 2
                break
            if count / len_product < 1 - min_threshold and intersection == 0:
                return 0, 0

    union = intersection + zeros

    if union == 0:
        return 0, 0

    return intersection / union, intersection


def get_chromosome_lengths(bam_path):
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    return {bam_file.get_reference_name(tid): l for tid, l in enumerate(bam_file.lengths) if l > 1000000}


def different_lengths_or_alignments(itv1, itv2, qlen_diff, diff):
    if (min(itv1.qlen2, itv2.qlen2) / max(itv1.qlen2, itv2.qlen2)) >= 1 - qlen_diff:
        return False
    if (min(itv1.n_alignments, itv2.n_alignments) / max(itv1.n_alignments, itv2.n_alignments)) >= 1 - diff:
        return False
    return True


def query_interval_trees(interval_trees, data, overlap_cutoff, jaccard_threshold, edge_threshold, qlen_diff, diff):
    min_threshold = min(jaccard_threshold)
    query_intervals = defaultdict(list)
    for itv in data:
        query_intervals[itv.qname].append(itv)
    G = nx.Graph()
    seen_edges = set([])
    match = set([])
    l1_comparisons = np.zeros(100000)
    l2_comparisons = np.zeros(100000)
    for query_key, list1 in query_intervals.items():
        edges = 0
        for itv in list1:
            overlap_intervals = interval_trees[itv.chrom].search_interval(itv.start, itv.end)
            for ol_start, ol_end, o_data in overlap_intervals:
                if o_data.qname == query_key:
                    continue
                b = tuple(sorted((o_data.qname, query_key)))
                if b in seen_edges:
                    continue
                if different_lengths_or_alignments(itv, o_data, qlen_diff, diff):
                    seen_edges.add(b)
                    continue
                # add counter
                list2 = query_intervals[o_data.qname]

                j, n_i = overall_jaccard_similarity_optimized(list1, list2, l1_comparisons, l2_comparisons,
                                                              overlap_cutoff, min_threshold)
                if n_i == 0:
                    continue
                target = jaccard_threshold[n_i - 1] if n_i - 1 < len(jaccard_threshold) else jaccard_threshold[-1]
                if j >= target:
                    match.add((query_key, o_data.qname, j))
                    G.add_edge(query_key, o_data.qname)
                    edges += 1
                seen_edges.add(b)
                if edges >= edge_threshold:
                    break

    match_df = pd.DataFrame(match, columns=['query1', 'query2', 'jaccard_similarity'])
    return match_df, G


def get_subgraphs(G):
    # Identify groups of connected components (or subgraphs)
    sub_graphs = nx.connected_components(G)
    sub_graphs = list(sub_graphs)
    return sub_graphs


def choose_alignment(bed_file):
    qname_grouped = bed_file.groupby('qname')
    avg_scores = qname_grouped['alignment_score'].mean()

    bed_file['avg_alignment_score'] = bed_file['qname'].map(avg_scores)

    cluster_grouped = bed_file.groupby('cluster')

    selected_reads = []

    for cluster_id, group in cluster_grouped:
        # Find the read with the highest average alignment score in the cluster
        max_alignment_read = group.loc[group['avg_alignment_score'].idxmax()]['qname']
        selected_reads.append(max_alignment_read)

    selected_reads_df = bed_file[bed_file['qname'].isin(selected_reads)]

    return selected_reads_df