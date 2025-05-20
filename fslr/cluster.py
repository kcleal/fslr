import pandas as pd
import networkx as nx
import numpy as np
import pysam
# from sortedintersect import IntervalSet
from superintervals import IntervalSet
from collections import defaultdict, namedtuple


IntervalItem = namedtuple('interval_item',
                          ['chrom', 'start', 'end', 'aln_size', 'qname', 'n_alignments', 'qlen2', 'middle', 'index'])


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


def rename_chromosomes(bed_file, chromosome_lengths, chromosome_mask):
    chromosome_names = sorted(set(bed_file['chrom'].unique().tolist()),
                              key=lambda x: int(x[3:]) if x[3:].isdigit() and x[:3] == 'chr' else float('inf'))
    chromosome_to_numeric_map = {name: i + 1 for i, name in enumerate(chromosome_names)}

    chr_lengths = {chromosome_to_numeric_map.get(k): v for k, v in chromosome_lengths.items()}
    bed_file['chrom'] = bed_file['chrom'].apply(chromosome_to_numeric_map.get)
    chromosome_mask = [chromosome_to_numeric_map.get(x) if x != 'subtelomere' else x for x in chromosome_mask]

    return bed_file, chr_lengths, chromosome_mask, chromosome_to_numeric_map


def chrom_to_str(bed_df, chromosome_to_numeric_map):
    num_to_string_map = {value: key for key, value in chromosome_to_numeric_map.items()}
    bed_df['chrom'] = bed_df['chrom'].map(num_to_string_map)
    return bed_df


def calc_coverage(bed_file, chromosome_lengths):
    coverage = {}

    grouped_by_chr = bed_file.groupby('chrom')

    for chr, group in grouped_by_chr:
        if chr not in chromosome_lengths:
            continue
        c = np.zeros(chromosome_lengths[chr] + 1)

        np.add.at(c, group['rstart'].values, 1)
        np.add.at(c, group['rend'].values, -1)

        coverage[chr] = np.cumsum(c)

    return coverage


def filter_high_coverage(data, bed_file, chromosome_lengths, threshold):
    cov = calc_coverage(bed_file, chromosome_lengths)
    new_alignments = []
    for aln in data:
        if cov[aln.chrom][aln.middle] > threshold:
            continue
        new_alignments.append(aln)
    return new_alignments


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
    chromosome_lengths = {key: value for key, value in chromosome_lengths.items() if value > 1000000}
    for a in read_alignments:
        if a.chrom in mask:
            continue
        if 'subtelomere' in mask:
            if a.chrom in chromosome_lengths and \
               (a.start < threshold or chromosome_lengths[a.chrom] - a.end < threshold):
                continue

        new_alignments.append(a)
    if len(read_alignments) == 1 and before >= 4:
        return []
    return new_alignments


def prepare_data(bed_df, cluster_mask, chromosome_lengths, threshold=500_000):
    # need to make sure rend > rstart for sortedintersect, and intervals are sorted
    bed_df['start'] = np.minimum(bed_df['rstart'], bed_df['rend'])
    bed_df['end'] = np.maximum(bed_df['rstart'], bed_df['rend'])
    bed_df['middle'] = bed_df['aln_size'] // 2 + bed_df['start']
    bed_df = bed_df.sort_values('start')
    columns = ['chrom', 'start', 'end', 'aln_size', 'qname', 'n_alignments', 'qlen2', 'middle']
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
    for k, v in interval_tree.items():
        v.index()
    return interval_tree


def calculate_overlap(interval1, interval2):
    overlap = max(0, (min(interval1.end, interval2.end) - max(interval1.start, interval2.start)))
    reciprocal_overlap = min(overlap / interval1.aln_size, overlap / interval2.aln_size)
    return reciprocal_overlap


# see the sizes of comparisons
def overall_jaccard_similarity(l1, l2, l2_comparisons, percentage, min_threshold):

    if not l1 or not l2:
        return 0, 0
    len1 = len(l1)
    len2 = len(l2)
    len_product = len1 * len2

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
                l2_comparisons[j] = 1
                intersection += 1
                zeros -= 2
                break
            if count / len_product < 1 - min_threshold and intersection == 0:
                continue

    union = intersection + zeros

    if union == 0:
        return 0, 0

    return intersection / union, intersection


def get_chromosome_lengths(bam_path):
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    return {bam_file.get_reference_name(tid): l for tid, l in enumerate(bam_file.lengths)}


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
    l2_comparisons = np.zeros(100000)

    for query_key, list1 in query_intervals.items():
        edges = 0
        for itv in list1:
            # overlap_intervals = interval_trees[itv.chrom].search_interval(itv.start, itv.end)
            overlap_intervals = interval_trees[itv.chrom].find_overlaps(itv.start, itv.end)
            for o_data in overlap_intervals:
                if o_data.qname == query_key:
                    continue
                b = tuple(sorted((o_data.qname, query_key)))
                if b in seen_edges:
                    continue
                seen_edges.add(b)
                if different_lengths_or_alignments(itv, o_data, qlen_diff, diff):
                    continue

                # add counter
                list2 = query_intervals[o_data.qname]

                j, n_i = overall_jaccard_similarity(list1, list2, l2_comparisons, overlap_cutoff, min_threshold)
                if n_i == 0:
                    continue
                target = jaccard_threshold[n_i - 1] if n_i - 1 < len(jaccard_threshold) else jaccard_threshold[-1]
                if j >= target:
                    match.add((query_key, o_data.qname, j))
                    G.add_edge(query_key, o_data.qname)
                    edges += 1
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
