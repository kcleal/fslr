import pandas as pd
from ncls import NCLS
import networkx as nx
import subprocess
import click
import numpy as np
import pysam
import time
from functools import wraps
import sortedcontainers
from sortedintersect import IntervalSet
from collections import defaultdict, namedtuple


IntervalItem = namedtuple('interval_item', ['chrom', 'start', 'end', 'qname', 'n_alignments', 'qlen2', 'index'])


def measure_time(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Execution time of {func.__name__}: {execution_time} seconds")
        return result
    return wrapper


def keep_fillings(bed_file):
    # make a qlen2 column: qlen-len(bread)
    # Identify unique values in the specified column
    bed_file.sort_values(by=['qname', 'qstart'])
    unique_values = bed_file['qname'].unique()

    filtered_df = bed_file.copy()

    for value in unique_values:
        first_occurrence = filtered_df[filtered_df['qname'] == value].index[0]
        last_occurrence = filtered_df[filtered_df['qname'] == value].index[-1]
        filtered_df.loc[first_occurrence:last_occurrence, 'qlen2'] = filtered_df.loc[first_occurrence, 'qlen'] - filtered_df.loc[first_occurrence, 'aln_size'] - filtered_df.loc[last_occurrence, 'aln_size']
        #remove the breads
        filtered_df = filtered_df.drop([first_occurrence, last_occurrence])

    # filter for mapq
    # delete line with mapq<threshold - as default don't delete any of the lines
    # filtered = filtered_df[filtered_df.mapq >= mapq_cutoff]

    # Reset the index of the resulting DataFrame
    filtered_df = filtered_df.reset_index(drop=True)

    return filtered_df


def mask_sequences2(read_alignments, mask, chromosome_lengths, threshold=500_000):
    if not mask:
        return read_alignments
    new_alignments = []
    before = len(read_alignments)
    for a in read_alignments:
        if a.chrom in mask:
            continue
        if 'subtelomere' in mask:
            if a.start < threshold or (a.chrom in chromosome_lengths and chromosome_lengths[a.chrom] - a.end < threshold):
                continue
        new_alignments.append(a)
    if len(read_alignments) == 1 and before >= 4:
        return []
    return new_alignments


def prepare_data(bed_df, cluster_mask, chromosome_lengths, threshold=500_000):
    # need to make sure rend > rstart for sortedintersect, and intervals are sorted
    bed_df['start'] = np.minimum(bed_df['rstart'], bed_df['rend'])
    bed_df['end'] = np.maximum(bed_df['rstart'], bed_df['rend'])
    bed_df.sort_values('start', inplace=True)
    columns = ['chrom', 'start', 'end', 'qname', 'n_alignments', 'qlen2']
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


# this function is redundant? seems to do the same job as calculate_overlap
# filter intervals that don't overlap at least as much as the user input value, return indexes
# def filter_overlap(overlaps, percentage, query_start, query_end):
#     filtered_overlaps = set([])
#     for i in overlaps:
#         if ((query_end - query_start) / (max(i[1], query_end) - min(i[0], query_start))) >= percentage:
#             filtered_overlaps.add(i[2])
#     return filtered_overlaps


def calculate_overlap(interval1, interval2):
    x_size = interval1.end - interval1.start
    y_size = interval2.end - interval2.start
    overlap = max(0, (min(interval1.end, interval2.end) - max(interval1.start, interval2.start)))
    reciprocal_overlap = min(overlap / x_size, overlap / y_size)
    return reciprocal_overlap


def overall_jaccard_similarity(l1, l2, percentage=0.8):
    if not l1 or not l2:
        return 0
    intersection = 0
    l1_comparisons = [0] * len(l1)
    l2_comparisons = [0] * len(l2)  # 0 for no intersection, 1 for an intersection
    for i, interval1 in enumerate(l1):
        for j, interval2 in enumerate(l2):
            if interval1.chrom == interval2.chrom and calculate_overlap(interval1, interval2) >= percentage:
                intersection += 1
                l1_comparisons[i] = 1
                l2_comparisons[j] = 1
                break
    intersections = max(sum(l1_comparisons), sum(l2_comparisons))
    union = intersections + l1_comparisons.count(0) + l2_comparisons.count(0)
    if not union:
        return 0
    return intersections / union, intersections


def get_chromosome_lengths(bam_path):
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    return {bam_file.get_reference_name(tid): l for tid, l in enumerate(bam_file.lengths) if l > 1000000}


def different_lengths_or_alignments(itv1, itv2, qlen_diff, diff):
    if (min(itv1.qlen2, itv2.qlen2) / max(itv1.qlen2, itv2.qlen2)) >= 1 - qlen_diff:
        return False
    if (min(itv1.n_alignments, itv2.n_alignments) / max(itv1.n_alignments, itv2.n_alignments)) >= 1 - diff:
        return False
    return True


@measure_time
def query_interval_trees(interval_trees, data, overlap_cutoff, jaccard_threshold, edge_threshold, qlen_diff, diff):
    query_intervals = defaultdict(list)
    for itv in data:
        query_intervals[itv.qname].append(itv)
    G = nx.Graph()
    seen_edges = set([])
    match = set([])

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

                list2 = query_intervals[o_data.qname]
                j, union = overall_jaccard_similarity(list1, list2, overlap_cutoff)
                target = jaccard_threshold[union] if union < len(jaccard_threshold) else jaccard_threshold[-1]
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


@measure_time
def choose_alignment(bed_file):
    # calculate average alignment scores for each read
    qname_grouped = bed_file.groupby('qname')
    avg_scores = qname_grouped['alignment_score'].mean()

    # Map the average alignment scores back to the bed_file using qname
    bed_file['avg_alignment_score'] = bed_file['qname'].map(avg_scores)

    # Group reads by cluster identifier
    cluster_grouped = bed_file.groupby('cluster')

    # list to store reads with highest score in each cluster
    selected_reads = []

    for cluster_id, group in cluster_grouped:
        # Find the read with the highest average alignment score in the cluster
        max_alignment_read = group.loc[group['avg_alignment_score'].idxmax()]['qname']
        selected_reads.append(max_alignment_read)

    selected_reads_df = bed_file[bed_file['qname'].isin(selected_reads)]

    return selected_reads_df
