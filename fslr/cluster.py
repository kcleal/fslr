import pandas as pd
from ncls import NCLS
import networkx as nx
import subprocess
import click
import pysam
import time
from functools import wraps
from sortedintersect import IntervalSet
from collections import defaultdict, namedtuple


IntervalItem = namedtuple('interval_item', ['chrom', 'rstart', 'rend', 'index'])
QueryItem = namedtuple('query_item', ['qlen2', 'n_alignments'])
IndexItem = namedtuple('index_item', ['qname', 'n_alignments', 'qlen2'])


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


def build_interval_trees(bed_df):
    interval_tree = defaultdict(lambda: IntervalSet(with_data=True))
    b = bed_df.sort_values(['chrom', 'rstart'])
    for chrom, start, end, index in zip(b['chrom'], b['rstart'], b['rend'], b.index):
        interval_tree[chrom].add(min(start, end), max(start, end), index)
    return interval_tree


# filter intervals that don't overlap at least as much as the user input value, return indexes
def filter_overlap(overlaps, percentage, query_start, query_end):
    filtered_overlaps = set([])
    for ii in overlaps:
        if ((query_end - query_start) / (max(int(ii[1]), query_end) - min(int(ii[0]), query_start))) >= percentage:
            filtered_overlaps.add(int(ii[2]))
    return filtered_overlaps


def calculate_overlap(interval1, interval2):
    x1, x2 = interval1[1], interval1[2]
    x_size = x2 - x1
    y1, y2 = interval2[1], interval2[2]
    y_size = y2 - y1
    overlap = max(0, (min(x2, y2) - max(x1, y1)))
    reciprocal_overlap = min(overlap / x_size, overlap / y_size)
    return reciprocal_overlap


def overall_jaccard_similarity(list1, list2, percentage):
    # Calculate overall intersection and union
    intersection = 0
    union = max(len(list1), len(list2))
    if union == 0:
        overall_jaccard_similarity = 0
    else:
        for interval1 in list1:
            for interval2 in list2:
                # Check if intervals are on the same chromosome
                if interval1[0] == interval2[0]:
                    if calculate_overlap(interval1, interval2) >= percentage:
                        intersection += 1
        overall_jaccard_similarity = intersection / union
    return overall_jaccard_similarity


def get_chromosome_lengths(bam_path):
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    return {bam_file.get_reference_name(tid): l for tid, l in enumerate(bam_file.lengths) if l > 1000000}


def mask_sequences2(read_alignments, mask, chromosome_lengths, threshold=500_000):
    if not mask:
        return read_alignments
    new_alignments = set([])
    before = len(read_alignments)
    for a in read_alignments:
        if a.chrom in mask:
            continue
        if 'subtelomere' in mask:
            if a.rstart < threshold or (a.chrom in chromosome_lengths and chromosome_lengths[a.chrom] - a.rend < threshold):
                continue
        new_alignments.add(a)
    if len(read_alignments) == 1 and before >= 4:
        return set([])
    return new_alignments


@measure_time
def query_interval_trees(interval_trees, bed, chromosome_lengths, cutoff, jaccard_threshold, edge_threshold, qlen_diff,
                         diff, mask):
    # prepare data tables
    query_intervals = {qname: set(IntervalItem(*i) for i in zip(g['chrom'], g['rstart'], g['rend'], g.index))
                       for qname, g in bed.groupby('qname')}
    query_intervals = {k: mask_sequences2(v, mask, chromosome_lengths) for k, v in query_intervals.items()}
    query_intervals = {k: v for k, v in query_intervals.items() if v}
    query_data = {qname: QueryItem(d.iloc[0]['qlen2'], d.iloc[0]['n_alignments']) for qname, d in bed.groupby('qname')}

    index_data = {index: IndexItem(name, n_aligns, qlen2) for index, name, n_aligns, qlen2 in
                  zip(bed.index, bed['qname'], bed['n_alignments'], bed['qlen2'])}

    G = nx.Graph()
    seen_edges = set([])
    match = set([])
    bad_list = set([])  # too big qlen difference, too big alignment number difference, too low jaccard similarity

    for query_key, set1 in query_intervals.items():
        edges = 0
        current_data = query_data[query_key]

        for itv in set1:
            overlap_intervals = interval_trees[itv.chrom].search_interval(min(itv.rstart, itv.rend), max(itv.rstart, itv.rend))
            overlaps2 = set([])

            for o in overlap_intervals:
                qlen2 = index_data[o[2]].qlen2
                if (min(current_data.qlen2, qlen2) / max(current_data.qlen2, qlen2)) >= 1 - qlen_diff:
                    overlaps2.add(o)
                else:
                    bad = index_data[o[2]].qname
                    bad_list.add(tuple(sorted((bad, query_key))))  # qname pairs where qlen2 was too different
                    seen_edges.add(tuple(sorted((bad, query_key))))

            # only keep intervals that overlap above a specified threshold
            filtered_overlap_indexes = filter_overlap(overlaps2, cutoff, itv.rstart, itv.rend)

            # check the difference between the number of alignments
            filtered_overlap2 = set([])
            for idx in filtered_overlap_indexes:
                alignment_num = index_data[idx].n_alignments
                if (min(current_data.n_alignments, alignment_num) / max(current_data.n_alignments, alignment_num)) >= 1 - diff:
                    filtered_overlap2.add(idx)
                else:
                    bad = index_data[idx].qname
                    bad_list.add(tuple(sorted((bad, query_key))))
                    seen_edges.add(tuple(sorted((bad, query_key))))

            # retrieve the query names of reads that passed
            for q in [index_data[i].qname for i in filtered_overlap2]:
                if q == query_key or tuple(sorted((q, query_key))) in seen_edges:
                    continue
                set2 = query_intervals[q]
                if not len(set2):
                    continue
                j = overall_jaccard_similarity(set1, set2, cutoff)
                if j >= jaccard_threshold:
                    match.add((query_key, q, j))
                    G.add_edge(query_key, q)
                    edges += 1
                else:
                    bad_list.add((q, query_key))
                seen_edges.add(tuple(sorted((q, query_key))))
                if edges >= edge_threshold:
                    break

    match_df = pd.DataFrame(match, columns=['query1', 'query2', 'jaccard_similarity'])
    return match_df, G, bad_list


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
