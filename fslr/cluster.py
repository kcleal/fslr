import pandas as pd
from ncls import NCLS
import networkx as nx
import subprocess
import click
import pysam
import time
from functools import wraps


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
    # transform df to dictionary that has chromosomes as keys
    interval_dict = {}
    for index, row in bed_df[['chrom', 'rstart', 'rend']].iterrows():
        key = row['chrom']
        if key not in interval_dict:
            interval_dict[key] = []
        interval_dict[key].append([row['rstart'], row['rend'], index])
    # make an interval tree for each chromosome
    interval_trees = {}
    for chrom, chrom_intervals in interval_dict.items():
        sorted_intervals = sorted(chrom_intervals, key=lambda x: (x[0], x[1], x[2]))
        starts = [interval[0] for interval in sorted_intervals]
        ends = [interval[1] for interval in sorted_intervals]
        idx = [interval[2] for interval in sorted_intervals]
        tree = NCLS(starts, ends, ids=idx)
        interval_trees[chrom] = tree
    return interval_trees

# filter intervals that don't overlap at least as much as the user input value
def filter_overlap(overlaps, percentage, query_start, query_end):
    filtered_overlaps = set([])
    interval = [i for i in overlaps]
    for ii in interval:
        if ((query_end - query_start) / (max(int(ii[1]), query_end) - min(int(ii[0]), query_start))) >= percentage:
            filtered_overlaps.add(int(ii[2])) #this will return the indexes of the overlaps
    return filtered_overlaps

#needed for calculating Jaccard similarity
def calculate_overlap(interval1, interval2):
    len_1 = interval1[2] - interval1[1]
    len_2 = max(interval1[2], interval2[2]) - min(interval1[1], interval2[1])
    percent = len_1 / len_2
    return percent

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

    # Calculate overall Jaccard similarity
        overall_jaccard_similarity = intersection / union

    return overall_jaccard_similarity

# get chromosome length from the header of the alignment file -> find subtelomere regions
def get_chromosome_lengths(bam_file):
    chromosome_lengths = {}
    command = ["samtools", "view", "-H", bam_file]
    header = subprocess.check_output(command, universal_newlines=True)
    for line in header.split('\n'):
        if line.startswith('@SQ'):
            fields = line.split('\t')
            for field in fields:
                if field.startswith('SN:'):
                    chrom = field.split(':')[1]
                elif field.startswith('LN:'):
                    length = int(field.split(':')[1])
            chromosome_lengths[chrom] = length
    chromosome_lengths = {key: value for key, value in chromosome_lengths.items() if value > 1000000}
    return chromosome_lengths

def mask_sequences(read, mask, chromosome_lengths):

    mask_set = set(mask)

    # Filter out rows based on mask
    delete = read['chrom'].str.contains('|'.join(mask_set))
    read = read[~delete]

    # Check for subtelomere in mask
    if 'subtelomere' in mask_set:
        subtelomere = (
                (read['chrom'].isin(chromosome_lengths.keys())) &
                ((read['rstart'] > (read['chrom'].map(chromosome_lengths) - 500000)) |
                 (read['rend'] < 500000))
        )
        read = read[~subtelomere]
    return read

@measure_time
def query_interval_trees(interval_trees, bed, chromosome_lengths, cutoff, jaccard_threshold, edge_threshold, qlen_diff, diff, mask):
    # transform df to dictionary that has the qnames as keys
    query_dict = {}
    for index, row in bed[['qname', 'chrom', 'rstart', 'rend']].iterrows():
        key = row['qname']
        if key not in query_dict:
            query_dict[key] = []
        query_dict[key].append((row['chrom'], row['rstart'], row['rend'], index))
    # Create a graph
    G = nx.Graph()
    seen_edges = set([])  # once I found an edge I don't want to loop through those qname pairs again
    match = set([])  # these will be added to the graph
    edges = 0
    bad_list = set([]) # too big qlen difference, too big alignment number difference, too low jaccard similarity
    # loop through the reads
    for query_key, query_values in query_dict.items():
        # loop through the alignments within the reads
        for chr_value, start_value, end_value, index in query_values:
            tree = interval_trees[chr_value] # get the interval tree
            overlaps = tree.find_overlap(start_value, end_value) # find overlapping intervals
            # filter for qlen2 (qlen2 = the length of the fillings)
            overlaps2 = set([])
            overlap_intervals = [o for o in overlaps]
            for o in overlap_intervals:
                qlen = int(bed.loc[o[2], 'qlen2'])
                if (int(bed.loc[index, 'qlen2']) * (1 + qlen_diff)) >= qlen and qlen >= (
                        int(bed.loc[index, 'qlen2']) * (1 - qlen_diff)):
                    overlaps2.add(o)
                else:
                    bad = bed.loc[o[2], 'qname']
                    bad_list.add(tuple(sorted((bad, query_key)))) # qname pairs where qlen2 was too different
                    seen_edges.add(tuple(sorted((bad, query_key))))
            # only keep intervals that overlap above a specified threshold
            filtered_overlap = filter_overlap(overlaps2, cutoff, start_value, end_value)
            # check the difference between the number of alignments
            filtered_overlap2 = set([])
            for aln in filtered_overlap:
                alignment_num = bed.loc[aln, 'n_alignments']
                if (int(bed.loc[index, 'n_alignments']) * (1 + diff)) >= alignment_num and alignment_num >= (
                        int(bed.loc[index, 'n_alignments']) * (1 - diff)):
                    filtered_overlap2.add(aln)
                else:
                    bad = bed.loc[aln, 'qname']
                    bad_list.add(tuple(sorted((bad, query_key))))
                    seen_edges.add(tuple(sorted((bad, query_key))))
            # retrieve the query names of reads that passed
            query_name_of_overlaps = bed.loc[list(filtered_overlap2), 'qname']
            # loop through these qnames
            for q in query_name_of_overlaps:
                # check if q, query_key are already in seen or not
                if q != query_key and tuple(sorted((q, query_key))) not in seen_edges:
                    # mask regions in case n_alignments > 3
                    if bed[bed['qname'] == query_key].iloc[0]['n_alignments'] > 3 and mask != None:
                        q1 = mask_sequences(bed[bed['qname'] == query_key], mask, chromosome_lengths)
                    else:
                        q1 = bed[bed['qname'] == query_key]
                    set1 = set(q1[['chrom', 'rstart', 'rend']].itertuples(index=False, name=None))
                    if bed[bed['qname'] == query_key].iloc[0]['n_alignments'] > 3 and mask != None:
                        q2 = mask_sequences(bed[bed['qname'] == query_key], mask, chromosome_lengths)
                    else:
                        q2 = bed[bed['qname'] == query_key]
                    set2 = set(q2[['chrom', 'rstart', 'rend']].itertuples(index=False, name=None))
                    # calculate jaccard sim. between each read
                    j = overall_jaccard_similarity(set1, set2, cutoff)
                    # check if jaccard sim. is above the threshold and add to the graph if it is
                    if j >= jaccard_threshold:
                        match.add((query_key, q, j))
                        G.add_node(q)
                        if edges == 0:
                            G.add_node(query_key)
                        G.add_edge(query_key, q)
                        edges += 1
                    else:
                        bad_list.add((q, query_key))
                    seen_edges.add(tuple(sorted(((q, query_key)))))

                if edges >= edge_threshold:
                    break  # breaks the inner loop for the current query

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






