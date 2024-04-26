import pandas as pd
from ncls import NCLS
import networkx as nx
import subprocess
import click
import pysam

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

# calculates jacc. sim. between 2 queries ("sandwiches")
# how many of the reads in the 2 queries overlap - intersection
# union - the number of fillings in the bigger query
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
    for e in mask:
        delete = read['chrom'].str.contains(e)
        read = read[~delete]
    if 'subtelomere' in mask:
        subtelomere = []
        for index, row in read.iterrows():
            chrom = row['chrom']
            if row['rstart'] > (chromosome_lengths[chrom] - 500000) or row['rend'] < 500000:
                subtelomere.add(index)
        read = read.drop(subtelomere)
    return read


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
    for query_key, query_values in query_dict.items():
        for chr_value, start_value, end_value, index in query_values:
            tree = interval_trees[chr_value]
            overlaps = tree.find_overlap(start_value, end_value)
            # filter for qlen2
            overlaps2 = set([])
            overlap_intervals = [o for o in overlaps]
            for o in overlap_intervals:
                qlen = int(bed.loc[o[2], 'qlen2'])
                if (int(bed.loc[index, 'qlen2']) * (1 + qlen_diff)) >= qlen and qlen >= (
                        int(bed.loc[index, 'qlen2']) * (1 - qlen_diff)):
                    overlaps2.add(o)
                else:
                    bad = bed.loc[o[2], 'qname']
                    bad_list.add(tuple(sorted((bad, query_key)))) # qname pairs where qlen was too different
                    seen_edges.add(tuple(sorted((bad, query_key))))

            filtered_overlap = filter_overlap(overlaps2, cutoff, start_value, end_value)

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
            query_name_of_overlaps = bed.loc[list(filtered_overlap2), 'qname']

            for q in query_name_of_overlaps:
                # check if q, query_key are already in seen or not
                if q != query_key and tuple(sorted((q, query_key))) not in seen_edges:
                    q1 = bed[bed['qname'] == query_key]
                    # mask TALENs and subtelomeres in case n_alignments > 3
                    if q1.iloc[0]['n_alignments'] > 3 and mask != None:
                        q1 = mask_sequences(q1, mask, chromosome_lengths)
                    set1 = set(q1[['chrom', 'rstart', 'rend']].itertuples(index=False, name=None))
                    q2 = bed[bed['qname'] == q]
                    if q2.iloc[0]['n_alignments'] > 3 and mask != None:
                        q2 = mask_sequences(q2, mask, chromosome_lengths)
                    set2 = set(q2[['chrom', 'rstart', 'rend']].itertuples(index=False, name=None))
                    j = overall_jaccard_similarity(set1, set2, cutoff)
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


def make_consensus_seq(subg, out, name, bed_file, primer_list):
    # make fasta files for each cluster -> use it to make a consensus sequence
    rows_list = []
    num = 0
    for clust in subg:
        seq_df = bed_file[bed_file['qname'].isin(clust)].dropna(subset=['seq'])[['qname', 'seq']]
        n_reads = len(clust)
        # Open the file for writing
        fasta_file_path = f'{out}/cluster/consensus_seq/{name}.cluster{num}.n_reads{n_reads}.fa'
        with open(fasta_file_path, 'w') as fasta_file:
            # Get the total number of rows in seq_df
            total_rows = seq_df.shape[0]
            # Iterate through DataFrame rows and write to the FASTA file
            for index, row in seq_df.iterrows():
                qname = row['qname']
                seq = row['seq']
                # Write sequence header
                fasta_file.write(f'>{qname}\n')
                # Write sequence data without trailing newline if it's the last row
                if index == total_rows - 1:
                    fasta_file.write(f'{seq}')
                else:
                    # Write sequence data with newline for non-last rows
                    fasta_file.write(f'{seq}\n')

        # create consensus sequence
        a = (f'abpoa {out}/cluster/consensus_seq/{name}.cluster{num}.n_reads{n_reads}.fa ' \
             f'| sed "s/Consensus_sequence/cluster:{num}.n_reads:{n_reads}/g" ' \
             f'> {out}/cluster/consensus_seq/{name}.cluster{num}.n_reads{n_reads}.cons.fa')
        # write stdout to a log file instead of the screen
        with open(f"{out}/cluster/abpoa_logfile.txt", "a") as log_file:
            subprocess.run(a, shell=True, stdout=log_file, stderr=subprocess.STDOUT)

        # add consensus sequence to the df
        with open(f'{out}/cluster/consensus_seq/{name}.cluster{num}.n_reads{n_reads}.cons.fa', 'r') as f:
            sequence = ""
            for line in f:
                if not line.startswith(">"):
                    sequence += line.strip()

        # df about the clusters
        df = bed_file[bed_file['qname'].isin(clust)]


        # add purity of the cluster
        row = [num, n_reads]
        for p in primer_list:
            row.append((df['qname'].str.contains(p).sum()/df.shape[0]))

        row.append(sequence)
        rows_list.append(row)
        num += 1

    col = ['cluster', 'n_reads']
    col.extend(primer_list)
    col.append('cons_seq')
    cluster_df = pd.DataFrame(rows_list, columns=col)
    cluster_df.to_csv(f'{out}/cluster/{name}.cluster.purity.csv', index= False)

def delete_alignments(input_bam, output_bam, alignments_to_delete):
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        header = infile.header
        with pysam.AlignmentFile(output_bam, "wb", header=header) as outfile:
            for alignment in infile:
                # Check if the alignment's query name is in the list of names to delete
                if alignment.query_name not in alignments_to_delete:
                    outfile.write(alignment)

def merge_bam_files(input_bam1, input_bam2, output_bam):
    with pysam.AlignmentFile(output_bam, "wb", header=pysam.AlignmentFile(input_bam1, "rb").header) as outfile:
        with pysam.AlignmentFile(input_bam1, "rb") as infile1, pysam.AlignmentFile(input_bam2, "rb") as infile2:
            for alignment in infile1:
                outfile.write(alignment)
            for alignment in infile2:
                outfile.write(alignment)


