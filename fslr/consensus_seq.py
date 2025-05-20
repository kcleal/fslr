import click
import pysam
import subprocess
import pandas as pd

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