====
fslr
====

Identify split-read alignments from a Fusion-Seq long read experiment.

Table of contents
-----------------

- `Installation`_
- `Usage`_
- `Options`_
- `Outputs`_
- `How it works`_

Installation
------------
Dependencies::

    bwa
    tantan
    samtools
    dodi
    abpoa

Install using::

    pip install -r requirements.txt
    pip install .

Usage
-----
Default usage:
::

    fslr --name samp1 \
         --out samp1_out \
         --ref T2T.fa \
         --primers 21q1,17p6 \
         --basecalled basecalled/samp1/pass \
         --procs 16
         --mask subtelomere,TALEN

Skip clustering:
::
    fslr --name samp1 \
         --out samp1_out \
         --ref T2T.fa \
         --primers 21q1,17p6 \
         --basecalled basecalled/samp1/pass \
         --procs 16 \
         --skip-clustering

Skip alignment:

    To skip the alignment step a [name].mappings.bed, [name].bwa_dodi.bam file is needed in the output folder.

::

    fslr --name samp1 \
         --out samp1_out \
         --ref T2T.fa \
         --primers 21q1 \
         --procs 16 \
         --skip-alignment \
         --mask subtelomere,TALEN

Options
-------

+---------------------------+------------------------------------------------------------------------------------------+
| Option                    | Description                                                                              |
+===========================+==========================================================================================+
| `--name`                  | Sample name.                                                                             |
+---------------------------+------------------------------------------------------------------------------------------+
| `--out`                   | Output folder.                                                                           |
+---------------------------+------------------------------------------------------------------------------------------+
| `--ref`                   | Reference genome.                                                                        |
+---------------------------+------------------------------------------------------------------------------------------+
| `--basecalled`            | Folder of basecalled reads in fastq format to analyse.                                   |
+---------------------------+------------------------------------------------------------------------------------------+
| `--primers`               | Comma-separated list of primer names. Make sure these are listed in primers.csv.         |
+---------------------------+------------------------------------------------------------------------------------------+
| `--trim-threshold`        | Threshold in range 0-1. Fraction of maximum primer alignment score; primer sites with    |
|                           | lower scores are labelled False.                                                         |
+---------------------------+------------------------------------------------------------------------------------------+
| `--keep-temp`             | Keep temp files.                                                                         |
+---------------------------+------------------------------------------------------------------------------------------+
| `--regions`               | Target regions in bed form to perform biased mapping.                                    |
+---------------------------+------------------------------------------------------------------------------------------+
| `--bias`                  | Multiply alignment score by bias if alignment falls within target regions.               |
+---------------------------+------------------------------------------------------------------------------------------+
| `--procs`                 | Number of processors to use.                                                             |
+---------------------------+------------------------------------------------------------------------------------------+
| `--skip-alignment`        | Skip alignment step.                                                                     |
+---------------------------+------------------------------------------------------------------------------------------+
| `--skip-interval-cluster` | Skip clustering step.                                                                    |
+---------------------------+------------------------------------------------------------------------------------------+
| `--jaccard-cutoff`        | Jaccard similarity index, a number between 0-1, below which reads won't be considered in |
|                           | the same cluster.                                                                        |
+---------------------------+------------------------------------------------------------------------------------------+
| `--overlap`               | A number between 0 and 1. Zero means two reads don't overlap at all, while 1 means the   |
|                           | start and end of the reads is identical.                                                 |
+---------------------------+------------------------------------------------------------------------------------------+
| `--n-alignmentdiff`       | How much the number of alignments in one cluster can differ. Fraction in the range 0-1.  |
+---------------------------+------------------------------------------------------------------------------------------+
| `--qlen-diff`             | Max difference in query length. Fraction in the range 0-1.                               |
+---------------------------+------------------------------------------------------------------------------------------+
| `--mask`                  | Comma separated list of regions/chromosomes to be excluded from the clustering e.g.:     |
|                           | subtemoleric regions, TALEN.                                                             |
+---------------------------+------------------------------------------------------------------------------------------+

Outputs
-------
Default usage
=============

Out folder:

* .without_primers.fq: Contains sequences of reads without identifiable primers.
* .mappings.bed: A text file that stores genomic regions as coordinates associated with the split-reads.
* .mappings.cluster.bed: Contains the same information about the reads as .mappings.bed with two additional columns; cluster and n_reads. The cluster column stores the cluster id-s of the reads. The n_reads column shows the number of reads within a cluster.
* .mappings_merged.bed: This file contains genomic regions of all the "singletons" from the initial alignment and the re-aligned consensus sequences.
* .bwa_dodi.bam: Alignment file after the initial alignment step.
* .bwa_dodi_cluster_merged.bam: Alignment file containing the "singletons" and the consensus sequences.
* .bai: Index files.
* .filter_counts_summary.csv: Contains information about the filtered reads.

Out/cluster folder:

* .cluster.consensus.fa: Consensus sequences of each cluster.
* .cluster.without_primers.fq: Consensus sequences without an identified primer.
* abpoa_logfile.txt: Messages (standard output) created by abPOA while generating the consensus sequences.
* .cluster.purity.csv: List the cluster id-s, the number of reads within a cluster, the consensus sequences and the proportion of reads within a cluster that have a specific primer.

Skip clustering
===============

Out folder:

* .without_primers.fa: Contains sequences of reads without identifiable primers.
* .bwa_dodi.bam: A compressed binary file that contains the aligned reads.
* .bwa_dodi.bai: Index file.
* .mappings.bed: A text file that stores genomic regions as coordinates associated with the split-reads.



How it works
------------

1. Filter reads:

    Remove repetitive sequences, junk sequences and concatemers from the input files.

2. Find reads with primers:

    Identify primers at the end of the reads and exclude any read from further analysis that doesn't have at least one
    primer at one end. The result of 1. and 2. is summarised in [name].filter_counts_summary.csv.

3. Align to the reference genome and choose the best alignments:

    Reads are aligned to the user specified reference genome using bwa mem. Out of the possible alignments the best are
    then selected using dodi.
    A BAM and BED file is saved at this stage; [name].bwa_dodi.bam, [name].mappings.bed.

4. Cluster the reads:

    The purpose of the clustering step is to identify highly similar reads that are potentially the result of the same
    event getting amplified prior to the sequencing.
    It works by constructing a graph based on the level of overlapping intervals and utilizing Jaccard-similarity
    measures.
    A [name].mappings.cluster.bed file is created that shows which reads and alignments are in the same cluster.

5. Choose a representative read for each cluster:

    Calculate an average alignment score for each read. Choose a representative read in each cluster with the highest
    average alignment score.



Reference masking and biased alignments
---------------------------------------

A bed file can be added using the --reference-mask option. These regions will be used to
create a masked reference where regions outside bed regions will be converted to N's.

Reads will be mapped to the masked reference, in addition to the main reference. All mappings will then
be subsequently processed by dodi do obtain the final set of alignments per read.

Note, currently only one region per chromosome is supported when using --reference-mask option.

Biased mapping is also supported by supplying a bed file using the --regions option. Alignments that overlap a target
region will have a bias added to their alignment scores during processing with dodi, making them more likely
to be chosen as a final output alignment. Output alignments retain their original alignment scores.