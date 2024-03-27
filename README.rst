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
Normal usage:
::

    fslr --name samp1 \
         --out samp1_out \
         --ref T2T.fa \
         --primers 21q1 \
         --basecalled basecalled/samp1/pass \
         --procs 16

Skip clustering:
::
    fslr --name samp1 \
         --out samp1_out \
         --ref T2T.fa \
         --primers 21q1 \
         --basecalled basecalled/samp1/pass \
         --procs 16 \
         --skip-interval-cluster

Skip alignment:

    To skip the alignment step a [name].mappings.bed file is needed in the output folder.

::

    fslr --name samp1 \
         --out samp1_out \
         --ref T2T.fa \
         --primers 21q1 \
         --procs 16 \
         --skip-alignment

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

Outputs
-------
Out folder:

    * .without_primers.fa: Contains sequences of reads without identifiable primers.
    * .primers_labelled.fq: Contains sequences of uniquely labelled reads that have at least one identified primer.
    * .bwa_dodi.bam: A compressed binary file that contains the aligned sequences.
    * .bwa_dodi.bai: Index file.
    * .mappings.bed: A text file that stores genomic regions as coordinates associated with the split-reads.
    * .mappings.cluster.bed: The .mappings.bed file supplemented with information about the clusters.
    * .filter_counts_summary.csv: Contains summary information about the filtered reads and the clusters.

Out/cluster folder:

    * .cluster.specifications.csv: A text file listing the identified clusters and their attributes.
    * .cluster.consensus.fa: Consensus sequence of each cluster.
    * .cluster.without_primers.fa: Consensus sequences without primers.
    * .primers_labelled.fq: Uniquely labelled consensus sequences that have at least one identified primer.
    * .bwa_dodi_cluster.bam: A compressed binary file that contains the aligned consensus sequences.
    * .bwa_dodi_cluster.bai: Index file.


How it works
------------

1. Filter reads

    Remove repetitive sequences, junk sequences and concatemers

2. Find reads with primers

3. Align to the reference genome and choose the best alignments

4. Cluster the reads:

    It works by constructing a graph based on the level of overlapping intervals and utilizing Jaccard-similarity measures. This step is needed to determine whether the alignments at a specific position come from different fusion events or if they actually come from the same event sequenced multiple times due to one single fusion event being amplified.
