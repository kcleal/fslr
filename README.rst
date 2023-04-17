====
fslr
====

Identify split-read alignments from a Fusion-Seq long read experiment.


Installation
------------
Dependencies::

    bwa
    tantan
    minimap2
    samtools
    dodi
    abpoa

Install using::

    pip install -r requirements.txt
    pip install .


Usage
-----

::

    fslr --name samp1 \
         --out samp1_out \
         --ref T2T.fa \
         --primers 21q1 \
         --basecalled basecalled/samp1/pass \
         --procs 16
