import pysam
from subprocess import run
import sys

def make_indexed_ref(mask_bed, outname, current_ref):
    current_ref = pysam.FastaFile(current_ref)
    print(f'Making masked reference genome {outname}', file=sys.stderr)
    with open(mask_bed) as f, open(outname, 'w') as ref:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            ref.write(f'>{line[0]}\n')
            start = int(line[1])
            end = int(line[2])
            if start > 0:
                ref.write('N'*start)
            ref.write(current_ref.fetch(line[0], start, end))
    run(f'bwa index {outname}', shell=True)
