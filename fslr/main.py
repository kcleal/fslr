import click
import pandas as pd
import os
import glob
from fslr import filter_junk_from_fq, find_reads_with_primers, collect_mapping_info
import multiprocessing
import subprocess
import sys


@click.command()
@click.option('--name', required=True, help='Sample name')
@click.option('--out', required=True, help='Output folder')
@click.option('--ref', required=True, help='Reference genome')
@click.option('--basecalled', required=True, help='Folder of basecalled reads in fastq format to analyse')
@click.option('--primers', required=True, help='Comma-separated list of primer names. Make sure these are listed in primers.csv')
@click.option('--procs', required=False, default=8, help='Processors to use')
@click.option('--keep-temp', required=False, is_flag=True, flag_value=True, help='Keep temp files')
def pipeline(**args):

    with multiprocessing.Manager() as manager:

        lock = manager.Lock()

        primers_d = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + '/primers.csv')
        args['primers'] = args['primers'].split(',')
        ps = set(primers_d['primer_name'])
        basename = f'{args["out"]}/{args["name"]}'

        for p in args['primers']:
            if p not in ps:
                raise ValueError('Input primer name not in primers.csv', p, ps)

        primers = {k: v for k, v in zip(primers_d['primer_name'], primers_d['primer_seq']) if k in args['primers']}
        primers_target = {k: v for k, v in zip(primers_d['primer_name'], primers_d['primer_alignment_target']) if k in args['primers']}

        if not os.path.exists(args['out']):
            os.mkdir(args['out'])

        filter_counts = manager.dict()
        filter_counts['name'] = args['name']
        filter_counts['total_kept'] = 0
        filter_counts['total_dropped'] = 0
        filter_counts['conactemers_dropped'] = 0
        filter_counts['junk_seqs_dropped'] = 0

        v = []
        for k in list(primers.keys()) + ['False']:
            for k2 in list(primers.keys()) + ['False']:
                if k == 'False':
                    p1 = 'False'
                else:
                    p1 = k + 'F'
                if k2 == 'False':
                    p2 = 'False'
                else:
                    p2 = k2 + 'R'
                v.append("_".join((p1, p2)))
                v.append("_".join((p2, p1)))

        for i in sorted(v):
            filter_counts[i] = 0

        print('Filtering reads: ', args['basecalled'], file=sys.stderr)
        fs = glob.glob(f'{args["basecalled"]}/*.fq.gz') + glob.glob(f'{args["basecalled"]}/*.fq') + glob.glob(
                             f'{args["basecalled"]}/*.fastq.gz') + glob.glob(f'{args["basecalled"]}/*.fastq') + glob.glob(f'{args["basecalled"]}/*.fasta') \
            + glob.glob(f'{args["basecalled"]}/*.fa') + glob.glob(f'{args["basecalled"]}/*.fasta.gz') + glob.glob(f'{args["basecalled"]}/*.fa.gz')
        print('Input files:', fs, file=sys.stderr)

        jobs = []
        for pth in fs:
            jobs.append((pth, primers, args['out'], args['name'], filter_counts, lock, args['keep_temp']))
        if args['procs'] > 1:
            with multiprocessing.Pool(args['procs']) as p:
                p.map(filter_junk_from_fq.func, jobs)
        else:
            for j in jobs:
                filter_junk_from_fq.func(j)

        jobs = []
        for pth in glob.glob(f'{args["out"]}/*filtered_junk.fq'):
            jobs.append((pth, primers_target, lock, filter_counts, args['keep_temp']))
        if args['procs'] > 1:
            with multiprocessing.Pool(args['procs']) as p:
                p.map(find_reads_with_primers.func, jobs)
        else:
            for j in jobs:
                find_reads_with_primers.func(j)

        print('Filter counts: ', filter_counts, file=sys.stderr)
        with open(basename + '.filter_counts.csv', 'w') as fc:
            fc.write(','.join([str(k) for k in filter_counts.keys()]) + '\n')
            fc.write(','.join([str(k) for k in filter_counts.values()]) + '\n')

        subprocess.run(f"cat {args['out']}/*.no_primers.fq > {basename}.without_primers.fq", shell=True)
        subprocess.run(f"rm {args['out']}/*.no_primers.fq", shell=True)

        c = "cat {basename}.*.primers_labelled.fq | " \
            "bwa mem -c 1000 -A2 -B3 -O5 -E2 -T0 -L0 -D 0.25 -r 1.25 -d 200 -k 11 -a -t{procs} {ref} - |" \
            "dodi --paired False -c 1 -u 21 --ol-cost 2 --max-overlap 50000 - |" \
            "samtools view -bh - |" \
            "samtools sort -o {basename}.bwa_dodi.bam; " \
            "samtools index {basename}.bwa_dodi.bam".format(basename=basename, procs=args['procs'], ref=args['ref'])

        subprocess.run(c, shell=True)

        if not args['keep_temp']:
            for pl in glob.glob(f"{basename}.*.primers_labelled.fq"):
                os.remove(pl)

        assert len(glob.glob(f"{basename}.bwa_dodi.bam")) == 1

        collect_mapping_info.mapping_info(f"{basename}.bwa_dodi.bam", f"{basename}.mappings.bed")

        print('fslr finished')
