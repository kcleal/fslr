import click
import pandas as pd
import os
import glob
from fslr import filter_junk_from_fq, find_reads_with_primers, collect_mapping_info, cluster_seqs
import multiprocessing
import subprocess
import sys

__version__ = '0.2'

@click.command()
@click.option('--name', required=True, help='Sample name')
@click.option('--out', required=True, help='Output folder')
@click.option('--ref', required=True, help='Reference genome')
@click.option('--basecalled', required=True, help='Folder of basecalled reads in fastq format to analyse')
@click.option('--primers', required=True, help='Comma-separated list of primer names. Make sure these are listed in primers.csv')
@click.option('--trim-threshold', required=False, help='Threshold in range 0-1. Fraction of maximum primer alignment score; primer sites with lower scores'
                                                      ' are labelled False',
              default=0.4, type=float, show_default=True)
@click.option('--cluster-fraction', required=False, help='The proportion of data to try and cluster, range 0-1',
              default='auto', show_default=True)
@click.option('--procs', required=False, default=1, show_default=True, help='Processors to use')
@click.option('--keep-temp', required=False, is_flag=True, flag_value=True, help='Keep temp files')
@click.option('--skip-alignment', required=False, is_flag=True, flag_value=True, help='Skip alignment step')
@click.option('--cluster', required=False, is_flag=True, flag_value=True, help='Skip clustering step')
@click.version_option(__version__)
def pipeline(**args):

    with multiprocessing.Manager() as manager:

        lock = manager.Lock()

        primers_d = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + '/primers.csv')
        args['primers'] = args['primers'].split(',')
        ps = set(primers_d['primer_name'])
        basename = f'{args["out"]}/{args["name"]}'
        print('Basename: ', basename, file=sys.stderr)

        for p in args['primers']:
            if p not in ps:
                raise ValueError('Input primer name not in primers.csv', p, ps)

        primers = {k: v for k, v in zip(primers_d['primer_name'], primers_d['primer_seq']) if k in args['primers']}
        primers_target = {k: v for k, v in zip(primers_d['primer_name'], primers_d['primer_alignment_target']) if k in args['primers']}

        filter_counts = manager.dict()
        filter_counts['name'] = args['name']
        filter_counts['total_kept'] = 0
        filter_counts['total_dropped'] = 0
        filter_counts['conactemers_dropped'] = 0
        filter_counts['junk_seqs_dropped'] = 0
        filter_counts['False_False'] = 0

        if not args['skip_alignment']:

            if not os.path.exists(args['out']):
                os.mkdir(args['out'])

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
                    pass
            else:
                for j in jobs:
                    filter_junk_from_fq.func(j)
                    pass

            jobs = []
            for pth in glob.glob(f'{args["out"]}/*filtered_junk.fq'):
                jobs.append((pth, primers_target, lock, filter_counts, args['keep_temp'], args['trim_threshold']))
            if args['procs'] > 1:
                with multiprocessing.Pool(args['procs']) as p:
                    p.map(find_reads_with_primers.func, jobs)
                    pass
            else:
                for j in jobs:
                    find_reads_with_primers.func(j)
                    pass

            print('Filter counts: ', filter_counts, file=sys.stderr)

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
                pass
                for pl in glob.glob(f"{basename}.*.primers_labelled.fq"):
                    os.remove(pl)

            assert len(glob.glob(f"{basename}.bwa_dodi.bam")) == 1

            collect_mapping_info.mapping_info(f"{basename}.bwa_dodi.bam", f"{basename}.mappings.bed")

        if args['cluster']:

            cf = args['cluster_fraction']
            if cf != 'auto':
                cf = float(cf)
            cluster_seqs.find_targets(filter_counts, "{basename}.bwa_dodi.bam".format(basename=basename), cluster_fraction=cf)
            cluster_seqs.cluster_paf(basename, args['procs'])

            os.remove("{basename}.minimap2_cluster.paf".format(basename=basename))

            c = "cat {basename}_clusters/.*.cons.fa | " \
                "bwa mem -c 1000 -A2 -B3 -O5 -E2 -T0 -L0 -D 0.25 -r 1.25 -d 200 -k 11 -a -t{procs} {ref} - |" \
                "dodi --paired False -c 1 -u 21 --ol-cost 2 --max-overlap 50000 - |" \
                "samtools view -bh - |" \
                "samtools sort -o {basename}.bwa_dodi.cons.bam; " \
                "samtools index {basename}.bwa_dodi.cons.bam".format(basename=basename, procs=args['procs'], ref=args['ref'])

            subprocess.run(c, shell=True)

            assert len(glob.glob(f"{basename}.bwa_dodi.cons.bam")) == 1

            collect_mapping_info.mapping_info(f"{basename}.bwa_dodi.cons.bam", f"{basename}.mappings.cons.bed")

        #
        with open(basename + '.filter_counts.csv', 'w') as fc:
            fc.write(','.join([str(k) for k in filter_counts.keys()]) + '\n')
            fc.write(','.join([str(k) for k in filter_counts.values()]) + '\n')

        print('fslr finished')
