import click
import pandas as pd
import glob
from fslr import filter_junk_from_fq, find_reads_with_primers, collect_mapping_info, cluster
import multiprocessing
import subprocess
import sys
import os
import pkg_resources

__version__ = pkg_resources.require("fslr")[0].version
file_path = os.path.dirname(os.path.realpath(__file__))


@click.command()
@click.option('--name', required=True, help='Sample name')
@click.option('--out', required=True, help='Output folder')
@click.option('--ref', required=True, help='Reference genome')
@click.option('--basecalled', required=False, help='Folder of basecalled reads in fastq format to analyse')
@click.option('--primers', required=True, help='Comma-separated list of primer names. Make sure these are listed in primers.csv')
@click.option('--trim-threshold', required=False, help='Threshold in range 0-1. Fraction of maximum primer alignment score; primer sites with lower scores'
                                                      ' are labelled False',
              default=0.4, type=float, show_default=True)
@click.option('--keep-temp', required=False, is_flag=True, flag_value=True, help='Keep temp files')
@click.option('--regions', required=False, type=click.Path(exists=True), help='Target regions in bed form to perform biased mapping')
@click.option('--bias', required=False, default=1.05, show_default=True, type=float, help='Multiply alignment score by bias if alignment falls within target regions')
@click.option('--procs', required=False, default=1, show_default=True, help='Number of processors to use')
@click.option('--skip-alignment', required=False, is_flag=True, help='Skip alignment step')
@click.option('--skip-interval-cluster', required=False, is_flag=True, help='Skip clustering step')
@click.option('--jaccard-cutoff', required=False, default=0.7, show_default=True, help="Jaccard similarity index, a number between 0-1, below which reads won't be considered in the same cluster")
@click.option('--overlap', required=False, default=0.8, show_default=True, help="A number between 0 and 1. Zero means two reads don't overlap at all, while 1 means the start and end of the reads is identical.")
@click.option('--n-alignmentdiff', default=0.25, required=False, show_default=True, help='How much the number of alignments in one cluster can differ. Fraction in the range 0-1.')
@click.option('--qlen-diff', default=0.04, required=False, show_default=True, help="Max difference in query length. Fraction in the range 0-1.")
@click.version_option(__version__)

def pipeline(**args):

    with multiprocessing.Manager() as manager:

        lock = manager.Lock()

        basename = f'{args["out"]}/{args["name"]}'
        print('Basename: ', basename, file=sys.stderr)

        filter_counts = manager.dict()
        filter_counts['name'] = args['name']
        filter_counts['total_kept'] = 0
        filter_counts['total_dropped'] = 0
        filter_counts['conactemers_dropped'] = 0
        filter_counts['junk_seqs_dropped'] = 0
        filter_counts['False_False'] = 0

        primers_d = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + '/primers.csv')
        args['primers'] = args['primers'].split(',')
        ps = set(primers_d['primer_name'])

        bias_params = "" if not args['regions'] else f"--bias {args['bias']} --include {args['regions']}"

        for p in args['primers']:
            if p not in ps:
                raise ValueError('Input primer name not in primers.csv', p, ps)

        primers = {k: v for k, v in zip(primers_d['primer_name'], primers_d['primer_seq']) if k in args['primers']}
        primers_target = {k: v for k, v in zip(primers_d['primer_name'], primers_d['primer_alignment_target']) if
                          k in args['primers']}

        if not os.path.exists(args['out']):
            os.mkdir(args['out'])

        if not args['skip_alignment']:

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
                jobs.append((pth, primers_target, lock, filter_counts, args['keep_temp'], args['trim_threshold']))
            if args['procs'] > 1:
                with multiprocessing.Pool(args['procs']) as p:
                    p.map(find_reads_with_primers.func, jobs)

            else:
                for j in jobs:
                    find_reads_with_primers.func(j)


            print('Filter counts: ', filter_counts, file=sys.stderr)

            subprocess.run(f"cat {args['out']}/*.no_primers.fq > {basename}.without_primers.fq", shell=True)
            subprocess.run(f"rm {args['out']}/*.no_primers.fq", shell=True)

            c = "cat {basename}.*.primers_labelled.fq | " \
                "bwa mem -c 1000 -A2 -B3 -O5 -E2 -T0 -L0 -D 0.25 -r 1.25 -d 200 -k 11 -a -t{procs} {ref} - |" \
                "dodi {bias_params} --paired False -c 1 -u 21 --ol-cost 2 --max-overlap 50000 - |" \
                "samtools view -bh - |" \
                "samtools sort -o {basename}.bwa_dodi.bam; " \
                "samtools index {basename}.bwa_dodi.bam".format(bias_params=bias_params,
                                                                basename=basename,
                                                                procs=args['procs'],
                                                                ref=args['ref'])

            subprocess.run(c, shell=True)

            if not args['keep_temp']:
                pass
                for pl in glob.glob(f"{basename}.*.primers_labelled.fq"):
                    os.remove(pl)

            assert len(glob.glob(f"{basename}.bwa_dodi.bam")) == 1

            collect_mapping_info.mapping_info(f"{basename}.bwa_dodi.bam",
                                              f"{basename}.mappings.bed",
                                              args['regions'])

            with open(f'{basename}.filter_counts_summary.csv', 'w') as fc:
                fc.write('Filter counts:' + '\n')
                fc.write(','.join([str(k) for k in filter_counts.keys()]) + '\n')
                fc.write(','.join([str(k) for k in filter_counts.values()]) + '\n')

        if not args['skip_interval_cluster']:

            filter_counts = manager.dict()
            filter_counts['name'] = f'{args["name"]}_cluster'
            filter_counts['total_kept'] = 0
            filter_counts['total_dropped'] = 0
            filter_counts['conactemers_dropped'] = 0
            filter_counts['junk_seqs_dropped'] = 0
            filter_counts['False_False'] = 0

            if not os.path.exists(args['out']):
                os.mkdir(args['out'])
            if not os.path.exists(f'{args["out"]}/cluster'):
                os.mkdir(f'{args["out"]}/cluster')

            print('Making clusters')
            # read in bed file
            bed_file = pd.read_csv(f'{basename}.mappings.bed', sep='\t')

            # arguments
            jaccard_cutoff = args['jaccard_cutoff']
            overlap = args['overlap']
            edge_threshold = 3
            qlen_diff = args['qlen_diff']
            n_alignments_diff = args['n_alignmentdiff']

            # delete the "breads", make qlen2 column == qlen without the breads
            filtered = cluster.keep_fillings(bed_file)
            # build interval trees for each chr
            interval_tree = cluster.build_interval_trees(filtered)
            # find queries that are similar and add them to a graph
            match_data, network, no_match = cluster.query_interval_trees(interval_tree, filtered, overlap, jaccard_cutoff, edge_threshold, qlen_diff, n_alignments_diff)
            # extract the subgraphs from the network
            subgraphs = cluster.get_subgraphs(network)
            if len(list(subgraphs)) == network.number_of_nodes():
                print("No clusters were found.")
                return
            subg_df = pd.DataFrame(subgraphs)
            subg_df = subg_df.T
            # subg_df.to_csv(f'{clustername}.cluster.qnames.csv', index=False)
            subg_long = pd.melt(subg_df, var_name='clusters', value_name='qname').dropna()
            subg_long['clusters'] = pd.to_numeric(subg_long['clusters'], errors='coerce')
            regions = bed_file[bed_file['qname'].isin(subg_long['qname'])]
            regions = regions[['chrom', 'rstart', 'rend']]
            bed_file = bed_file.merge(subg_long, on='qname', how='left')
            bed_file.clusters.fillna(bed_file.qname, inplace=True)
            bed_file.to_csv(f'{basename}.mappings.cluster.bed', index=False, sep='\t')
            regions.to_csv(f'{args["out"]}/cluster/regions.bed', index=False, sep='\t', header=False)


            print('Creating consensus sequences')
            if not os.path.exists(f'{args["out"]}/cluster/consensus_seq'):
                os.mkdir(f'{args["out"]}/cluster/consensus_seq')

            primer = args['primers']
            primer = [value.strip() for value in primer]
            cluster_specs = cluster.make_consensus_seq(subgraphs, args["out"], args["name"], bed_file, primer)

            cat = (f'cat {args["out"]}/cluster/consensus_seq/{args["name"]}.cluster*.n_reads*.cons.fa > {args["out"]}/cluster/{args["name"]}.cluster.consensus.fa ; rm -rf {args["out"]}/cluster/consensus_seq/')
            subprocess.run(cat, shell=True)

            # run the cons seqs through the pipeline again
            print(f'Filtering reads: {args["out"]}/cluster/{args["name"]}.cluster.consensus.fa', file=sys.stderr)
            fs = glob.glob(f'{args["out"]}/cluster/{args["name"]}.cluster.consensus.fa')

            jobs = []
            for pth in fs:
                jobs.append((pth, primers_target, lock, filter_counts, args['keep_temp'], args['trim_threshold']))
            if args['procs'] > 1:
                with multiprocessing.Pool(args['procs']) as p:
                    p.map(find_reads_with_primers.func, jobs)
            else:
                for j in jobs:
                    find_reads_with_primers.func(j)

            print('Filter counts: ', filter_counts, file=sys.stderr)

            subprocess.run(f'cat {args["out"]}/cluster/*.no_primers.fq > {args["out"]}/cluster/{args["name"]}.cons.without_primers.fq', shell=True)
            subprocess.run(f'rm {args["out"]}/cluster/*.no_primers.fq', shell=True)


            d = "cat {out}/cluster/*.primers_labelled.fq | " \
                "bwa mem -c 1000 -A2 -B3 -O5 -E2 -T0 -L0 -D 0.25 -r 1.25 -d 200 -k 11 -a -t{procs} {ref} - |" \
                "dodi {bias_params} --paired False -c 1 -u 21 --ol-cost 2 --max-overlap 50000 - |" \
                "samtools view -bh - |" \
                "samtools sort -o {out}/cluster/{name}.bwa_dodi_cluster.bam; " \
                "samtools view {out}/{name}.bwa_dodi.bam -b -h -o {out}/{name}.bwa_dodi_removed.bam -U {out}/{name}.bwa_dodi_unremoved.bam -L {out}/cluster/regions.bed; " \
                "samtools merge -o {out}/{name}.bwa_dodi_cluster_merged.bam {out}/{name}.bwa_dodi_unremoved.bam {out}/cluster/{name}.bwa_dodi_cluster.bam; " \
                "samtools index {out}/{name}.bwa_dodi_cluster_merged.bam".format(bias_params=bias_params,
                                                                                 out=args['out'],
                                                                                 name = args['name'],
                                                                                 procs=args['procs'],
                                                                                 ref=args['ref'])

            subprocess.run(d, shell=True)
            subprocess.run(f'rm {args["out"]}/cluster/{args["name"]}.bwa_dodi_cluster.bam {args["out"]}/{args["name"]}.bwa_dodi_unremoved.bam {args["out"]}/{args["name"]}.bwa_dodi.bam {args["out"]}/{args["name"]}.bwa_dodi.bam.bai {args["out"]}/cluster/regions.bed', shell=True)
            # add to filter_counts
            file_list = glob.glob(f'{basename}.filter_counts_summary.csv')

            if file_list:
                # File exists, open it and append a new line
                with open(f'{basename}.filter_counts_summary.csv', 'a') as fc:
                    fc.write(','.join([str(k) for k in filter_counts.values()]) + '\n')
            else:
                # File does not exist, create it and add a new line
                with open(f'{basename}.filter_counts_summary.csv', 'w') as fc:
                    fc.write('Filter counts:' + '\n')
                    fc.write(','.join([str(k) for k in filter_counts.keys()]) + '\n')
                    fc.write(','.join([str(k) for k in filter_counts.values()]) + '\n')

            # need to add a summary file
            with open(f'{basename}.filter_counts_summary.csv', 'a') as fc:
                fc.write('Summary:' + '\n')


        print('fslr finished')
