import pysam 
import pandas as pd
from collections import defaultdict
from importlib.metadata import version


def get_query_pos_from_cigartuples(r):
    # Infer the position on the query sequence of the alignment using cigar string
    start = 0
    query_length = r.infer_read_length()  # Note, this also counts hard-clips
    end = query_length
    if r.cigartuples[0][0] == 4 or r.cigartuples[0][0] == 5:
        start += r.cigartuples[0][1]
    if r.cigartuples[-1][0] == 4 or r.cigartuples[-1][0] == 5:
        end -= r.cigartuples[-1][1]
    return start, end, query_length


def mapping_info(f, outf, regions_path, primers):
    flsr_version = version("fslr")

    af = pysam.AlignmentFile(f, 'r')
    d = defaultdict(list)
    for a in af.fetch(until_eof=True):
        if not a.flag & 4:
            d[a.qname].append(a)

    regions = defaultdict(list)
    if regions_path:
        with open(regions_path, 'r') as rgns:
            for line in rgns:
                l = line.strip().split('\t')
                chrom = l[0]
                start = int(l[1])
                end = int(l[2])
                regions[chrom].append(pd.Interval(left=start, right=end))

    res = []
    no = 0
    yes = 0
    for qname, v in d.items():
        flag = [(index, i) for index, i in enumerate(v) if not i.flag & 2304]
        if len(flag) > 1:  # todo check bug in dodi, not currently setting primary alignment flag properly
            flag = [flag[flag.index(max(flag, key=lambda x: x[1].get_tag('AS')))]]

        if len(flag) != 1:
            print('Error in ', f, 'flag problem', len(flag), [i.flag for i in v])
            quit()
        pri_index, pri_read = flag[0]
        primary_reverse = bool(pri_read.flag & 16)
        seq = pri_read.get_forward_sequence()
        n_aligns = len(v)
        any_seq = False

        temp = []
        for index, a in enumerate(v):
            qstart, qend, qlen = get_query_pos_from_cigartuples(a)
            align_reverse = bool(a.flag & 16)
            if primary_reverse != align_reverse:
                start_temp = qlen - qend
                qend = start_temp + qend - qstart
                qstart = start_temp
            pri = index == pri_index
            if not pri:
                no += 1
            else:
                yes += 1
                any_seq = len(seq) if seq else 0

            chrom = af.get_reference_name(a.rname)
            start = a.reference_start + 1
            end = a.reference_end
            t = pd.Interval(start, end)
            if regions and chrom in regions and any(t.overlaps(q) for q in regions[chrom]):
                overlaps = 1
            else:
                overlaps = 0

            rd = {'qname': a.qname,
                 'n_alignments': n_aligns,
                 'chrom': chrom,
                 'rstart': start,
                 'rend': end,
                 'strand': '-' if align_reverse else '+',
                 'qstart': qstart,
                 'qend': qend,
                 'qlen': qlen,
                 'aln_size': qend - qstart,
                 'mapq': a.mapq,
                 'alignment_score': a.get_tag('AS'),
                 'seq': seq if pri else '',
                 'fslr_version': flsr_version,
                 'inferred_by_primer': 0,
                 }

            if regions:
                rd['overlaps_region'] = overlaps

            temp.append(rd)

        if not any_seq:
            print('missing', qname, [(len(vv.seq), vv.infer_query_length()) if vv.seq else vv.infer_query_length() for vv in v])
            quit()

        if len(temp) > 1:
            res += temp
            continue

        # Add in inferred alignments at the ends 'missing bread'
        temp = sorted(temp, key=lambda r: r['qstart'])
        qstart_gap = temp[0]['qstart']
        qend_gap = temp[0]['qlen'] - temp[0]['qend']
        p_names = temp[0]['qname'].split('.')[-1].split('_')
        p1, p2 = [i.rstrip('FR') for i in p_names]

        if qstart_gap > 5 and qend_gap > 5:  # gaps at both ends, skip these
            res += temp
            continue
        # This was considered but deemed unnecessary as primers have already been identified
        # != 'False' and len(primers[p2]):
        if p1 != 'False':
            temp[0]['n_alignments'] = 2
            p1_rd = {'qname': temp[0]['qname'],
                     'n_alignments': 2,
                     'chrom': p1,
                     'rstart': 0,
                     'rend': 0,
                     'strand': '-' if p_names[0][-1] == "R" else '+',
                     'qstart': 0,
                     'qend': len(primers[p1]),
                     'qlen': qlen,
                     'aln_size': 0,
                     'mapq': 0,
                     'alignment_score': 0,
                     'seq': '',
                     'fslr_version': flsr_version,
                     'inferred_by_primer': 1,
                     }
            temp = [p1_rd, temp[0]]
        elif p2 != 'False':
            temp[0]['n_alignments'] = 2
            p2_rd = {'qname': temp[0]['qname'],
                     'n_alignments': 2,
                     'chrom': p2,
                     'rstart': 0,
                     'rend': 0,
                     'strand': '-' if p_names[1][-1] == "R" else '+',
                     'qstart': qlen - len(primers[p2]),
                     'qend': qlen,
                     'qlen': qlen,
                     'aln_size': 0,
                     'mapq': 0,
                     'alignment_score': 0,
                     'seq': '',
                     'fslr_version': flsr_version,
                     'inferred_by_primer': 1,
                     }
            temp = [temp[0], p2_rd]

        res += temp

    df = pd.DataFrame.from_records(res).sort_values(['qname', 'qstart'])

    bad_anchors = []
    # flag reads with small anchoring alignments
    for grp, d in df.groupby('qname'):
        aln_s = list(d['aln_size'])
        if aln_s[0] < 50 or aln_s[-1] < 50:
            bad_anchors += [1] * len(d)
        else:
            bad_anchors += [0] * len(d)
    df['short_anchor<50bp'] = bad_anchors

    df = df.sort_values(['n_alignments', 'qname', 'qstart'], ascending=[False, True, True])

    cols = ['chrom', 'rstart', 'rend', 'qname', 'n_alignments', 'aln_size', 'qstart', 'qend', 'strand', 'mapq', 'qlen',
             'alignment_score', 'short_anchor<50bp', 'fslr_version', 'inferred_by_primer', 'seq']
    if regions:
        cols.append('overlaps_region')
    df = df[cols]
    df.to_csv(outf, index=False, sep="\t")


if __name__ == '__main__':
    import argparse
    parse = argparse.ArgumentParser()
    parse.add_argument('--bam', help='bam file to assess')
    parse.add_argument('--out', help='out put bed file')
    args = parse.parse_args()
    mapping_info(args.bam, args.out)
    print('Done')
