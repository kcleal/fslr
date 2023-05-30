import sys

min_cov = float(sys.argv[1])
i = 0
c = 0
for r in sys.stdin:
    i += 1
    d = r.strip().split('\t')
    if d[0] == d[5]:
        continue  # same name

    matching_bases = int(d[9])
    ql = int(d[1])
    qs = int(d[2])
    qe = int(d[3])
    tl = int(d[6])
    ts = int(d[7])
    te = int(d[8])
    query_coverage = (qe - qs) / ql
    target_coverage = (te - ts) / tl

    qname = d[0].split('.')[0]
    tname = d[5].split('.')[0]

    # t = ('ccc26fc6-89a8-b49a-20d2-36c488b3e6be', '934a9293-459f-4146-ae0d-fc3b37c06e94')
    # if qname in t and tname in t:
    #     print(r, query_coverage, target_coverage, matching_bases / ql, matching_bases / tl)

    if query_coverage < min_cov or \
            target_coverage < min_cov or \
            (matching_bases / ql) < 0.5 or \
            (matching_bases / tl) < 0.5:
        continue
    sys.stdout.write(r)
    c += 1

print(f'Input alignmnets: {i}, kept alignments {c}', file=sys.stderr)
