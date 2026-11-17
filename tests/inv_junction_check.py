import sys
copy_f, tsv = sys.argv[1:3]
site_bp, rows = None, []
for ln in open(copy_f):
    if ln.startswith('# site_bp:'): site_bp = [int(x) for x in ln.split(':')[1].split()]
    elif ln[0] != '#': rows.append(ln.split())
n = bad = 0
for ln in open(tsv):
    c = ln.split('\t')
    if c[0] == 'sample_a': continue
    i, j = int(c[0].split('#')[0][3:]), int(c[2].split('#')[0][3:])
    sl = [s for s, b in enumerate(site_bp) if int(c[5]) <= b <= int(c[6])]
    n += 1
    bad += sum(rows[i][s] == rows[j][s] for s in sl) < 0.5 * len(sl)
print(("PASS" if n and not bad else "FAIL") + f": inversion junctions: {bad} of {n} segments copy different founders")
