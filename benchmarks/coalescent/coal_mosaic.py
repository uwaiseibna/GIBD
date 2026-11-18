#!/usr/bin/env python3
import argparse, gzip, random, sys
import msprime, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("gfa")
ap.add_argument("--out", required=True)
ap.add_argument("--ref", default="CHM13")
ap.add_argument("--samples", type=int, default=20)
ap.add_argument("--ne", type=float, default=200)
ap.add_argument("--census", type=float, default=10)
ap.add_argument("--cm", type=float, default=10.0, help="cM over the reference path")
ap.add_argument("--min-truth-cM", type=float, default=0.0)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--donor-list", help="file of path names allowed as donors")
a = ap.parse_args()
rng = random.Random(a.seed)

op = gzip.open if a.gfa.endswith(".gz") else open
seg_lines, link_lines, nlen, paths = [], [], {}, []
with op(a.gfa, "rt") as f:
    for line in f:
        t = line.rstrip("\n").split("\t")
        if t[0] == "S":
            nlen[t[1]] = len(t[2]); seg_lines.append(line)
        elif t[0] == "L":
            link_lines.append(line)
        elif t[0] == "P":
            paths.append((t[1], [(s[:-1], s[-1]) for s in t[2].split(",")]))
ref_i = next(i for i, (n, _) in enumerate(paths) if n.startswith(a.ref))
ref_name, ref_walk = paths[ref_i]
refpos, seen, pos = {}, {}, 0
for n, o in ref_walk:
    seen[n] = seen.get(n, 0) + 1
    refpos.setdefault(n, pos); pos += nlen[n]
L = pos
refpos = {n: p for n, p in refpos.items() if seen[n] == 1}
ref_order = sorted(refpos, key=refpos.get)
ref_fwd = {n for n, o in ref_walk if o == "+"}

keep = set(open(a.donor_list).read().split()) if a.donor_list else None
donors = []
for i, (name, walk) in enumerate(paths):
    if i == ref_i or (keep is not None and name not in keep):
        continue
    cnt = {}
    for n, o in walk:
        cnt[n] = cnt.get(n, 0) + 1
    anc, last = {}, -1
    for k, (n, o) in enumerate(walk):
        if n in refpos and n in ref_fwd and o == "+" and cnt[n] == 1 and refpos[n] > last:
            anc[n] = k; last = refpos[n]
    if anc and min(refpos[n] for n in anc) < 0.02 * L and last > 0.98 * L:
        donors.append((name, walk, anc))
print(f"ref {ref_name} {L} bp, {len(donors)} usable donors of {len(paths)-1}", file=sys.stderr)

dem = msprime.Demography.isolated_model([a.ne])
dem.add_census(time=a.census)
rate = a.cm / 100 / L
ts = msprime.sim_ancestry(samples=a.samples, demography=dem, sequence_length=L,
                          recombination_rate=rate, random_seed=a.seed)
census = np.where(ts.tables.nodes.time == a.census)[0]
samp = list(ts.samples())
ed = ts.tables.link_ancestors(samp, census)
pieces = {s: [] for s in samp}
spans = {}
for l, r, p, c in zip(ed.left, ed.right, ed.parent, ed.child):
    pieces[int(c)].append((float(l), float(r), int(p)))
    spans.setdefault(int(p), []).append((float(l), float(r)))
for s in samp:
    pieces[s].sort()

ev = sorted((l, r, p) for p, iv in spans.items() for l, r in iv)
active, color = [], {}
for l, r, p in ev:
    active = [(r2, p2) for r2, p2 in active if r2 > l]
    used = {color[p2] for _, p2 in active if p2 in color}
    if p not in color:
        free = [d for d in range(len(donors)) if d not in used]
        if not free:
            sys.exit(f"more than {len(donors)} concurrent census lineages")
        color[p] = rng.choice(free)
    active.append((r, p))
for p, iv in spans.items():
    for l, r in iv:
        for l2, r2, p2 in ev:
            if p2 != p and l2 < r and l < r2 and color[p2] == color[p]:
                sys.exit("donor collision")
print(f"{ts.num_trees} trees, {len(color)} census nodes, "
      f"{len(set(color.values()))} donors used", file=sys.stderr)

def nxt_common(da, db, b, ka):
    for n in ref_order[bisect(b):]:
        if n in da[2] and n in db[2] and da[2][n] >= ka:
            return n
    return None

import bisect as _b
ref_starts = [refpos[n] for n in ref_order]
def bisect(b):
    return _b.bisect_left(ref_starts, b)

out_paths, eff = [], {}
for s in samp:
    pc = []
    for l, r, p in pieces[s]:
        if pc and pc[-1][2] == color[p]:
            pc[-1] = (pc[-1][0], r, color[p], pc[-1][3]); continue
        pc.append((l, r, color[p], p))
    walk, cur, ka, cuts = [], pc[0][2], 0, [0.0]
    for (l, r, d, p) in pc[1:]:
        n = nxt_common(donors[cur], donors[d], l, ka)
        if n is None:
            sys.exit("no common anchor")
        walk += donors[cur][1][ka:donors[cur][2][n]]
        cuts.append(refpos[n]); cur, ka = d, donors[d][2][n]
    walk += donors[cur][1][ka:]
    cuts.append(L)
    eff[s] = [(cuts[i], cuts[i + 1], pc[i][2]) for i in range(len(pc))]
    out_paths.append((s, walk))

ind = {s: ts.node(s).individual for s in samp}
hap = {}
for s in samp:
    i = ind[s]
    hap[s] = 1 + list(ts.individual(i).nodes).index(s)
def nm(s): return f"SIM{ind[s]:03d}", hap[s]

with open(a.out + ".gfa", "w") as g:
    g.write("H\tVN:Z:1.0\n")
    g.writelines(seg_lines); g.writelines(link_lines)
    w = lambda walk: "".join((">" if o == "+" else "<") + n for n, o in walk)
    g.write(f"W\t{a.ref}\t0\tref\t0\t{L}\t{w(ref_walk)}\n")
    for s, walk in out_paths:
        sm, h = nm(s)
        g.write(f"W\t{sm}\t{h - 1}\tref\t0\t{sum(nlen[n] for n, _ in walk)}\t{w(walk)}\n")
with open(a.out + ".map", "w") as m:
    for k, bp in enumerate(np.linspace(0, L, 201).astype(int)):
        m.write(f"ref\tm{k}\t{bp * a.cm / L:.6f}\t{bp}\n")
with open(a.out + ".donors", "w") as f:
    for s in samp:
        sm, h = nm(s)
        for l, r, d in eff[s]:
            f.write(f"{sm}\t{h}\t{int(l)}\t{int(r)}\t{donors[d][0]}\n")
nt = 0
with open(a.out + ".truth", "w") as f:
    for i, x in enumerate(samp):
        for y in samp[i + 1:]:
            iv = []
            for l1, r1, d1 in eff[x]:
                for l2, r2, d2 in eff[y]:
                    if d1 == d2 and max(l1, l2) < min(r1, r2):
                        iv.append([max(l1, l2), min(r1, r2)])
            iv.sort(); mg = []
            for l, r in iv:
                if mg and l <= mg[-1][1]:
                    mg[-1][1] = max(mg[-1][1], r)
                else:
                    mg.append([l, r])
            for l, r in mg:
                cm = (r - l) * a.cm / L
                if cm >= a.min_truth_cM:
                    (sa, ha), (sb, hb) = nm(x), nm(y)
                    f.write(f"{sa}\t{ha}\t{sb}\t{hb}\t{int(l)}\t{int(r)}\t{cm:.4f}\n"); nt += 1
print(f"wrote {a.out}.gfa ({len(out_paths)} sample paths), {nt} truth segments", file=sys.stderr)
