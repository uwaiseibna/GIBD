#!/usr/bin/env python3
import sys
truth, calls = sys.argv[1], sys.argv[2]
min_cm = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
ref = sys.argv[4] if len(sys.argv) > 4 else "CHM13"
def key(a, ha, b, hb):
    x, y = (a, int(ha)), (b, int(hb))
    return (x, y) if x <= y else (y, x)
T = {}
for line in open(truth):
    a, ha, b, hb, s, e, cm = line.split()
    T.setdefault(key(a, ha, b, hb), []).append((int(s), int(e), float(cm)))
C = {}
for line in open(calls):
    t = line.rstrip("\n").split("\t")
    if t[0] == "sample_a" or t[0] == ref or t[2] == ref:
        continue
    C.setdefault(key(t[0], int(t[1]), t[2], int(t[3])), []).append((int(t[5]), int(t[6]), float(t[7]), t[8:]))
ov = lambda a, b: max(0, min(a[1], b[1]) - max(a[0], b[0]))
tt = [(k, s) for k, v in T.items() for s in v if s[2] >= min_cm]
hit = sum(1 for k, s in tt if any(ov(s, c) >= 0.5 * (s[1] - s[0]) for c in C.get(k, [])))
cc = [(k, c) for k, v in C.items() for c in v]
tp = sum(1 for k, c in cc if any(ov(c, s) > 0 for s in T.get(k, [])))
tbp = sum(s[1] - s[0] for _, s in tt)
cov = sum(sum(ov(s, c) for c in C.get(k, [])) for k, s in tt)
cbp = sum(c[1] - c[0] for _, c in cc)
inb = sum(sum(ov(c, s) for s in T.get(k, [])) for k, c in cc)
fp_pairs = sorted({k for k, c in cc if not any(ov(c, s) > 0 for s in T.get(k, []))})
print(f"truth>={min_cm}cM: {len(tt)} segs  recall(>=50% covered) {hit}/{len(tt)} = {hit/max(1,len(tt)):.3f}  bp-recall {cov/max(1,tbp):.3f}")
print(f"calls: {len(cc)} segs  precision(overlap truth) {tp}/{len(cc)} = {tp/max(1,len(cc)):.3f}  bp-precision {inb/max(1,cbp):.3f}  FP pairs {len(fp_pairs)}")
for k in fp_pairs[:8]:
    print("  FP", k, C[k][:2])
print("ROW", len(tt), hit, len(cc), tp, tbp, cov, cbp, inb)
