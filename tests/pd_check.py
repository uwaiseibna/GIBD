#!/usr/bin/env python3
"""pd_check.py — regression check for mismatch-tolerant (--max-gap-cM) chaining.

Reads a sim_panel.py panel with within-IBD PRIVATE DELETIONS (indels invisible to a
SNP-only caller) and the two graph outputs (chained on the bp gap and on the cM gap). Asserts:
  - cM-gap RECALL of true founder-IBD == 1.0 and strictly beats bp-gap (which breaks
    at the private-deletion reference-bp gaps it cannot bridge),
  - both modes PRECISION == 1.0 against the full founder ground truth (no false IBD).
Exits 0 (PASS) / 1 (FAIL). Uses only <out>.truth and <out>.copy (no java/hap-ibd).

Usage: pd_check.py <out>.truth <out>.copy <bpgap.tsv> <cmgap.tsv> [min_cM]
"""
import sys, re
truth_f, copy_f, bp_f, cm_f = sys.argv[1:5]
MINCM = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5

def norm(i, j): return (i, j) if i <= j else (j, i)

truth = []
with open(truth_f) as f:
    for ln in f:
        if ln[0] == '#' or ln.startswith('hap_i'): continue
        c = ln.split()
        if float(c[7]) >= MINCM:
            truth.append((norm(int(c[0]), int(c[1])), int(c[4]), int(c[5])))

site_bp, rows = None, []
with open(copy_f) as f:
    for ln in f:
        if ln.startswith('# site_bp:'):
            site_bp = [int(x) for x in ln.split(':')[1].split()]
        elif ln[0] != '#':
            rows.append([int(x) for x in ln.split()])
S = len(site_bp)

def load(path):
    segs = {}
    with open(path) as f:
        for ln in f:
            if ln.startswith('sample_a'): continue
            c = ln.rstrip('\n').split('\t')
            if len(c) < 8: continue
            ma, mb = re.search(r'smp(\d+)', c[0]), re.search(r'smp(\d+)', c[2])
            if not ma or not mb: continue
            segs.setdefault(norm(int(ma.group(1)), int(mb.group(1))), []).append(
                (int(c[5]), int(c[6])))
    return segs

def recall(segs):
    hit = sum(1 for pr, s, e in truth
              if any(min(e, b) > max(s, a) for a, b in segs.get(pr, [])))
    return hit / len(truth)

def precision(segs):
    tp = fp = 0
    for pr, lst in segs.items():
        i, j = pr
        for (sbp, ebp) in lst:
            sl = [s for s in range(S) if sbp <= site_bp[s] <= ebp]
            if not sl: continue
            ag = sum(1 for s in sl if rows[i][s] == rows[j][s]) / len(sl)
            tp, fp = (tp + 1, fp) if ag >= 0.5 else (tp, fp + 1)
    return tp / max(tp + fp, 1)

bp, cm = load(bp_f), load(cm_f)
r_bp, r_cm = recall(bp), recall(cm)
p_bp, p_cm = precision(bp), precision(cm)
ok = (abs(r_cm - 1.0) < 1e-9 and r_cm > r_bp + 0.05
      and abs(p_bp - 1.0) < 1e-9 and abs(p_cm - 1.0) < 1e-9)
msg = (f"cM-gap recall {r_cm:.3f} (> bp-gap {r_bp:.3f}); "
       f"precision bp={p_bp:.3f} cM={p_cm:.3f}")
if ok:
    print("PASS: mismatch-tolerant chaining recovers indel-fragmented IBD -- " + msg)
    sys.exit(0)
print("FAIL: " + msg); sys.exit(1)
