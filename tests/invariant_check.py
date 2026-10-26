#!/usr/bin/env python3
"""invariant_check.py — guard the gibd merge-output contract.

Checks the properties a user is promised for ANY threshold choice:
  no-overlap : for each unordered haplotype pair, no reported segment overlaps
               another (containment is a special case, so this also forbids one
               segment nested inside another).

Usage:
  invariant_check.py overlap <graph.tsv>
  invariant_check.py monotone <looser.tsv> <tighter.tsv>   # looser must be a
      span-superset of tighter (loosening a selection threshold never drops a call)

Exits 0 (PASS) / 1 (FAIL).
"""
import sys
from collections import defaultdict


def load(path):
    pairs = defaultdict(list)
    with open(path) as f:
        for ln in f:
            if ln.startswith('sample_a'):
                continue
            c = ln.rstrip('\n').split('\t')
            if len(c) < 12:
                continue
            a, b = c[0] + '/' + c[1], c[2] + '/' + c[3]
            pairs[tuple(sorted((a, b)))].append((int(c[5]), int(c[6])))
    return pairs


def check_overlap(path):
    pairs = load(path)
    bad = 0
    for k, segs in pairs.items():
        segs.sort()
        for i in range(len(segs) - 1):
            if segs[i + 1][0] < segs[i][1]:
                bad += 1
    n = sum(len(v) for v in pairs.values())
    if bad:
        print(f"FAIL: {bad} overlapping/contained segment pair(s) in {n} segments")
        return 1
    print(f"PASS: no overlap/containment ({n} segments)")
    return 0


def check_monotone(looser, tighter):
    lo, ti = load(looser), load(tighter)
    missing = 0
    for k, segs in ti.items():
        for (s, e) in segs:
            if not any(min(e, b) > max(s, a) for a, b in lo.get(k, [])):
                missing += 1
    if missing:
        print(f"FAIL: {missing} tighter-threshold segment(s) not covered by the "
              f"looser run (non-monotone)")
        return 1
    print(f"PASS: looser run is a span-superset of the tighter run (monotone)")
    return 0


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'overlap':
        sys.exit(check_overlap(sys.argv[2]))
    elif mode == 'monotone':
        sys.exit(check_monotone(sys.argv[2], sys.argv[3]))
    else:
        print("unknown mode", mode)
        sys.exit(2)
