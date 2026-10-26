#!/usr/bin/env python3
# Soundness verifier for cyclic-region sweep output (catches FALSE IBD).
#
# For each reported match (hap_a, hap_b, inverted, len_bp), the two haplotypes must
# ACTUALLY share a subpath of >= that bp length: a forward common substring when
# inverted=0, or a common substring of a vs reverse-complement(b) when inverted=1.
# We bound every match by the TRUE longest common substring (in bp, over oriented
# nodes) of the appropriate orientation. A match longer than that bound is IBD the
# sweep fabricated -> FAIL. (Under-calling/missing matches is conservative = OK.)
#
# Usage: verify_mixed.py <graph.gfa> <sweep.tsv> [label]
import sys

gfa, tsv = sys.argv[1], sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else tsv

nlen, paths, order = {}, {}, []
with open(gfa) as f:
    for line in f:
        t = line.rstrip('\n').split('\t')
        if t[0] == 'S':
            nlen[int(t[1])] = len(t[2])
        elif t[0] == 'P':
            seq = [(int(s[:-1]), s[-1]) for s in t[2].split(',')]
            paths[len(order)] = seq          # path index == GBWT hap id (P-line order)
            order.append(t[1])

def rc(seq):                                  # reverse-complement of an oriented walk
    return [(nid, '-' if o == '+' else '+') for nid, o in reversed(seq)]

def lcs_bp(A, B):                             # longest common contiguous oriented run, in bp
    nB = len(B); prev = [0]*(nB+1); best = 0
    for i in range(len(A)):
        curr = [0]*(nB+1)
        ai = A[i]
        for j in range(nB):
            if ai == B[j]:
                curr[j+1] = prev[j] + nlen[ai[0]]
                if curr[j+1] > best: best = curr[j+1]
        prev = curr
    return best

# Bounds per pair (memoized).
fwd_bound, rc_bound = {}, {}
def bounds(a, b):
    if (a, b) not in fwd_bound:
        fwd_bound[(a, b)] = lcs_bp(paths[a], paths[b])
        rc_bound[(a, b)]  = lcs_bp(paths[a], rc(paths[b]))
    return fwd_bound[(a, b)], rc_bound[(a, b)]

violations, n = [], 0
with open(tsv) as f:
    header = f.readline()
    for line in f:
        c = line.split('\t')
        if len(c) < 6: continue
        a, b, inv, lbp = int(c[0]), int(c[1]), int(c[2]), int(c[5])
        if a > b: a, b = b, a
        n += 1
        fb, rb = bounds(a, b)
        bound = rb if inv else fb
        if lbp > bound:
            violations.append((a, b, inv, lbp, bound))

print(f"[{label}] {n} reported matches checked against true longest-common-substring bounds")
if violations:
    print(f"  FAIL: {len(violations)} FALSE-IBD violation(s) (reported len_bp > genuine shared substring):")
    for a, b, inv, lbp, bound in violations[:10]:
        print(f"    {order[a]} vs {order[b]} inverted={inv}: reported {lbp} bp > genuine {bound} bp")
    sys.exit(1)
print("  PASS: every reported match is backed by a genuine shared subpath (no false IBD)")
