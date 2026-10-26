#!/usr/bin/env python3
# Orientation-aware IBS runs (merge mode). Thirty 2 kb reference nodes on a 1 cM/kb map; nodes 14-16
# (bp 26000-32000) are an inversion carried by inv_a and inv_b (one locus). Nodes 31, 32 are alternative
# alleles of nodes 12 and 18 (a discordance of 2 kb, more than max-gap 1000 bp).
#   ref, fwd_a  1..30 forward
#   inv_a/b     1..13,16-,15-,14-,17..30
#   n12         fwd with node 12 -> 31            (discordance at 22000-24000)
#   n12_18      fwd with 12 -> 31 and 18 -> 32    (discordances at 22000-24000 and 34000-36000)
#   n18         fwd with 18 -> 32                 (second carrier of node 32, so min-mac 2 keeps it)
# Expected (extend/seed/output as in the Makefile):
#   short flank: inv_a x n12 (extend 10, seed 10): 0-22000 forward, and ONE run 24000-60000 = 2 kb stub +
#     inversion + right flank, inverted=1 OPP:26000-32000 (the 2 kb stub < min-extend joins the run);
#   no flank seed: inv_a x fwd_a (seed 30, extend 10): flanks 26 and 28 cM, neither a seed; the run 0-60000
#     (60 cM) seeds -> one segment, inverted=1, OPP;
#   run below min-extend: inv_a x n12_18 (extend 12): stub 24000-26000 + inversion + stub 32000-34000 = 10 cM
#     < 12 -> not an extension: 0-22000 and 36000-60000, inverted=0, no OPP.
import random
random.seed(9)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))
N = 30
S = {i: seq(2000) for i in range(1, N + 3)}
fwd = [(i, '+') for i in range(1, N + 1)]
inv = [(i, '+') for i in range(1, 14)] + [(16, '-'), (15, '-'), (14, '-')] + [(i, '+') for i in range(17, N + 1)]
def sub(p, k, alt): return [(alt, '+') if n == k else (n, o) for n, o in p]
paths = {'ref': fwd, 'inv_a': inv, 'inv_b': inv, 'fwd_a': fwd, 'n12': sub(fwd, 12, 31),
         'n12_18': sub(sub(fwd, 12, 31), 18, 32), 'n18': sub(fwd, 18, 32)}
with open('orirun.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for i, s in S.items(): f.write(f'S\t{i}\t{s}\n')
    for name, p in paths.items():
        f.write(f"P\t{name}\t{','.join(f'{n}{o}' for n, o in p)}\t*\n")
with open('orirun.map', 'w') as f:
    for k in range(0, 2 * N + 1):
        f.write(f'chr1\t.\t{k * 2.0:.1f}\t{k * 2000}\n')
