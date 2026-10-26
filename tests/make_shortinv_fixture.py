#!/usr/bin/env python3
# Short inversion inside long IBD (IBD segments, inverted-match floor).
# Thirty 2 kb reference nodes on a 1 cM/kb map; nodes 14-16 (bp 26000-32000) are
# an inversion carried by two haplotypes (one locus). With --min-extend-cM 10 the
# forward floor is 10 kb, longer than the 6 kb inversion: the inverted match is
# found only because inverted matches have their own floor (max(100, max-gap-bp)).
#   ref      1..30 forward (reference path, --ref-path 0)
#   inv_a/b  1..13,16-,15-,14-,17..30   inversion carriers
#   fwd_a/b  1..30                      reference orientation
# Expected: inv_a x fwd_a = ONE segment 0-60000, inverted=1, OPP:26000-32000
# (IBD across the orientation difference); inv_a x inv_b and fwd_a x fwd_b = one
# segment 0-60000, inverted=0.
import random
random.seed(5)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))
N = 30
S = {i: seq(2000) for i in range(1, N + 1)}
fwd = [(i, '+') for i in range(1, N + 1)]
inv = [(i, '+') for i in range(1, 14)] + [(16, '-'), (15, '-'), (14, '-')] + [(i, '+') for i in range(17, N + 1)]
paths = {'ref': fwd, 'inv_a': inv, 'inv_b': inv, 'fwd_a': fwd, 'fwd_b': fwd}
with open('shortinv.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for i, s in S.items(): f.write(f'S\t{i}\t{s}\n')
    for name, p in paths.items():
        f.write(f"P\t{name}\t{','.join(f'{n}{o}' for n, o in p)}\t*\n")
with open('shortinv.map', 'w') as f:
    for k in range(0, 2 * N + 1):
        f.write(f'chr1\t.\t{k * 2.0:.1f}\t{k * 2000}\n')
