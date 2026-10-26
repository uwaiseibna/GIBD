#!/usr/bin/env python3
# Shared inversion allele fixture (merge-mode column inv_allele).
# Ten 2 kb reference nodes; nodes 4-6 (bp 6000-12000) are an inversion carried by
# two haplotypes, so the panel has one inversion locus there. Node 8 is inverted
# on ONE haplotype only (private, not a locus).
#   ref      1..10 forward (reference path, --ref-path 0)
#   inv_a/b  1,2,3,6-,5-,4-,7..10   both inverted   -> INV:6000-12000
#   fwd_a/b  1..10                  both reference  -> REF:6000-12000
#   inv_rc   inv_a stored reverse-complemented (strand, not orientation): its
#            own direction is reverse, so 4+,5+,6+ is its inverted run
#   priv     1..7,8-,9,10           private inversion (1 carrier: no locus)
#   ins_f    1..8,11,9,10          2 kb insertion (node 11) after node 8
#   ins_rc   ins_f stored reverse-complemented
#   idup_a/b 1..8,3-,2-,9,10         inverted DUPLICATION of nodes 2-3 (not in place)
# Expected: inv_a x inv_b = INV only; inv_a x inv_rc = INV (opposite strands);
# fwd_a x fwd_b = REF only (no entry at the private node 8); inv_a x fwd_a and
# inv_rc x fwd_a = inverted=1, OPP:6000-12000 (opposite alleles at the locus); ins_f x ins_rc = insertion=1;
# idup_a x idup_b = insertion=1, inverted=0 and NO INV entry (an inverted copy is
# not an inversion and not an orientation difference).
import random
random.seed(11)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))
N = 10
S = {i: seq(2000) for i in range(1, N + 2)}
fwd = [(i, '+') for i in range(1, N + 1)]
inv = [(1, '+'), (2, '+'), (3, '+'), (6, '-'), (5, '-'), (4, '-')] + [(i, '+') for i in range(7, N + 1)]
rc = [(n, '-' if o == '+' else '+') for n, o in reversed(inv)]
priv = [(i, '+') for i in range(1, 8)] + [(8, '-'), (9, '+'), (10, '+')]
insf = [(i, '+') for i in range(1, 9)] + [(11, '+'), (9, '+'), (10, '+')]
rcf = lambda p: [(n, '-' if o == '+' else '+') for n, o in reversed(p)]
paths = {'ref': fwd, 'inv_a': inv, 'inv_b': inv, 'fwd_a': fwd, 'fwd_b': fwd, 'inv_rc': rc, 'priv': priv,
         'ins_f': insf, 'ins_rc': rcf(insf),
         'idup_a': [(i, '+') for i in range(1, 9)] + [(3, '-'), (2, '-'), (9, '+'), (10, '+')],
         'idup_b': [(i, '+') for i in range(1, 9)] + [(3, '-'), (2, '-'), (9, '+'), (10, '+')]}
with open('invallele.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for i, s in S.items(): f.write(f'S\t{i}\t{s}\n')
    for name, p in paths.items():
        f.write(f"P\t{name}\t{','.join(f'{n}{o}' for n, o in p)}\t*\n")
with open('invallele.map', 'w') as f:
    f.write('chr1\t.\t0.0\t0\nchr1\t.\t20.0\t20000\n')
