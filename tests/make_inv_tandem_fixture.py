#!/usr/bin/env python3
# Inverted-tandem (both-orientation) multi-node cycle fixture.
# Same repeat unit A(2) B(3) C(4) as make_tandem_fixture.py, but the array is
# carried in BOTH orientations across haplotypes: some haps traverse the unit
# forward (2+ 3+ 4+, back-edge 4->2), some carry it INVERTED (4- 3- 2-, chained
# 2- -> 4-). So node ids 2/3/4 appear as '+' in one stored thread and '-' in
# another, as in inverted repeats of real graphs (e.g. the AMY1 locus). gibd must
# report the inverted matches across the array (checked in `make test`).
import random
random.seed(7)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))

L=2000; UNIT=600; R=2000
S = {1:seq(L), 2:seq(UNIT), 3:seq(UNIT), 4:seq(UNIT), 5:seq(R)}

def fwd(copies):              # L + (A B C)+ * copies + R        (forward array)
    nodes=[1]
    for _ in range(copies): nodes += [(2,'+'),(3,'+'),(4,'+')]
    nodes=[(1,'+')]+nodes[1:]+[(5,'+')]
    return nodes

def inv(copies):              # L + (C B A)- * copies + R        (inverted array)
    nodes=[(1,'+')]
    for _ in range(copies): nodes += [(4,'-'),(3,'-'),(2,'-')]
    nodes += [(5,'+')]
    return nodes

def fmt(nodes): return ','.join(f'{n}{o}' for n,o in nodes)

paths = {
    'h1_fwd2': fmt(fwd(2)),   # forward, 2 copies
    'h2_fwd2': fmt(fwd(2)),   # identical to h1 -> forward HBD
    'h3_inv2': fmt(inv(2)),   # SAME array inverted -> inversion vs h1/h2
    'h4_inv2': fmt(inv(2)),   # identical to h3 -> HBD between the inverted haps
    'h5_fwd3': fmt(fwd(3)),   # forward, 3 copies (phase/copy-number offset)
    'h6_inv3': fmt(inv(3)),   # inverted, 3 copies (offset on the inverted strand)
}

with open('inv_tandem.gfa','w') as f:
    f.write('H\tVN:Z:1.0\n')
    for i,s in S.items(): f.write(f'S\t{i}\t{s}\n')
    # forward unit links (their reverses, added by the bidirectional GBWT, carry
    # the inverted traversals); L-lines are informational — gfa2gbwt reads P-lines.
    for a,b in [(1,2),(2,3),(3,4),(4,2),(4,5)]:
        f.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for name,p in paths.items():
        f.write(f'P\t{name}\t{p}\t*\n')

print("wrote inv_tandem.gfa with paths:")
for n,p in paths.items(): print(" ",n,p)
