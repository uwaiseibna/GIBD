#!/usr/bin/env python3
# Forward-coherent multi-node tandem-CNV fixture (NO inverted repeat).
# Unit = A(2) B(3) C(4) traversed a variable number of times by each haplotype,
# ALL in + orientation. Union edges L->A, A->B, B->C, C->A, C->R make {A,B,C} a
# multi-node SCC with zero both-orientation traversal — the case unrollCycles
# should resolve to a DAG (SCC fallback 0) if the obstruction is purely about
# inverted repeats.
import random
random.seed(7)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))

L=2000; UNIT=600; R=2000
S = {1:seq(L), 2:seq(UNIT), 3:seq(UNIT), 4:seq(UNIT), 5:seq(R)}

def path(copies):  # L + (A B C)*copies + R, all forward
    nodes=[1]
    for _ in range(copies): nodes += [2,3,4]
    nodes += [5]
    return ','.join(f'{n}+' for n in nodes)

paths = {
    'h1_2copy': path(2),
    'h2_3copy': path(3),
    'h3_2copy': path(2),   # identical to h1 -> full IBD expected
    'h4_1copy': path(1),
    'h5_2copy': path(2),   # identical to h1/h3
}

with open('tandem.gfa','w') as f:
    f.write('H\tVN:Z:1.0\n')
    for i,s in S.items(): f.write(f'S\t{i}\t{s}\n')
    for a,b in [(1,2),(2,3),(3,4),(4,2),(4,5)]:
        f.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for name,p in paths.items():
        f.write(f'P\t{name}\t{p}\t*\n')
print("wrote tandem.gfa with paths:")
for n,p in paths.items(): print(" ",n,p)
