#!/usr/bin/env python3
# Mixed-orientation cycle fixture — the case AMY1 exposed that --unroll-coherent's
# per-RUN majority classification cannot make sigma-equivariant: a single maximal
# cyclic run that contains BOTH canonical (+) and anti-canonical (-) steps.
#
# Unit A(2) B(3) C(4), cyclic via back-edge C->A. "Hairpin" haplotypes traverse a
# forward copy then immediately an inverted copy (2+3+4+ 4-3-2-) -> one contiguous
# cyclic run with mixed orientation (a palindromic / inverted-tandem read-through).
# Two hairpins turn at DIFFERENT points, so any layer mis-numbering could collide
# unrelated occurrences -> the verifier (verify_mixed.py) checks SOUNDNESS: no
# reported match may exceed the true longest common substring (fwd, or vs RC for an
# inversion) of the pair. That directly catches FALSE IBD, the dangerous direction.
import random
random.seed(7)
def seq(n): return ''.join(random.choice('ACGT') for _ in range(n))

L=2000; UNIT=600; R=2000
S = {1:seq(L), 2:seq(UNIT), 3:seq(UNIT), 4:seq(UNIT), 5:seq(R)}

def fmt(nodes): return ','.join(f'{n}{o}' for n,o in nodes)
A=(2,'+'); B=(3,'+'); C=(4,'+'); Ai=(2,'-'); Bi=(3,'-'); Ci=(4,'-'); Lf=(1,'+'); Rf=(5,'+')

paths = {
    'h0_F2':  fmt([Lf, A,B,C, A,B,C, Rf]),                 # forward, 2 copies
    'h1_F2':  fmt([Lf, A,B,C, A,B,C, Rf]),                 # = h0 -> HBD
    'h2_I2':  fmt([Lf, Ci,Bi,Ai, Ci,Bi,Ai, Rf]),          # inverted, 2 copies
    'h3_HP1': fmt([Lf, A,B,C, Ci,Bi,Ai, Rf]),             # fwd->inv hairpin (+->- crossing)
    'h4_HP2': fmt([Lf, A,B,C, A,B,C, Ci,Bi,Ai, Rf]),      # fwd->inv hairpin, later turn
    'h5_RHP': fmt([Lf, Ci,Bi,Ai, A,B,C, Rf]),             # inv->fwd hairpin (-->+ crossing):
                                                          # together with h3/h4 makes crossings
                                                          # go BOTH ways -> AMY1-like residual
                                                          # that coherent cannot fully resolve
}

with open('mixed.gfa','w') as f:
    f.write('H\tVN:Z:1.0\n')
    for i,s in S.items(): f.write(f'S\t{i}\t{s}\n')
    for a,b in [(1,2),(2,3),(3,4),(4,2),(4,5)]:
        f.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for name,p in paths.items():
        f.write(f'P\t{name}\t{p}\t*\n')

print("wrote mixed.gfa with paths:")
for n,p in paths.items(): print(" ",n,p)
