#!/usr/bin/env python3
"""
sim_panel.py — simulate a many-haplotype variation graph with a segregating
non-reference INSERTION and a segregating INVERSION, for population-scale
IBD-moat measurement against a real hap-ibd run.

Population model: K founders, each a random SNP haplotype over S sites. Each
sample haplotype is a Li-Stephens mosaic that copies a founder and switches with
probability p_switch per site (recombination). Copying the same founder over a
region => true IBD between samples there.
  INSERTION: a subset of founders carry a multi-kb insertion at one locus, so a
    sample copying a carrier founder across it carries the insertion. Both the
    graph and hap-ibd detect the spanning IBD (flank SNPs); the graph additionally
    LABELS co-inheritance (insertion=1). => annotation contribution.
  INVERSION: an orientation-only flip of a region, carried by a random subset of
    samples INDEPENDENTLY of the SNP mosaic (--inv-freq). Two samples can thus
    share the founder background across the region (SNP-IBD) yet differ in
    orientation -- the pair hap-ibd BREAKS on (mirrored ref-coord alleles
    mismatch) but the graph keeps whole via forward-vs-reverse matching
    (inverted=1). => DETECTION moat.

Graph (proper variation graph, so SNP alt-alleles are tiny bubble nodes that are
NOT mistaken for insertions):
    filler F_s (shared backbone, 100 bp) -- [R_s | A_s] (1 bp SNP bubble) -- ...
    reference path (P-line 0) = all ref alleles, no insertion.
    off-reference nodes: A_s (1 bp, < min-insertion-bp) and the INS node (>=50 bp).

Outputs:
    <out>.gfa   variation graph (S-lines + P-lines: ref path 0, then N samples)
    <out>.map   PLINK genetic map (chrom marker cM bp) for cM lengths
    <out>.truth per-PAIR true IBD spans over the insertion locus (founder-IBD
                ground truth; carrier=1 => spans the insertion). For recall.
    <out>.vcf   SNP-only diploid VCF (individual k = sample-haps 2k,2k+1) to run a
                REAL hap-ibd on the identical panel. The INS is absent by design.

Three-way comparison:
    GRAPH (full):  gfa2gbwt <out>.gfa <out> ;
        gibd --gfa <out>.gfa --ref-path 0 --map <out>.map ... <out>.gbwt
    LINEAR view :  same + --no-insertion        (graph-internal linear baseline)
    REAL hap-ibd:  hap-ibd gt=<out>.vcf map=<out>.map out=<out>.hapibd
    Recall each vs <out>.truth's carrier pairs (the insertion-spanning IBD).
    Map sim hap h <-> hap-ibd (ind, hap):  h = 2*ind + (hap-1).
"""
import sys, random, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--sites", type=int, default=1000)
ap.add_argument("--founders", type=int, default=8)
ap.add_argument("--samples", type=int, default=500)     # haplotypes
ap.add_argument("--p-switch", type=float, default=0.02) # recomb per site
ap.add_argument("--maf", type=float, default=0.3)
ap.add_argument("--filler-bp", type=int, default=100)
ap.add_argument("--ins-bp", type=int, default=2000)
ap.add_argument("--ins-site", type=int, default=-1)     # default: middle
ap.add_argument("--carrier-founders", type=int, default=3)  # of K carry the INS
# --- deletion: a reference-PRESENT segment (a dedicated DEL node on the ref path)
#     that a subset of founders LACK. A sample copying such a founder OMITS the DEL
#     node => its walk skips del_bp of reference (a reference-coordinate deletion).
#     hap-ibd sees this as ONE biallelic marker (ref/alt) and ignores its length;
#     the graph sees the skipped reference span (a bp gap) and types it deletion=1.
ap.add_argument("--del-bp", type=int, default=2000)
ap.add_argument("--del-site", type=int, default=-1)     # default: 3S//4
ap.add_argument("--del-carrier-founders", type=int, default=3)  # of K LACK the seg
# --- PRIVATE (per-haplotype) deletions: scattered reference blocks that each
#     haplotype INDEPENDENTLY omits with prob --priv-del-rate. Unlike the founder
#     --del above, these are NOT founder-linked, so two truly-IBD haplotypes (same
#     founder) DIFFER wherever one privately deleted a block the other kept -- a
#     within-IBD-segment indel discordance, the case hap-ibd tolerates (it never
#     sees these: they are EXCLUDED from the VCF, like real indels dropped from a
#     SNP-only callset) but an exact-match graph sweep BREAKS on. The reference-bp
#     gap of such a block (>= --priv-del-bp) exceeds --max-gap-bp, so bp-gap
#     chaining cannot bridge it; cM-gap chaining (the indel is cM-negligible) can.
#     This is the controlled scenario that distinguishes --max-gap-cM from bp-gap.
ap.add_argument("--priv-del-rate", type=float, default=0.0)   # per hap per block (0=off)
ap.add_argument("--priv-del-bp", type=int, default=2000)      # block size (> max-gap-bp)
ap.add_argument("--priv-del-sites", type=int, default=0)      # # scattered blocks
ap.add_argument("--cm-per-mb", type=float, default=1.3)
# --- inversion: an orientation-only flip of a region, carried by a random subset
#     of samples INDEPENDENTLY of the SNP mosaic. Independence is deliberate: it
#     lets two samples share the founder background across the region (SNP-IBD) yet
#     differ in orientation -- the pair a linear method BREAKS on (mirrored alleles
#     mismatch) but the graph keeps whole via forward-vs-reverse matching.
ap.add_argument("--inv-freq", type=float, default=0.3)   # fraction of samples inverted
ap.add_argument("--inv-lo-site", type=int, default=-1)   # region start site (default S//8)
ap.add_argument("--inv-hi-site", type=int, default=-1)   # region end   site (default 3S//8)
ap.add_argument("--seed", type=int, default=7)
ap.add_argument("--out", default="panel_sim")
a = ap.parse_args()
random.seed(a.seed)

S, K, N = a.sites, a.founders, a.samples
ins_site = a.ins_site if a.ins_site >= 0 else S // 2
carriers = set(range(a.carrier_founders))         # founders {0..c-1} carry INS
del_site = a.del_site if a.del_site >= 0 else 3 * S // 4
# deletion carriers = the LAST d founders LACK the reference segment (a DIFFERENT
# founder set than the INS carriers, so the two SVs are not confounded).
del_carriers = set(range(K - a.del_carrier_founders, K))
inv_lo = a.inv_lo_site if a.inv_lo_site >= 0 else S // 8
inv_hi = a.inv_hi_site if a.inv_hi_site >= 0 else 3 * S // 8
inv_mid = (inv_lo + inv_hi) // 2
def mirror(s): return inv_lo + inv_hi - s        # reference-coord allele scramble
samp_inv = [random.random() < a.inv_freq for _ in range(N)]  # per-sample inversion bit

# scattered private-deletion block sites (evenly spread, avoiding the special loci
# so the private indels do not confound the INS/DEL/inversion truth).
_special = {ins_site, del_site} | set(range(inv_lo, inv_hi + 1))
pd_sites = []
if a.priv_del_sites > 0:
    for k in range(a.priv_del_sites):
        s = int((k + 0.5) * S / a.priv_del_sites)
        if 0 <= s < S and s not in _special and s not in pd_sites:
            pd_sites.append(s)
pd_set = set(pd_sites)
# per-haplotype, per-block omission bit: samp_pd_del[h] = {site: True if h omits it}
samp_pd_del = [{s: (random.random() < a.priv_del_rate) for s in pd_sites}
               for _ in range(N)]

# founder SNP alleles
fal = [[1 if random.random() < a.maf else 0 for _ in range(S)] for _ in range(K)]

# node ids: filler F_s = 10*(s+1); ref R_s = +1; alt A_s = +2; INS = 9000001
def F(s): return 10 * (s + 1)
def R(s): return 10 * (s + 1) + 1
def A(s): return 10 * (s + 1) + 2
INS = 9000001
DEL = 9000002   # reference-present segment that deletion-carriers OMIT
def PD(s): return 9100001 + s   # private-deletion block at site s (ref-present)

# bp coordinate of each filler/ref node along the reference (for the map + cM)
# layout: F_0 R_0 F_1 R_1 ... F_{S-1} R_{S-1} F_S ; SNP nodes are 1 bp.
bp = {}
pos = 0
for s in range(S):
    bp[F(s)] = pos; pos += a.filler_bp
    if s == del_site:                  # DEL is reference-present: consumes ref bp
        bp[DEL] = pos; pos += a.del_bp
    if s in pd_set:                    # private-del block is reference-present
        bp[PD(s)] = pos; pos += a.priv_del_bp
    bp[R(s)] = pos; pos += 1
bp[F(S)] = pos; pos += a.filler_bp
total_bp = pos

# ---- emit GFA ----
with open(a.out + ".gfa", "w") as g:
    # S-lines (sequence length is what the sweep uses; content is filler)
    for s in range(S):
        g.write(f"S\t{F(s)}\t{'A'*a.filler_bp}\n")
        g.write(f"S\t{R(s)}\tA\n")
        g.write(f"S\t{A(s)}\tC\n")
    g.write(f"S\t{F(S)}\t{'A'*a.filler_bp}\n")
    g.write(f"S\t{INS}\t{'G'*a.ins_bp}\n")
    g.write(f"S\t{DEL}\t{'T'*a.del_bp}\n")
    for s in pd_sites:
        g.write(f"S\t{PD(s)}\t{'C'*a.priv_del_bp}\n")

    # reference path (path 0): all ref alleles, no insertion, WITH the DEL segment
    # (the reference carries it; deletion-carriers omit it).
    ref_walk = []
    for s in range(S):
        ref_walk.append(f"{F(s)}+")
        if s == del_site:
            ref_walk.append(f"{DEL}+")
        if s in pd_set:
            ref_walk.append(f"{PD(s)}+")
        ref_walk.append(f"{R(s)}+")
    ref_walk.append(f"{F(S)}+")
    g.write("P\tref#0#sim\t" + ",".join(ref_walk) + "\t*\n")

    # samples: Li-Stephens mosaic
    truth = []   # (sample, carries_ins, copy-array) for IBD truth
    samp_alleles = []
    samp_carry = []
    samp_carry_del = []
    for i in range(N):
        cur = random.randrange(K)
        copy = []
        for s in range(S):
            if s > 0 and random.random() < a.p_switch:
                cur = random.randrange(K)
            copy.append(cur)
        alleles = [fal[copy[s]][s] for s in range(S)]
        carries = copy[ins_site] in carriers
        carries_del = copy[del_site] in del_carriers   # founder LACKS the segment
        samp_alleles.append(alleles); samp_carry.append(carries)
        samp_carry_del.append(carries_del)
        truth.append(copy)
        # forward token list; record the allele-node indices bounding the inv region
        toks = []
        rs = re = None
        for s in range(S):
            toks.append(f"{F(s)}+")
            if carries and s == ins_site:
                toks.append(f"{INS}+")            # insertion right after F(ins_site)
            if (not carries_del) and s == del_site:
                toks.append(f"{DEL}+")            # non-carriers keep the ref segment
            if s in pd_set and not samp_pd_del[i][s]:
                toks.append(f"{PD(s)}+")          # keep block unless privately deleted
            if s == inv_lo: rs = len(toks)        # index of allele(inv_lo)
            toks.append(f"{A(s)}+" if alleles[s] else f"{R(s)}+")
            if s == inv_hi: re = len(toks) - 1    # index of allele(inv_hi)
        toks.append(f"{F(S)}+")
        # inversion: reverse the block [allele(inv_lo) .. allele(inv_hi)] and flip
        # each node's orientation. Entry anchor F(inv_lo)+ and exit F(inv_hi+1)+
        # stay forward => clean breakpoints. (inv region is disjoint from ins_site.)
        if samp_inv[i] and rs is not None:
            block = toks[rs:re + 1]
            toks[rs:re + 1] = [t[:-1] + ('-' if t[-1] == '+' else '+')
                               for t in reversed(block)]
        g.write(f"P\tsmp{i}#1#sim\t" + ",".join(toks) + "\t*\n")

# ---- founder copy matrix (full ground truth, for precision scoring) ------------
# Row i = the founder each sample-haplotype i copies at every site; two haplotypes
# are truly IBD over any run where their rows agree. Bp of each site's filler is
# also written so a caller's bp interval can be tested against founder agreement.
with open(a.out + ".copy", "w") as c:
    c.write("# site_bp: " + " ".join(str(bp[F(s)]) for s in range(S)) + "\n")
    for i in range(N):
        c.write(" ".join(str(x) for x in truth[i]) + "\n")

# ---- genetic map: marker at each filler bp, uniform cM ----
rate = a.cm_per_mb / 1e6
with open(a.out + ".map", "w") as m:
    for s in range(S + 1):
        b = bp[F(s)]
        m.write(f"sim\tm{s}\t{b*rate:.8f}\t{b}\n")

# ---- per-pair IBD truth over the insertion locus -------------------------------
# Ground truth = pairs of sample haplotypes that co-inherit the SAME founder across
# ins_site (true IBD by descent). For each such pair we emit the maximal run around
# ins_site where BOTH copy that founder, as a reference bp/cM interval. carrier=1
# means the shared founder carries the INS => the IBD segment spans the insertion
# (the moat: graph keeps it whole, linear hap-ibd should break/under-credit it).
# Only pairs sharing the founder at ins_site are written (the locus-spanning set),
# grouped by founder so cost is Sum_f g_f^2, not N^2.
def ibd_span(ci, cj, f, anchor):
    lo = anchor
    while lo > 0 and ci[lo - 1] == f and cj[lo - 1] == f:
        lo -= 1
    hi = anchor
    while hi + 1 < S and ci[hi + 1] == f and cj[hi + 1] == f:
        hi += 1
    return lo, hi

# group samples by the founder they copy at ins_site
by_founder = {}
for i in range(N):
    by_founder.setdefault(truth[i][ins_site], []).append(i)

n_pairs = 0
with open(a.out + ".truth", "w") as t:
    t.write("# per-pair IBD spanning the insertion locus (founder-IBD ground truth)\n")
    t.write("hap_i\thap_j\tfounder\tcarrier\tstart_bp\tend_bp\tlen_sites\tcM\n")
    for f, members in by_founder.items():
        carrier = 1 if f in carriers else 0
        for ai in range(len(members)):
            for aj in range(ai + 1, len(members)):
                i, j = members[ai], members[aj]
                lo, hi = ibd_span(truth[i], truth[j], f, ins_site)
                sbp = bp[F(lo)]
                ebp = bp[F(hi + 1)]
                cm = (ebp - sbp) * rate
                t.write(f"{i}\t{j}\t{f}\t{carrier}\t{sbp}\t{ebp}\t{hi-lo+1}\t{cm:.6f}\n")
                n_pairs += 1

# ---- per-pair IBD truth over the DELETION locus --------------------------------
# Same construction as the insertion truth, grouped at del_site. carrier=1 means
# the shared founder LACKS the reference segment => both haplotypes skip it (the
# IBD segment spans a deletion junction; graph types deletion=1, hap-ibd sees a
# point marker and loses the deleted length).
by_founder_del = {}
for i in range(N):
    by_founder_del.setdefault(truth[i][del_site], []).append(i)
n_del_pairs = n_del_carrier = 0
with open(a.out + ".del_truth", "w") as t:
    t.write("# per-pair IBD spanning the deletion locus (founder-IBD ground truth)\n")
    t.write("hap_i\thap_j\tfounder\tcarrier\tstart_bp\tend_bp\tlen_sites\tcM\n")
    for f, members in by_founder_del.items():
        carrier = 1 if f in del_carriers else 0
        for ai in range(len(members)):
            for aj in range(ai + 1, len(members)):
                i, j = members[ai], members[aj]
                lo, hi = ibd_span(truth[i], truth[j], f, del_site)
                sbp = bp[F(lo)]
                ebp = bp[F(hi + 1)]
                cm = (ebp - sbp) * rate
                t.write(f"{i}\t{j}\t{f}\t{carrier}\t{sbp}\t{ebp}\t{hi-lo+1}\t{cm:.6f}\n")
                n_del_pairs += 1; n_del_carrier += carrier

# ---- linear VCF export (SNPs + the SV loci as biallelic markers) -----------------
# UPDATED: the insertion AND deletion now appear in the VCF as single biallelic
# records (REF/ALT carry the true SV sequence/length). This is what a real SV-aware
# callset would feed hap-ibd. The point of the native-vs-native comparison is that
# hap-ibd uses ONLY the POS + 0/1 genotype of each marker: it treats the SV as a
# POINT marker indistinguishable from a SNP, discarding the length/sequence/type.
# The IBD it then calls is correct; its representation of the SV is not.
# Diploid pairing: individual k = sample-haplotypes (2k, 2k+1). hap-ibd consumes
# this VCF + <out>.map (already PLINK 4-col). Map sim hap h back from hap-ibd
# output via: h = 2*ind + (hap-1).
n_indiv = N // 2   # drops a trailing unpaired haplotype if N is odd
# linear allele of hap h at ref site s: an inverted carrier has its region
# alleles SCRAMBLED to mirrored reference positions (what a reference-coordinate
# caller sees), so a carrier and a non-carrier of the same founder MISMATCH
# across the region -> hap-ibd breaks. The graph instead matches them
# forward-vs-reverse (inverted=1) and keeps the segment whole.
def vcf_allele(h, s):
    if samp_inv[h] and inv_lo <= s <= inv_hi:
        return samp_alleles[h][mirror(s)]
    return samp_alleles[h][s]

def gt_row(allele_of):  # allele_of(h) -> 0/1 for each haplotype h
    return "\t".join(f"{allele_of(2*k)}|{allele_of(2*k+1)}" for k in range(n_indiv))

# (pos, vcf_line) records; sorted by POS so the VCF is coordinate-ordered.
recs = []
for s in range(S):
    recs.append((bp[R(s)] + 1,
                 f"sim\t{bp[R(s)]+1}\tm{s}\tA\tC\t.\tPASS\t.\tGT\t"
                 + gt_row(lambda h, s=s: vcf_allele(h, s))))
# the INSERTION as ONE biallelic indel marker (REF=anchor, ALT=anchor+inserted seq).
# hap-ibd uses only POS + GT: to it this is a point marker identical to a SNP.
ins_pos = bp[F(ins_site)] + 1
recs.append((ins_pos,
             f"sim\t{ins_pos}\tINS_{ins_site}\tA\tA{'G'*a.ins_bp}\t.\tPASS\t"
             f"SVTYPE=INS;SVLEN={a.ins_bp}\tGT\t"
             + gt_row(lambda h: 1 if samp_carry[h] else 0)))
# the DELETION as ONE biallelic indel marker (REF=anchor+deleted seq, ALT=anchor).
del_pos = bp[DEL]   # anchor base just inside the reference segment
recs.append((del_pos,
             f"sim\t{del_pos}\tDEL_{del_site}\tT{'T'*(a.del_bp-1)}\tT\t.\tPASS\t"
             f"SVTYPE=DEL;SVLEN=-{a.del_bp}\tGT\t"
             + gt_row(lambda h: 1 if samp_carry_del[h] else 0)))
recs.sort(key=lambda r: r[0])
with open(a.out + ".vcf", "w") as v:
    v.write("##fileformat=VCFv4.2\n")
    v.write("##contig=<ID=sim>\n")
    v.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
    v.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
    v.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    v.write(cols + "".join(f"\tind{k}" for k in range(n_indiv)) + "\n")
    for _, line in recs:
        v.write(line + "\n")

# ---- per-pair INVERSION truth over the inversion region ------------------------
# Pairs that co-inherit the same founder across inv_mid (SNP-IBD over the region).
# moat=1 means they DIFFER in inversion status (one inverted, one not): a linear
# method breaks (mirrored alleles mismatch) but the graph matches forward-vs-reverse
# (inverted=1) and keeps the IBD whole. moat=0 = same status (both tools recover).
inv_bp_lo, inv_bp_hi = bp[F(inv_lo)], bp[F(inv_hi + 1)]
by_founder_inv = {}
for i in range(N):
    by_founder_inv.setdefault(truth[i][inv_mid], []).append(i)
n_inv_pairs = n_inv_moat = 0
with open(a.out + ".inv_truth", "w") as t:
    t.write(f"# inversion region sites [{inv_lo},{inv_hi}] ref bp "
            f"[{inv_bp_lo},{inv_bp_hi}]; per-pair SNP-IBD over inv_mid\n")
    t.write("hap_i\thap_j\tfounder\tmoat\tinv_i\tinv_j\tstart_bp\tend_bp\tlen_sites\tcM\n")
    for f, members in by_founder_inv.items():
        for ai in range(len(members)):
            for aj in range(ai + 1, len(members)):
                i, j = members[ai], members[aj]
                lo, hi = ibd_span(truth[i], truth[j], f, inv_mid)
                sbp, ebp = bp[F(lo)], bp[F(hi + 1)]
                cm = (ebp - sbp) * rate
                moat = 1 if samp_inv[i] != samp_inv[j] else 0
                t.write(f"{i}\t{j}\t{f}\t{moat}\t{int(samp_inv[i])}\t{int(samp_inv[j])}"
                        f"\t{sbp}\t{ebp}\t{hi-lo+1}\t{cm:.6f}\n")
                n_inv_pairs += 1; n_inv_moat += moat

n_carry = sum(samp_carry)
n_carrier_pairs = sum(1 for f in by_founder if f in carriers
                      for _ in range(len(by_founder[f]) * (len(by_founder[f]) - 1) // 2))
print(f"wrote {a.out}.gfa : ref + {N} samples, {S} sites, {total_bp} bp ref")
print(f"  insertion {a.ins_bp} bp at site {ins_site} (ref bp {bp[F(ins_site)]}), "
      f"carried by {n_carry}/{N} samples ({a.carrier_founders}/{K} founders)")
print(f"  cM span {total_bp*rate:.3f} cM ; map {a.out}.map")
print(f"  inversion sites [{inv_lo},{inv_hi}] ref bp [{inv_bp_lo},{inv_bp_hi}] "
      f"({(inv_bp_hi-inv_bp_lo)*rate:.3f} cM), carried by {sum(samp_inv)}/{N} samples "
      f"(freq {a.inv_freq})")
print(f"  deletion  {a.del_bp} bp at site {del_site} (ref bp {bp[DEL]}), "
      f"LACKED by {sum(samp_carry_del)}/{N} samples "
      f"({a.del_carrier_founders}/{K} founders)")
if pd_sites:
    _tot = sum(sum(1 for v in d.values() if v) for d in samp_pd_del)
    print(f"  PRIVATE deletions: {len(pd_sites)} blocks x {a.priv_del_bp} bp at sites "
          f"{pd_sites}, rate {a.priv_del_rate} -> {_tot} private-del events "
          f"(NOT in VCF; within-IBD indel discordance for the cM-gap test)")
print(f"  {a.out}.truth     : {n_pairs} INS locus-spanning IBD pairs "
      f"({n_carrier_pairs} carrier/insertion-spanning)")
print(f"  {a.out}.del_truth : {n_del_pairs} DEL locus-spanning IBD pairs "
      f"({n_del_carrier} carrier/deletion-spanning)")
print(f"  {a.out}.inv_truth : {n_inv_pairs} INV region SNP-IBD pairs "
      f"({n_inv_moat} moat = differ in orientation, linear breaks)")
print(f"  {a.out}.vcf       : {n_indiv} diploid individuals (SNPs + INS + DEL "
      f"as biallelic markers, for linear hap-ibd)")
