# GIBD: SV-aware identity-by-descent in pangenome graphs

GIBD finds identity-by-descent (IBD) segments between pairs of haplotypes in a
pangenome graph. It labels each segment with the structural variants (SVs) that
both haplotypes carry.

GIBD reads the haplotypes as paths in a
[GBWT](https://github.com/jltsiren/gbwt) index. It enumerates the
pairwise-maximal exact matches (gMEMs) between all pairs of haplotypes in a single
pass over the index, visiting each shared subpath once no matter how many
haplotypes share it. It then joins the gMEMs of each pair into IBD segments with
the seed-and-extend model of hap-ibd (Zhou et al., 2020). Finally, it labels each
segment with:

- **insertions** and **deletions** that both haplotypes carry, relative to a
  reference genome in the graph;
- **inversions** that both haplotypes carry (`INV`), and inversion loci where both
  haplotypes have the reference orientation (`REF`);
- **inverted segments**: segments in which only one of the two haplotypes carries
  an inversion (`inverted=1`, `OPP`). The two haplotypes traverse the same nodes
  in opposite orientations there, so GIBD joins the flanking matches across the
  inversion instead of breaking the segment.

When every node has a bounded number of neighbors, gMEM enumeration runs in
O(n + z) time, where n is the total length of the haplotype paths (both
orientations) and z is the number of matching pairs of occurrences.

## Citation

If you use GIBD, please cite:

> Uwaise Islam, Franco Marsico, Erik Garrison, Christina Boucher.
> *GIBD: Pairwise-Maximal Matching on Pangenome Graphs for IBD Detection.*
> 2026. (Preprint link to be added.)

```bibtex
@article{islam2026gibd,
  title   = {{GIBD}: Pairwise-Maximal Matching on Pangenome Graphs for {IBD} Detection},
  author  = {Islam, Uwaise and Marsico, Franco and Garrison, Erik and Boucher, Christina},
  year    = {2026},
  note    = {Preprint}
}
```

## Installation

GIBD is written in C++17 and depends on [GBWT](https://github.com/jltsiren/gbwt)
and the [vgteam fork of sdsl-lite](https://github.com/vgteam/sdsl-lite).
It has been tested on Linux with GCC 11.5.

Requirements: `git`, `g++` (C++17, OpenMP), `cmake`, `make`, `libzstd`
(with `pkg-config`), and Python 3 for the tests.

```bash
git clone https://github.com/<user>/gibd.git
cd gibd
scripts/install_deps.sh        # builds sdsl-lite and GBWT into deps/ (paper commits)
make                           # builds gibd and gfa2gbwt
make test                      # 21 checks on synthetic graphs; each prints PASS
```

If sdsl-lite and GBWT are already installed elsewhere:

```bash
make GBWT_DIR=/path/to/gbwt SDSL_PREFIX=/path/to/sdsl/install/prefix
```

### Programs

| Program         | Purpose |
|-----------------|---------|
| `gibd`          | Enumerates gMEMs, joins them into IBD/HBD segments and labels SVs. |
| `gfa2gbwt`      | Builds a bidirectional GBWT from the P- and W-lines of a GFA, keeping sample and haplotype names (PanSN). |

## Quick start

This example simulates 40 haplotypes with insertions, deletions and an inversion,
builds the GBWT and reports labeled IBD segments.

```bash
python3 tests/sim_panel.py --sites 2000 --founders 6 --samples 40 --p-switch 0.002 \
    --filler-bp 1000 --cm-per-mb 2.0 --inv-freq 0.3 --seed 1 --out demo
./gfa2gbwt demo.gfa demo                       # -> demo.gbwt
./gibd --gfa demo.gfa --ref-path 0 --map demo.map --min-markers 0 \
    --min-seed-cM 1 --min-extend-cM 0.5 --min-output-cM 1 \
    -o demo.ibd.tsv demo.gbwt
```

## Input

GIBD needs three inputs to report IBD segments:

1. **A GBWT** of the haplotypes (`*.gbwt`), and the **GFA** of the same graph
   (`--gfa`) for node lengths. Node IDs in the two must agree.
2. **A reference genome** in the graph, which places each match on the chromosome
   and is the baseline for SV labels: `--ref-prefix GRCh38#0#chr22` (every path
   whose PanSN name starts with the prefix) or `--ref-path <path id>`.
3. **A genetic map** in PLINK format (`chrom  marker  cM  bp`) on the coordinates of
   the reference genome (`--map`).

Without the reference genome and the map, `gibd` reports only the gMEMs.

### From an HPRC graph (GBZ)

Extract one chromosome with all of its haplotypes using
[vg](https://github.com/vgteam/vg):

```bash
vg chunk -x hprc-v2.0-mc-grch38.gbz --gbz --contig chr22 -b chr22   # -> chr22_0_chr22.gbz
vg gbwt -Z chr22_0_chr22.gbz --metadata                            # check: hundreds of haplotypes
vg convert -f chr22_0_chr22.gbz > chr22.gfa
./gfa2gbwt chr22.gfa chr22                                          # -> chr22.gbwt (keeps sample names)

./gibd --gfa chr22.gfa --ref-prefix 'GRCh38#0#chr22' --map plink.chr22.GRCh38.map \
    -o chr22.ibd.tsv chr22.gbwt
```

`vg chunk -T` loses haplotypes on large regions; use the component chunk above.

## Parameters

The segment parameters follow hap-ibd and default to its values.

| Option | Default | Meaning |
|--------|---------|---------|
| `--min-seed-cM` | 2.0 | Minimum length of the gMEM that starts a segment. |
| `--min-extend-cM` | 1.0 | Minimum length of a gMEM that extends a segment. |
| `--min-output-cM` | 2.0 | Minimum reported segment length. |
| `--max-gap-bp` | 1000 | Maximum gap between consecutive gMEMs of a segment (`--max-gap-cM` for a genetic gap). |
| `--min-markers` | 100 | Variant sites a starting gMEM must span (an extending gMEM: the same count scaled by `min-extend/min-seed`). A variant site is a reference node traversed by some but not all haplotypes. |
| `--min-mac` | 2 | Nodes or edges traversed by fewer haplotypes are masked before matching. |
| `--min-insertion-bp` / `--min-deletion-bp` / `--min-inversion-bp` | 50 | Minimum SV length for a label. |
| `--min-bp` | 1000 | Minimum gMEM length in bases when only gMEMs are reported. |
| `--no-inversion` | off | Drop inverted matches (the view of a linear method). |
| `-o` | stdout | Output file. |

Run `./gibd` without arguments for the full list.

## Output

With a reference genome and a map, `gibd` writes one IBD segment per line (the
first eight columns follow hap-ibd):

| Column | Meaning |
|--------|---------|
| `sample_a`, `hap_a`, `sample_b`, `hap_b` | The two haplotypes. |
| `chrom`, `start_bp`, `end_bp`, `cM` | Segment position on the reference genome and genetic length. |
| `type` | `IBD`, or `HBD` for the two haplotypes of one sample. |
| `inverted` | 1 if only one haplotype carries an inversion inside the segment. |
| `insertion`, `deletion` | 1 if both haplotypes carry an insertion or deletion inside the segment. |
| `inv_allele` | Inversion loci in the segment, as `KIND:start-end` on the reference: `INV` both haplotypes inverted, `REF` both in reference orientation, `OPP` opposite orientations; `.` if none. |

Without a map, `gibd` writes one gMEM per line: `hap_a hap_b inverted end_node
len_nodes len_bp type a_start a_end b_start`, followed by `insertion deletion`
when a reference genome is given.

## Testing

`make test` builds small graphs with insertions, deletions, inversions, tandem and
inverted repeats and cycles, and runs in `test_out/`. It checks the gMEMs against a
brute-force enumeration of all occurrence pairs (`tests/expected/`), the SV labels
against the planted variants, and the rules for joining gMEMs into segments.
