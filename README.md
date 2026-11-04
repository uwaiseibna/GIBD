# GIBD: structural variant aware IBD detection on pangenome graphs

GIBD finds identity-by-descent (IBD) segments between all pairs of haplotypes in
a pangenome graph and labels each segment with the structural variants (SVs)
that the two haplotypes carry.

GIBD reads the haplotypes as paths in a [GBWT](https://github.com/jltsiren/gbwt)
index. A single traversal of the GBWT enumerates the pairwise-maximal exact
matches (gMEMs) between all pairs of haplotypes and visits each shared subpath
once, however many haplotypes share it. On graphs of bounded degree this takes
O(n + occ) time, where n is the number of node visits on the haplotype paths in
both orientations and occ is the number of gMEMs. GIBD then joins the gMEMs of
each pair into IBD segments with the seed-and-extend model of hap-ibd
(Zhou et al., 2020) and labels each segment with

- **insertions** and **deletions** that both haplotypes carry, relative to a
  reference genome in the graph;
- **inversions** that both haplotypes carry (`INV`), and inversion loci where
  both haplotypes have the reference orientation (`REF`);
- **inversions that only one haplotype carries** (`inverted=1`, `OPP`). The two
  haplotypes traverse the same nodes in opposite orientations there, so GIBD
  joins the matches on both sides into one segment.

## Citation

> Uwaise Islam, Franco Marsico, Erik Garrison, Christina Boucher.
> *Structural Variant Aware Identity-by-Descent Detection on Pangenome Graphs.*
> 2026. (Preprint link to be added.)

```bibtex
@article{islam2026gibd,
  title  = {Structural Variant Aware Identity-by-Descent Detection on Pangenome Graphs},
  author = {Islam, Uwaise and Marsico, Franco and Garrison, Erik and Boucher, Christina},
  year   = {2026},
  note   = {Preprint}
}
```

## Installation

GIBD is written in C++17 and depends on [GBWT](https://github.com/jltsiren/gbwt)
and the [vgteam fork of sdsl-lite](https://github.com/vgteam/sdsl-lite). It has
been tested on Linux with GCC 11.5.

Requirements: `git`, `g++` (C++17, OpenMP), `cmake`, `make`, `libzstd` (with
`pkg-config`), and Python 3 for the tests.

```bash
git clone https://github.com/uwaiseibna/GIBD.git
cd GIBD
scripts/install_deps.sh        # builds sdsl-lite and GBWT into deps/
make                           # builds gibd and gfa2gbwt
make test                      # 21 checks on small graphs; each prints PASS
```

If sdsl-lite and GBWT are already installed:

```bash
make GBWT_DIR=/path/to/gbwt SDSL_PREFIX=/path/to/sdsl/install/prefix
```

| Program    | Purpose |
|------------|---------|
| `gibd`     | Enumerates gMEMs, joins them into IBD and HBD segments, and labels SVs. |
| `gfa2gbwt` | Builds a bidirectional GBWT from the P- and W-lines of a GFA, keeping sample and haplotype names (PanSN). |

## Quick start

Simulate 40 haplotypes with an insertion, a deletion and an inversion, build the
GBWT, and report labeled IBD segments:

```bash
python3 tests/sim_panel.py --sites 2000 --founders 6 --samples 40 --p-switch 0.002 \
    --filler-bp 1000 --cm-per-mb 2.0 --inv-freq 0.3 --seed 1 --out demo
./gfa2gbwt demo.gfa demo                        # -> demo.gbwt
./gibd --gfa demo.gfa --ref-path 0 --map demo.map --min-markers 0 \
    --min-seed-cM 1 --min-extend-cM 0.5 --min-output-cM 1 \
    -o demo.ibd.tsv demo.gbwt
```

## Input

1. **A GBWT** of the haplotypes and the **GFA** of the same graph (`--gfa`), for
   node lengths. For large graphs, a table of node lengths (`--nlen`) is smaller
   and faster to read:
   `awk -F'\t' '$1=="S"{print $2"\t"length($3)}' graph.gfa > graph.nlen`
2. **A reference genome** in the graph, which places each match on the
   chromosome and is the baseline for SV labels: `--ref-prefix S` selects every
   path whose PanSN name starts with `S`, and `--ref-path N` selects path `N`.
3. **A genetic map** in PLINK format (`chrom marker cM bp`) on the coordinates of
   the reference genome (`--map`).

Without a map, `gibd` reports the gMEMs only.

### From a pangenome graph (GBZ)

Extract one chromosome with all of its haplotypes using
[vg](https://github.com/vgteam/vg), for example from the HPRC graph:

```bash
vg chunk -x hprc-v2.0-mc-grch38.gbz --gbz --contig chr22 -b chr22   # -> chr22_0_chr22.gbz
vg convert -f chr22_0_chr22.gbz > chr22.gfa
./gfa2gbwt chr22.gfa chr22                                          # -> chr22.gbwt
./gibd --gfa chr22.gfa --ref-prefix 'GRCh38#0#chr22' --map chr22.map \
    -o chr22.ibd.tsv chr22.gbwt
```

Use the component chunk above. `vg chunk -T` drops haplotypes on large regions.

### From a phased VCF

Build a graph of a phased, bgzipped and indexed VCF with vg. The contig names in
the VCF and the FASTA must agree, and the reference path is named after the
contig:

```bash
vg autoindex --workflow giraffe -r ref.fa -v phased.vcf.gz -p panel   # -> panel.giraffe.gbz
vg convert -f panel.giraffe.gbz > panel.gfa
./gfa2gbwt panel.gfa panel                                            # -> panel.gbwt
./gibd --gfa panel.gfa --ref-prefix chr20 --map chr20.map -o panel.ibd.tsv panel.gbwt
```

SVs with explicit sequences in the VCF become alleles of the graph, so GIBD
labels them like any other SV. `vg autoindex` also writes read-mapping indexes,
which GIBD does not use.

## Parameters

The segment parameters follow hap-ibd and default to its values. Run `./gibd`
without arguments for the full list.

| Option | Default | Meaning |
|--------|---------|---------|
| `--min-seed-cM` | 2.0 | Minimum length of the gMEM that starts a segment. |
| `--min-extend-cM` | 1.0 | Minimum length of a gMEM that extends a segment. |
| `--min-output-cM` | 2.0 | Minimum length of a reported segment. |
| `--max-gap-bp` | 1000 | Maximum gap between consecutive gMEMs of a segment (`--max-gap-cM` for a gap in cM). |
| `--min-markers` | 100 | Variant sites that a starting gMEM must span; an extending gMEM needs the same count scaled by `min-extend / min-seed`. A variant site is a reference node traversed by some but not all haplotypes. 0 turns it off. |
| `--min-mac` | 2 | Nodes and edges traversed by fewer haplotypes are masked before matching. SV labels still read the original paths. |
| `--min-insertion-bp`, `--min-deletion-bp`, `--min-inversion-bp` | 50 | Minimum SV length for a label. |
| `--sv-length-tolerance` F, `--sv-length-tolerance-bp` B | 0.10, 10 | Two haplotypes share an SV of the same type when their reference intervals overlap and their lengths differ by at most max(B bp, F × the longer length). |
| `--min-bp` | 1000 | Minimum gMEM length in bases when only gMEMs are reported. |
| `--no-inversion` | off | Ignore inverted gMEMs, so a segment cannot extend across an inversion. |
| `-o` | stdout | Output file. |

## Output

With a reference genome and a map, `gibd` writes one IBD segment per line. The
first eight columns follow hap-ibd.

| Column | Meaning |
|--------|---------|
| `sample_a`, `hap_a`, `sample_b`, `hap_b` | The two haplotypes. |
| `chrom`, `start_bp`, `end_bp`, `cM` | Position of the segment on the reference genome and its genetic length. |
| `type` | `IBD`, or `HBD` for the two haplotypes of one sample. |
| `inverted` | 1 if only one haplotype carries an inversion inside the segment. |
| `insertion`, `deletion` | 1 if both haplotypes carry an insertion or a deletion inside the segment. |
| `inv_allele` | Inversion loci in the segment as `KIND:start-end` on the reference: `INV` both haplotypes inverted, `REF` both in reference orientation, `OPP` opposite orientations; `.` if none. |

Without a map, `gibd` writes one gMEM per line: `hap_a hap_b inverted end_node
len_nodes len_bp type a_start a_end b_start`, followed by `insertion deletion`
when a reference genome is given.

## Testing

`make test` builds small graphs with insertions, deletions, inversions, tandem
and inverted repeats, and cycles in `test_out/`. It checks the gMEMs against a
brute-force enumeration of all occurrence pairs (`tests/expected/`), the SV
labels against the planted variants, and the rules for joining gMEMs into
segments.
