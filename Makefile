# Makefile for GIBD.
#
#   scripts/install_deps.sh     sdsl-lite + GBWT into deps/ (once)
#   make                        builds gibd and gfa2gbwt
#   make test                   synthetic fixtures in test_out/; every check prints PASS
#
# Dependencies installed elsewhere:
#   make GBWT_DIR=/path/to/gbwt SDSL_PREFIX=/path/to/sdsl/prefix

SHELL       := /bin/bash
GBWT_DIR    ?= deps/gbwt
SDSL_PREFIX ?= deps/sdsl

CXX        = g++
CXXFLAGS   = -std=c++17 -O2 -fopenmp -Wall -Wextra
INCLUDES   = -I$(GBWT_DIR)/include -I$(SDSL_PREFIX)/include -Isrc
LDFLAGS    = -L$(GBWT_DIR)/lib -L$(SDSL_PREFIX)/lib -lgbwt -lsdsl -lzstd -fopenmp
VPATH      = src

.PHONY: all test check clean

all: gibd gfa2gbwt

# gibd = a shared header + the enumerator + the CLI driver:
#   common.h       shared types/helpers (SeedOut, RefCoords, PathView, Options, ...)
#   enumerate.cpp  bidirectional gMEM enumerator over the move structures (move_*.h)
#   main.cpp       CLI, IBD segments, SV labels, output
gibd: enumerate.cpp main.cpp common.h mask_rare.h move_lf.h move_phi.h move_locate.h move_balance.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(filter %.cpp,$^) -o $@ $(LDFLAGS)

gfa2gbwt: gfa2gbwt.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS)

make_test_gbwt: make_test_gbwt.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS)

# The checks run inside test_out/. tests/expected/*.gmems.tsv hold the gMEMs of the
# fixtures from a brute-force enumeration of all occurrence pairs.
# Expected columns: hap_a(1) hap_b(2) inverted(3) a_start(4) a_end(5) b_start(6) len_nodes(7) len_bp(8) type(9)
# gibd columns:     hap_a(1) hap_b(2) inverted(3) end_node(4) len_nodes(5) len_bp(6) type(7) [a_start(8) a_end(9) b_start(10)]
test: all make_test_gbwt
	@mkdir -p test_out
	@$(MAKE) --no-print-directory -C test_out -f $(CURDIR)/Makefile check \
	  BIN=$(CURDIR) TESTS=$(CURDIR)/tests EXP=$(CURDIR)/tests/expected

check:
	$(BIN)/make_test_gbwt
	$(BIN)/gibd    --gfa test.gfa --min-bp 600 -o gibd.tsv test.gbwt
	@echo "--- Seed set parity (incl. inversions) ---"
	diff <(awk 'NR>1{print $$1,$$2,$$3,$$7,$$8,$$9}' $(EXP)/test.gmems.tsv | sort) \
	     <(awk 'NR>1{print $$1,$$2,$$3,$$5,$$6,$$7}' gibd.tsv  | sort) \
	  && echo "PASS: seed sets match" || echo "FAIL: seed sets differ"
	@echo "--- Positions (incl. inversions) ---"
	$(BIN)/gibd --gfa test.gfa --min-bp 600 -o gibd_loc.tsv test.gbwt
	diff <(awk 'NR>1{print $$1,$$2,$$3,$$4,$$5,$$6,$$7,$$8}'    $(EXP)/test.gmems.tsv | sort) \
	     <(awk 'NR>1{print $$1,$$2,$$3,$$8,$$9,$$10,$$5,$$6}'   gibd_loc.tsv | sort) \
	  && echo "PASS: positions match" || echo "FAIL: positions differ"
	@echo "--- Insertion labeling (--ref-path 0) ---"
	$(BIN)/gibd --gfa test.gfa --min-bp 600 --ref-path 0 -o gibd_ins.tsv test.gbwt
	@# (5,6) spans off-reference node 11 -> insertion=1; (0,3) is on-reference -> 0.
	@ins56=$$(awk '$$1==5 && $$2==6 {print $$11}' gibd_ins.tsv); \
	 ins03=$$(awk '$$1==0 && $$2==3 {print $$11}' gibd_ins.tsv); \
	 if [ "$$ins56" = "1" ] && [ "$$ins03" = "0" ]; then \
	   echo "PASS: insertion labels correct (5,6)=1 (0,3)=0"; \
	 else \
	   echo "FAIL: insertion labels (5,6)=$$ins56 (0,3)=$$ins03"; \
	 fi
	@echo "--- Error-gap bridging (--map) ---"
	@# Synthetic map: 1 cM per 100 bp over the reference span.
	printf 'chr1\t.\t0.0\t0\nchr1\t.\t30.0\t3000\n' > test_map.tsv
	@# gMEMs chain on the genetic gap; at 1 cM/100 bp, --max-gap-cM 10 == 1000 bp.
	$(BIN)/gibd --gfa test.gfa --min-bp 600 --ref-path 0 --map test_map.tsv --min-markers 0 \
	  --min-seed-cM 2 --min-output-cM 2 --max-gap-cM 10 -o gibd_bridge.tsv test.gbwt
	$(BIN)/gibd --gfa test.gfa --min-bp 600 --ref-path 0 --map test_map.tsv --min-markers 0 \
	  --min-seed-cM 2 --min-output-cM 2 --max-gap-cM 0 --min-mac 1 -o gibd_split.tsv test.gbwt
	$(BIN)/gibd --gfa test.gfa --min-bp 600 --ref-path 0 --map test_map.tsv --min-markers 0 \
	  --min-seed-cM 2 --min-output-cM 2 --max-gap-cM 0 -o gibd_mac.tsv test.gbwt
	@# Merge columns: sample_a(1) hap_a(2) sample_b(3) hap_b(4) chrom(5)
	@#               start_bp(6) end_bp(7) cM(8) type(9) inverted(10) insertion(11)
	@# In the fixture sample id == path id and each sample has one haplotype.
	@b=$$(awk 'NR>1 && $$1==9 && $$3==10' gibd_bridge.tsv | wc -l); \
	 s=$$(awk 'NR>1 && $$1==9 && $$3==10' gibd_split.tsv  | wc -l); \
	 bspan=$$(awk 'NR>1 && $$1==9 && $$3==10 {print $$6"-"$$7}' gibd_bridge.tsv); \
	 if [ "$$b" = "1" ] && [ "$$s" = "2" ]; then \
	   echo "PASS: (9,10) bridges to 1 segment [$$bspan] at gap=1000, splits to 2 at gap=0"; \
	 else \
	   echo "FAIL: (9,10) bridge=$$b (expected 1) split=$$s (expected 2)"; \
	 fi
	@m=$$(awk 'NR>1 && $$1==9 && $$3==10' gibd_mac.tsv | wc -l); \
	 if [ "$$m" = "1" ]; then \
	   echo "PASS: (9,10) private substitution bridged at gap=0 by --min-mac 2 (graph min-mac)"; \
	 else \
	   echo "FAIL: (9,10) --min-mac 2 at gap=0 gave $$m segment(s) (expected 1)"; \
	 fi
	@echo "--- Sub-seed-flank emission (floor lowered to min-extend-cM) ---"
	@# (0,11) right flank 7-8 is 400 bp (< min-bp 600). With --map the
	@# floor drops to min-extend-cM, so it is emitted and bridges to reach bp 1600.
	@n=$$(awk 'NR>1 && $$1==0 && $$3==11' gibd_bridge.tsv | wc -l); \
	 end=$$(awk 'NR>1 && $$1==0 && $$3==11 {print $$7}' gibd_bridge.tsv); \
	 if [ "$$n" = "1" ] && [ "$$end" = "1600" ]; then \
	   echo "PASS: (0,11) merges to 1 segment reaching bp $$end (sub-min-bp flank bridged)"; \
	 else \
	   echo "FAIL: (0,11) segments=$$n end_bp=$$end (expected 1 and 1600)"; \
	 fi
	@echo "--- SV-aware gap labeling ---"
	@# An SV is labeled only when BOTH haplotypes carry it.
	@# (5,6): BOTH carry node 11 => SHARED insertion => insertion=1.
	@# (0,5): only path5 carries node 11 (one-sided) => NOT shared => insertion=0.
	@# (9,10): nodes 13 vs 14 off-ref on BOTH sides, net 0 => SUBSTITUTION => insertion=0.
	@# (0,4): inversion => inverted=1, insertion=0.
	@i56=$$(awk 'NR>1 && $$1==5 && $$3==6 {print $$11}' gibd_bridge.tsv); \
	 i05=$$(awk 'NR>1 && $$1==0 && $$3==5 {print $$11}' gibd_bridge.tsv); \
	 i910=$$(awk 'NR>1 && $$1==9 && $$3==10 {print $$11}' gibd_bridge.tsv); \
	 v04=$$(awk 'NR>1 && $$1==0 && $$3==4 {print $$10}' gibd_bridge.tsv); \
	 i04=$$(awk 'NR>1 && $$1==0 && $$3==4 {print $$11}' gibd_bridge.tsv); \
	 if [ "$$i56" = "1" ] && [ "$$i05" = "0" ] && [ "$$i910" = "0" ] && [ "$$v04" = "1" ] && [ "$$i04" = "0" ]; then \
	   echo "PASS: (5,6) shared insertion=1; (0,5) one-sided insertion=0; (9,10) substitution insertion=0; (0,4) inverted=1 insertion=0"; \
	 else \
	   echo "FAIL: (5,6) ins=$$i56 (exp 1); (0,5) ins=$$i05 (exp 0); (9,10) ins=$$i910 (exp 0); (0,4) inv=$$v04 ins=$$i04 (exp 1,0)"; \
	 fi
	@echo "--- Substitution / inversion must NOT be labeled DELETION ---"
	@# (9,10): nodes 13 vs 14 => each hap carries its OWN sequence across the gap (a
	@# SUBSTITUTION; net reference-skip 0 on both sides), so the bridged [0-1600] segment
	@# must be deletion=0 (and insertion=0). The (0,4) inversion: a reverse traversal is
	@# NOT a skip, so deletion=0 (inverted=1 stands alone).
	@d910=$$(awk 'NR>1 && $$1==9 && $$3==10 {print $$12}' gibd_bridge.tsv); \
	 i910=$$(awk 'NR>1 && $$1==9 && $$3==10 {print $$11}' gibd_bridge.tsv); \
	 d04=$$(awk 'NR>1 && $$1==0 && $$3==4 {print $$12}' gibd_bridge.tsv); \
	 if [ "$$d910" = "0" ] && [ "$$i910" = "0" ] && [ "$$d04" = "0" ]; then \
	   echo "PASS: substitution (9,10) deletion=0 insertion=0; inversion (0,4) deletion=0 (no substitution/inversion mislabeled as a reference skip)"; \
	 else \
	   echo "FAIL: (9,10) del=$$d910 (exp 0) ins=$$i910 (exp 0); (0,4) del=$$d04 (exp 0)"; \
	 fi
	@echo "--- --no-inversion fidelity: inverted block is a DISCORDANCE barrier, not a bridgeable gap ---"
	@# (0,4) share forward flanks [0-600] and [1200-1600] around an inverted block
	@# [600-1200] (6 cM). Inversion-ON: one segment [0-1600] inverted=1. A faithful
	@# linear view (--no-inversion) must BREAK at the inversion even when max-gap-cM
	@# (10) far exceeds the block: the dropped inverted match is kept as a non-emitted
	@# barrier so the merge cannot silently bridge the two forward flanks.
	$(BIN)/gibd --gfa test.gfa --ref-path 0 --map test_map.tsv --min-markers 0 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 10 \
	  -o barrier_on.tsv test.gbwt 2>/dev/null
	$(BIN)/gibd --gfa test.gfa --ref-path 0 --map test_map.tsv --min-markers 0 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 10 \
	  --no-inversion -o barrier_ni.tsv test.gbwt 2>/dev/null
	@on=$$(awk 'NR>1 && $$1==0 && $$3==4' barrier_on.tsv | wc -l); \
	 onspan=$$(awk 'NR>1 && $$1==0 && $$3==4 && $$6==0 && $$7==1600 && $$10==1' barrier_on.tsv | wc -l); \
	 ni=$$(awk 'NR>1 && $$1==0 && $$3==4' barrier_ni.tsv | wc -l); \
	 nispan=$$(awk 'NR>1 && $$1==0 && $$3==4 && $$6==0 && $$7==1600' barrier_ni.tsv | wc -l); \
	 if [ "$$on" = "1" ] && [ "$$onspan" = "1" ] && [ "$$ni" = "2" ] && [ "$$nispan" = "0" ]; then \
	   echo "PASS: inversion-on (0,4) spans [0-1600] inverted=1; --no-inversion breaks into 2 flanks (barrier blocks the 6cM-block bridge under gap=10)"; \
	 else \
	   echo "FAIL: on=$$on (exp 1) onspan=$$onspan (exp 1); no-inv segs=$$ni (exp 2) spanning=$$nispan (exp 0)"; \
	 fi
	@echo "--- Cyclic fixtures: gMEMs == expected brute-force gMEMs (tandem, inverted tandem, mixed) ---"
	@# The expected file lists every occurrence pair; gibd keeps one witness per (subpath, pair, orientation), so compare on that
	@# key, and require every reported witness to be one of the expected occurrence pairs.
	python3 $(TESTS)/make_tandem_fixture.py >/dev/null
	python3 $(TESTS)/make_inv_tandem_fixture.py >/dev/null
	python3 $(TESTS)/make_mixed_fixture.py >/dev/null
	@for f in tandem inv_tandem mixed; do \
	   $(BIN)/gfa2gbwt $$f.gfa $$f 2>/dev/null; \
	   $(BIN)/gibd --gfa $$f.gfa --ref-path 0 --min-bp 1000 $$f.gbwt -o $$f.gibd.tsv 2>/dev/null; \
	   if diff <(awk 'NR>1{print $$1,$$2,$$3,$$7,$$8,$$9}' $(EXP)/$$f.gmems.tsv | sort -u) \
	           <(awk 'NR>1{print $$1,$$2,$$3,$$5,$$6,$$7}' $$f.gibd.tsv  | sort -u) >/dev/null \
	      && [ -z "$$(comm -13 <(awk 'NR>1{print $$1,$$2,$$3,$$4,$$5,$$6,$$7,$$8,$$9}' $(EXP)/$$f.gmems.tsv | sort) \
	                          <(awk 'NR>1{print $$1,$$2,$$3,$$8,$$9,$$10,$$5,$$6,$$7}' $$f.gibd.tsv | sort))" ]; then \
	     echo "PASS: $$f matches the expected gMEMs (key sets equal, every witness is a true occurrence pair)"; \
	   else echo "FAIL: $$f differs from the expected gMEMs"; fi; \
	 done
	@bp3800=$$(awk 'NR>1 && $$6==3800' tandem.gibd.tsv | wc -l); \
	 bp5600=$$(awk 'NR>1 && $$6==5600' tandem.gibd.tsv | wc -l); \
	 inv3600=$$(awk 'NR>1 && $$3==1 && $$6==3600' inv_tandem.gibd.tsv | wc -l); \
	 inv5400=$$(awk 'NR>1 && $$3==1 && $$6==5400' inv_tandem.gibd.tsv | wc -l); \
	 if [ "$$bp3800" -gt 0 ] && [ "$$bp5600" -gt 0 ] && [ "$$inv3600" -gt 0 ] && [ "$$inv5400" -gt 0 ]; then \
	   echo "PASS: copy-number-aware tandem IBD (3800+5600 bp) and inverted-repeat inversions (3600+5400 bp)"; \
	 else \
	   echo "FAIL: tandem bp3800=$$bp3800 bp5600=$$bp5600 inv_tandem inv3600=$$inv3600 inv5400=$$inv5400 (all exp >0)"; \
	 fi
	@if python3 $(TESTS)/verify_mixed.py mixed.gfa mixed.gibd.tsv mixed >/dev/null 2>&1; then \
	   echo "PASS: mixed-orientation cycles sound (every match backed by a genuine shared subpath)"; \
	 else \
	   echo "FAIL: FALSE IBD on mixed-orientation cycles:"; \
	   python3 $(TESTS)/verify_mixed.py mixed.gfa mixed.gibd.tsv mixed; \
	 fi
	@echo "--- Graph IBD: cM-gap chaining recovers IBD fragmented by within-IBD indels ---"
	@# sim panel with PRIVATE per-haplotype deletions (absent from any SNP-only VCF) =
	@# within-IBD-segment indel discordance. No cM-bridging (--max-gap-cM 0) fragments
	@# at the deletions and under-recalls; the DEFAULT cM-gap (indel = cM-negligible)
	@# recovers full recall at perfect precision vs the founder ground truth. Both runs
	@# are the SAME robust graph mode -- only the gap budget differs (a monotone knob).
	python3 $(TESTS)/sim_panel.py --sites 400 --founders 6 --samples 40 --p-switch 0.005 \
	  --filler-bp 1000 --cm-per-mb 2.0 --priv-del-rate 0.2 --priv-del-bp 3000 \
	  --priv-del-sites 12 --inv-freq 0 --seed 3 --out pd_test >/dev/null 2>&1
	$(BIN)/gfa2gbwt pd_test.gfa pd_test 2>/dev/null
	@# hap-ibd seed-and-extend; graph-scaled floors (--min-extend-cM = standard min-piece
	@# floor that suppresses conserved-backbone over-extension; --min-seed-cM = exact seed).
	@# baseline: no genetic-gap bridging (--max-gap-cM 0) -> breaks at each indel.
	$(BIN)/gibd --gfa pd_test.gfa --ref-path 0 --map pd_test.map --min-markers 0 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 0 \
	  pd_test.gbwt -o pd_test.bp.tsv 2>/dev/null
	@# cM-gap bridges the indels (gap budget 0.01 cM > the 3 kb dels = 0.006 cM).
	$(BIN)/gibd --gfa pd_test.gfa --ref-path 0 --map pd_test.map --min-markers 0 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 0.01 \
	  pd_test.gbwt -o pd_test.cm.tsv 2>/dev/null
	@python3 $(TESTS)/pd_check.py pd_test.truth pd_test.copy pd_test.bp.tsv pd_test.cm.tsv 0.5
	@echo "--- Shared inversion allele (inv_allele) and typing across contig strands ---"
	@# tests/make_invallele_fixture.py: nodes 4-6 inverted on two haplotypes (a locus),
	@# node 8 inverted on one (no locus); inv_rc / ins_rc are stored reverse-complemented.
	python3 $(TESTS)/make_invallele_fixture.py
	$(BIN)/gfa2gbwt invallele.gfa invallele >/dev/null 2>&1
	$(BIN)/gibd --gfa invallele.gfa --ref-path 0 --map invallele.map --min-markers 0 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 0.05 \
	  -o invallele.tsv invallele.gbwt 2>/dev/null
	@row() { awk -F'\t' -v a=$$1 -v b=$$2 'NR>1 && (($$1==a && $$3==b) || ($$1==b && $$3==a)) {print $$10, $$11, $$13}' invallele.tsv; }; \
	 r1=$$(row inv_a inv_b); r2=$$(row inv_a inv_rc); r3=$$(row fwd_a fwd_b); \
	 r4=$$(row inv_a fwd_a); r5=$$(row inv_rc fwd_a); r6=$$(row ins_f ins_rc); r7=$$(row idup_a idup_b); \
	 if [ "$$r1" = "0 0 INV:6000-12000" ] && [ "$$r2" = "0 0 INV:6000-12000" ] \
	    && [ "$$r3" = "0 0 REF:6000-12000" ] && [ "$$r4" = "1 0 OPP:6000-12000" ] && [ "$$r5" = "1 0 OPP:6000-12000" ] \
	    && [ "$${r6#0 1 }" != "$$r6" ] && [ "$${r7#0 1 }" != "$$r7" ] && [ "$${r7#*INV}" = "$$r7" ]; then \
	   echo "PASS: shared inversion allele INV/REF per locus (private inversion no locus), opposite-orientation pairs OPP per locus, typing across contig strands, inverted duplication = insertion not inversion"; \
	 else \
	   echo "FAIL: inv_allele: inv|inv '$$r1' inv|inv_rc '$$r2' fwd|fwd '$$r3' inv|fwd '$$r4' inv_rc|fwd '$$r5' ins|ins_rc '$$r6' idup|idup '$$r7'"; \
	 fi
	@echo "--- Robustness contract: any threshold runs, bounded, invariants hold ---"
	@# DEGENERATE params (extend 0 + negatives) must NOT crash, NOT flood/collapse:
	@# the effective floor is clamped (SAFE_EXTEND_CM), so the run is bounded and valid.
	$(BIN)/gibd --gfa pd_test.gfa --ref-path 0 --map pd_test.map --min-markers 0 \
	  --min-extend-cM 0 --min-seed-cM -1 --min-output-cM 0.5 --max-gap-cM 0.01 \
	  pd_test.gbwt -o pd_test.degen.tsv 2> pd_test.degen.err; \
	  rc=$$?; \
	  if [ "$$rc" != "0" ]; then echo "FAIL: degenerate params crashed (rc=$$rc)"; \
	  else python3 $(TESTS)/invariant_check.py overlap pd_test.degen.tsv >/dev/null \
	         && echo "PASS: degenerate params run + no overlap/containment" \
	         || echo "FAIL: degenerate-param output violates no-overlap invariant"; fi
	@# MONOTONICITY: a looser selection threshold (lower --min-output-cM) must be a
	@# span-superset of a tighter one (no previously-reported IBD dropped).
	$(BIN)/gibd --gfa pd_test.gfa --ref-path 0 --map pd_test.map --min-markers 0 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 0.01 \
	  pd_test.gbwt -o pd_test.loose.tsv 2>/dev/null
	$(BIN)/gibd --gfa pd_test.gfa --ref-path 0 --map pd_test.map --min-markers 0 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 1.0 --max-gap-cM 0.01 \
	  pd_test.gbwt -o pd_test.tight.tsv 2>/dev/null
	@python3 $(TESTS)/invariant_check.py monotone pd_test.loose.tsv pd_test.tight.tsv \
	  | sed 's/^PASS/PASS: min-output monotone (looser superset of tighter) --/'
	@echo "--- Marker floor (hap-ibd min-markers on graph variant sites) ---"
	@# A seed must span >= min-markers variant sites, an extension the scaled count;
	@# raising the floor only removes IBD (span-subset of the floor-free run), and a
	@# floor above the fixture's site count removes every segment.
	$(BIN)/gibd --gfa pd_test.gfa --ref-path 0 --map pd_test.map --min-markers 100 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 0.01 \
	  pd_test.gbwt -o pd_test.mk100.tsv 2> pd_test.mk100.err
	$(BIN)/gibd --gfa pd_test.gfa --ref-path 0 --map pd_test.map --min-markers 100000 \
	  --min-extend-cM 0.02 --min-seed-cM 0.05 --min-output-cM 0.5 --max-gap-cM 0.01 \
	  pd_test.gbwt -o pd_test.mkall.tsv 2>/dev/null
	@n0=$$(($$(wc -l < pd_test.loose.tsv)-1)); n100=$$(($$(wc -l < pd_test.mk100.tsv)-1)); \
	 nall=$$(($$(wc -l < pd_test.mkall.tsv)-1)); sites=$$(grep -o 'Marker sites: [0-9]*' pd_test.mk100.err); \
	 if [ "$$nall" = "0" ] && [ "$$n100" -lt "$$n0" ] && [ "$$n100" -gt "0" ] \
	    && python3 $(TESTS)/invariant_check.py monotone pd_test.loose.tsv pd_test.mk100.tsv >/dev/null; \
	 then echo "PASS: marker floor ($$sites) removes IBD monotonically: $$n0 -> $$n100 rows at 100, 0 above the site count"; \
	 else echo "FAIL: marker floor (rows $$n0 / $$n100 / $$nall)"; fi
	@echo "--- Short inversion inside long IBD (inverted-match floor below the forward floor) ---"
	@# Forward floor 10 kb > the 6 kb inversion: carrier x non-carrier must still be ONE
	@# segment 0-60000 crossing the inversion (inverted=1, OPP), not two flanks.
	python3 $(TESTS)/make_shortinv_fixture.py && $(BIN)/gfa2gbwt shortinv.gfa shortinv >/dev/null 2>&1
	$(BIN)/gibd --gfa shortinv.gfa --ref-path 0 --map shortinv.map --min-markers 0 \
	  --min-extend-cM 10 --min-seed-cM 10 --min-output-cM 20 shortinv.gbwt -o shortinv.tsv 2>/dev/null
	@got=$$(awk -F'\t' 'NR>1 && $$1!="ref" && $$3!="ref" {print $$1"x"$$3":"$$6"-"$$7":"$$10":"$$13}' shortinv.tsv | sort | tr '\n' ' '); \
	 want="fwd_axfwd_b:0-60000:0:REF:26000-32000 inv_axfwd_a:0-60000:1:OPP:26000-32000 inv_axfwd_b:0-60000:1:OPP:26000-32000 inv_axinv_b:0-60000:0:INV:26000-32000 inv_bxfwd_a:0-60000:1:OPP:26000-32000 inv_bxfwd_b:0-60000:1:OPP:26000-32000 "; \
	 [ "$$got" = "$$want" ] && echo "PASS: short inversion (6 kb < 10 kb floor) crossed: carrier x non-carrier = one segment, inverted=1, OPP" \
	   || echo "FAIL: short inversion: got [$$got]"
	@echo "--- Orientation-aware IBS runs (seed/extension on forward+inverted+forward identical up to orientation) ---"
	python3 $(TESTS)/make_orirun_fixture.py && $(BIN)/gfa2gbwt orirun.gfa orirun >/dev/null 2>&1
	$(BIN)/gibd --gfa orirun.gfa --ref-path 0 --map orirun.map --min-markers 0 --min-extend-cM 10 --min-seed-cM 10 \
	  --min-output-cM 20 orirun.gbwt -o orirun.A.tsv 2>/dev/null
	$(BIN)/gibd --gfa orirun.gfa --ref-path 0 --map orirun.map --min-markers 0 --min-extend-cM 10 --min-seed-cM 30 \
	  --min-output-cM 20 orirun.gbwt -o orirun.B.tsv 2>/dev/null
	$(BIN)/gibd --gfa orirun.gfa --ref-path 0 --map orirun.map --min-markers 0 --min-extend-cM 12 --min-seed-cM 12 \
	  --min-output-cM 20 orirun.gbwt -o orirun.C.tsv 2>/dev/null
	@q(){ awk -F'\t' -v q=$$2 'NR>1 && $$1=="inv_a" && $$3==q {printf "%s-%s:%s:%s ", $$6,$$7,$$10,$$13}' $$1; }; \
	 a=$$(q orirun.A.tsv n12); b=$$(q orirun.B.tsv fwd_a); c=$$(q orirun.C.tsv n12_18); \
	 [ "$$a" = "0-22000:0:. 24000-60000:1:OPP:26000-32000 " ] && [ "$$b" = "0-60000:1:OPP:26000-32000 " ] \
	   && [ "$$c" = "0-22000:0:. 36000-60000:0:. " ] \
	   && echo "PASS: orientation-aware IBS runs: short flank stub joins the run, run seeds without a flank seed, run < min-extend is not an extension" \
	   || echo "FAIL: orientation-aware IBS runs: A [$$a] B [$$b] C [$$c]"

clean:
	rm -f gibd gfa2gbwt make_test_gbwt
	rm -rf test_out
