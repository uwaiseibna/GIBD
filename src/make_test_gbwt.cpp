/*
 * make_test_gbwt.cpp  —  build a small synthetic GBWT for the gibd tests
 *
 * Graph: 12 nodes, each 200 bp.  Node 12 carries a self-loop (tandem repeat).
 * 9 paths (one haplotype per sample):
 *   path0 (sample 0): 1,2,3,4,5,6,7,8
 *   path1 (sample 1): 1,2,3,4,5,9,7,8        (diverges at node 6 → takes 9)
 *   path2 (sample 2): 1,2,10,4,5,6,7,8       (diverges at node 3 → takes 10)
 *   path3 (sample 3): 1,2,3,4,5,6,7,8        (same as path0)
 *   path4 (sample 4): 1,2,3,6-,5-,4-,7,8     (block 4-5-6 inverted)
 *   path5 (sample 5): 1,2,3,4,5,6,11,7,8     (node 11 inserted vs reference)
 *   path6 (sample 6): 1,2,3,4,5,6,11,7,8     (same as path5; shares insertion)
 *   path7 (sample 7): 1,2,3,4,5,6,12,12,12,7,8  (node 12 tandem ×3 — self-loop)
 *   path8 (sample 8): 1,2,3,4,5,6,12,12,12,7,8  (same as path7; equal-copy tandem)
 *   path9 (sample 9): 1,2,3,13,5,6,7,8          (node-4 substitution → 13: an "error")
 *   path10(sample10): 1,2,3,14,5,6,7,8          (node-4 substitution → 14: a different "error")
 *
 * Error-gap pair (9,10): they share 1-2-3 (left, bp[0,600]) and 5-6-7-8 (right,
 * bp[800,1600]) but diverge at the node-4 position (13 vs 14, off-reference).
 * On the reference this is a 200 bp discordance. The merge (merge mode, --map)
 * must report ONE segment bp[0,1600] when --max-gap-bp >= 200 and TWO segments
 * when --max-gap-bp = 0. Both flanks are >= min-bp 600, so this case is emitted
 * at the default floor (it validates the merge, not the lowered floor).
 *
 * Sub-seed-flank pair (0,11): path11 substitutes node 6 -> 15, so it shares
 * 1-2-3-4-5 (left, bp[0,1000], 10 cM) and only 7-8 (right, bp[1200,1600],
 * 400 bp = 4 cM) with the reference, gapped 200 bp at node 6. The 400 bp right
 * flank is BELOW the raw min-bp 600 floor, so it is dropped in raw mode but
 * emitted in merge mode (floor lowered to min-extend-cM). It must bridge so
 * (0,11) reports ONE segment reaching bp 1600 — this validates Step 1.
 *
 * At min_bp=600 (≥3 nodes × 200 bp), expected matches include:
 *   IBD (0,3): nodes 1-2-3-4-5-6-7-8         → 1600 bp
 *   inverted (0,4): nodes 4-5-6 (reverse)    → 600 bp
 *   IBD (5,6): nodes 1-2-3-4-5-6-11-7-8      → 1800 bp (spans insertion 11)
 *   IBD (7,8): nodes 1-2-3-4-5-6-12-12-12-7-8 → 2200 bp (spans tandem ×3)
 * tests/expected/test.gmems.tsv lists the full seed set; node 11 is off-reference
 * (--ref-path 0), so the (5,6) match is labeled insertion=1.
 *
 * Output: test.gbwt + test.gfa
 * Build:
 *   g++ -std=c++17 -O2 -I deps/gbwt/include -I $SDSL_PREFIX/include \
 *       -L deps/gbwt/lib -L $SDSL_PREFIX/lib \
 *       make_test_gbwt.cpp -o make_test_gbwt -lgbwt -lsdsl -lzstd
 */

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/gbwt.h>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main()
{
    constexpr int N_PATHS  = 12;
    constexpr int N_NODES  = 15;
    constexpr int NODE_LEN = 200;  // bp per node

    // vector_type is std::vector<short_type> where short_type = uint32_t;
    // Node::encode returns node_type (uint64_t), so cast explicitly.
    using vt = gbwt::vector_type::value_type;
    auto enc = [](int id) -> vt {
        return static_cast<vt>(gbwt::Node::encode(id, false));
    };
    auto encr = [](int id) -> vt {   // reverse-orientation node (for inversions)
        return static_cast<vt>(gbwt::Node::encode(id, true));
    };

    // ---- Define paths ----
    std::vector<gbwt::vector_type> paths = {
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(6),enc(7),enc(8) },   // path0
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(9),enc(7),enc(8) },   // path1
        { enc(1),enc(2),enc(10),enc(4),enc(5),enc(6),enc(7),enc(8) },  // path2
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(6),enc(7),enc(8) },   // path3
        { enc(1),enc(2),enc(3),encr(6),encr(5),encr(4),enc(7),enc(8) },// path4 (inv 4-5-6)
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(6),enc(11),enc(7),enc(8) }, // path5 (+11)
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(6),enc(11),enc(7),enc(8) }, // path6 (+11)
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(6),
          enc(12),enc(12),enc(12),enc(7),enc(8) },                          // path7 (12×3)
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(6),
          enc(12),enc(12),enc(12),enc(7),enc(8) },                          // path8 (12×3)
        { enc(1),enc(2),enc(3),enc(13),enc(5),enc(6),enc(7),enc(8) },       // path9 (node4→13)
        { enc(1),enc(2),enc(3),enc(14),enc(5),enc(6),enc(7),enc(8) },       // path10(node4→14)
        { enc(1),enc(2),enc(3),enc(4),enc(5),enc(15),enc(7),enc(8) },       // path11(node6→15)
    };

    // ---- Build DynamicGBWT (bidirectional) ----
    // For a bidirectional GBWT the text must interleave fwd and rev paths:
    //   fwd0 $ rev0 $ fwd1 $ rev1 $ ...
    // has_both_orientations=true signals the builder that this is already interleaved.
    gbwt::vector_type text;
    for (const auto& p : paths) {
        // Forward path
        text.insert(text.end(), p.begin(), p.end());
        text.push_back(static_cast<vt>(gbwt::ENDMARKER));
        // Reverse complement path (reverse node order, flip each node's orientation)
        gbwt::vector_type rev;
        gbwt::reversePath(p, rev);
        text.insert(text.end(), rev.begin(), rev.end());
        text.push_back(static_cast<vt>(gbwt::ENDMARKER));
    }
    gbwt::DynamicGBWT dynamic;
    dynamic.insert(text, /*has_both_orientations=*/true);

    // ---- Attach metadata: 4 samples, each with 1 haplotype (1 path each) ----
    // In a bidirectional GBWT each path contributes 2 sequences (fwd + rev),
    // but metadata path entries correspond to paths (not sequences).
    for (int i = 0; i < N_PATHS; ++i)
        dynamic.metadata.addPath(
            static_cast<gbwt::size_type>(i),  // sample i
            0,                                  // contig 0
            0,                                  // phase 0
            0);                                 // count 0
    dynamic.metadata.setSamples(N_PATHS);
    dynamic.metadata.setHaplotypes(1);
    dynamic.metadata.setContigs(1);
    dynamic.header.set(gbwt::GBWTHeader::FLAG_METADATA);

    // ---- Compress to GBWT and save ----
    gbwt::GBWT index(dynamic);
    {
        std::ofstream out("test.gbwt", std::ios::binary);
        if (!out) throw std::runtime_error("Cannot write test.gbwt");
        index.serialize(out);
    }
    std::cerr << "Wrote test.gbwt: " << index.sequences() << " sequences, "
              << index.sigma() << " alphabet size\n";

    // ---- Write GFA (S-lines only — no edges needed for node lengths) ----
    // Sequence is N_LEN A's (symbolic; tools only read the length).
    {
        std::ofstream gfa("test.gfa");
        if (!gfa) throw std::runtime_error("Cannot write test.gfa");
        std::string seq(NODE_LEN, 'A');
        for (int n = 1; n <= N_NODES; ++n)
            gfa << "S\t" << n << "\t" << seq << "\n";
    }
    std::cerr << "Wrote test.gfa: " << N_NODES
              << " nodes × " << NODE_LEN << " bp\n";

    return 0;
}
