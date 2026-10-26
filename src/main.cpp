/*
 * Command-line driver of GIBD.
 *
 * Loads the GBWT, the node lengths, the reference genome and the genetic map,
 * enumerates the gMEMs (enumerate.cpp), and then either writes the gMEMs or
 * joins them into IBD segments labeled with the SVs both haplotypes carry.
 */
#include "common.h"
#include "mask_rare.h"

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options] <haplotypes.gbwt>\n"
        << "\n"
        << "Input:\n"
        << "  --gfa FILE            GFA of the graph; node lengths are read from its S-lines\n"
        << "  --nlen FILE           node lengths as node_id<TAB>length, instead of --gfa\n"
        << "  --ref-prefix S        reference genome: every path whose PanSN name\n"
        << "                        (sample#hap#contig) starts with S, e.g. 'GRCh38#0#chr22'\n"
        << "  --ref-path N          reference genome: the path with id N\n"
        << "  --ref-start-bp N      chromosome position of the first base of a single\n"
        << "                        reference path [path metadata offset, else 0]\n"
        << "  --map FILE            PLINK genetic map (chrom marker cM bp) on the coordinates\n"
        << "                        of the reference genome; needs --ref-prefix or --ref-path\n"
        << "\n"
        << "IBD segments (with --map; the defaults are those of hap-ibd):\n"
        << "  --min-seed-cM F       minimum length of the gMEM that starts a segment [2.0]\n"
        << "  --min-extend-cM F     minimum length of a gMEM that extends a segment [1.0];\n"
        << "                        values below 0.005 are raised to 0.005\n"
        << "  --min-output-cM F     minimum length of a reported segment [2.0]\n"
        << "  --max-gap-bp N        maximum gap between consecutive gMEMs of a segment,\n"
        << "                        also around inverted gMEMs [1000]\n"
        << "  --max-gap-cM F        maximum gap in cM instead of --max-gap-bp [off]\n"
        << "  --min-markers N       variant sites a starting gMEM must span; an extending gMEM\n"
        << "                        needs N * min-extend / min-seed. A variant site is a\n"
        << "                        reference node walked by some but not all haplotypes.\n"
        << "                        0 = off [100]\n"
        << "  --min-mac N           nodes and edges walked by fewer than N haplotypes are\n"
        << "                        replaced by the most common route before matching;\n"
        << "                        SV labels still read the original walks. 1 = off [2]\n"
        << "\n"
        << "SV labels (with a reference genome):\n"
        << "  --min-insertion-bp N  minimum insertion length [50]\n"
        << "  --min-deletion-bp N   minimum deletion length [50]\n"
        << "  --min-inversion-bp N  minimum inversion length [50]\n"
        << "\n"
        << "gMEMs (without --map):\n"
        << "  --min-bp N            minimum gMEM length in bp [1000]\n"
        << "\n"
        << "Other:\n"
        << "  --no-inversion        ignore inverted gMEMs; a segment cannot be extended\n"
        << "                        across one\n"
        << "  -o FILE               output file [stdout]\n"
        << "\n"
        << "Output (TSV with a header line):\n"
        << "  with --map     sample_a hap_a sample_b hap_b chrom start_bp end_bp cM type\n"
        << "                 inverted insertion deletion inv_allele\n"
        << "  without --map  hap_a hap_b inverted end_node len_nodes len_bp type a_start\n"
        << "                 a_end b_start [insertion deletion]\n";
}

static Options parseArgs(int argc, char* argv[]) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto nxt = [&]() -> std::string {
            if (++i >= argc) throw std::runtime_error("Missing value for " + a);
            return argv[i];
        };
        if      (a == "--gfa")      o.gfa_file  = nxt();
        else if (a == "--nlen")     o.nlen_file = nxt();
        else if (a == "--map")          o.map_file     = nxt();
        else if (a == "--ref-start-bp") { o.ref_start_bp = std::stoll(nxt());
                                          o.ref_start_bp_set = true; }
        else if (a == "--min-seed-cM")   o.min_seed_cM = std::stod(nxt());
        else if (a == "--min-extend-cM") o.min_extend_cM = std::stod(nxt());
        else if (a == "--min-output-cM") o.min_output_cM = std::stod(nxt());
        else if (a == "--max-gap-bp")    o.max_gap_bp    = std::stoll(nxt());
        else if (a == "--max-gap-cM")  { o.max_gap_cM    = std::stod(nxt()); o.max_gap_cM_set = true; }
        else if (a == "--min-markers")   o.min_markers   = static_cast<uint32_t>(std::stoul(nxt()));
        else if (a == "--min-mac")       o.min_mac       = static_cast<uint32_t>(std::stoul(nxt()));
        else if (a == "--min-bp")   o.min_bp    = std::stoull(nxt());
        else if (a == "--ref-path") o.ref_path  = std::stoll(nxt());
        else if (a == "--ref-prefix") o.ref_prefix = nxt();
        else if (a == "--no-inversion") o.no_inversion = true;
        else if (a == "--min-insertion-bp") o.min_ins_bp = std::stoull(nxt());
        else if (a == "--min-deletion-bp") o.min_del_bp = std::stoull(nxt());
        else if (a == "--min-inversion-bp") o.min_inv_bp = std::stoull(nxt());
        else if (a == "-o")         o.output    = nxt();
        else if (a[0] != '-')       o.gbwt_file = a;
        else { std::cerr << "Unknown: " << a << "\n"; std::exit(1); }
    }
    if (o.gbwt_file.empty()) { usage(argv[0]); std::exit(1); }
    if (o.ref_path >= 0 && !o.ref_prefix.empty())
        throw std::runtime_error("--ref-path and --ref-prefix are mutually exclusive");
    const bool refp = (o.ref_path >= 0) || !o.ref_prefix.empty();
    if (!o.map_file.empty()) {
        if (!refp)
            throw std::runtime_error("--map needs a reference genome (--ref-prefix or --ref-path)");
        // Threshold knobs accept any value: clamp a negative length/gap to 0 ("no floor")
        // with a note rather than aborting. Only structural conflicts stay hard errors.
        auto clamp0 = [](double& v, const char* name) {
            if (v < 0.0) { std::cerr << "Note: " << name << " < 0 clamped to 0.\n"; v = 0.0; }
        };
        clamp0(o.min_seed_cM,   "--min-seed-cM");
        clamp0(o.min_extend_cM, "--min-extend-cM");
        clamp0(o.min_output_cM, "--min-output-cM");
        clamp0(o.max_gap_cM,    "--max-gap-cM");
        if (o.max_gap_bp < 0) { std::cerr << "Note: --max-gap-bp < 0 clamped to 0.\n"; o.max_gap_bp = 0; }
    }
    return o;
}

int main(int argc, char* argv[]) {
    const Options opt = parseArgs(argc, argv);

    gbwt::GBWT index;
    {
        std::ifstream in(opt.gbwt_file, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open: " + opt.gbwt_file);
        index.load(in);
    }
    if (index.empty()) { std::cerr << "Empty GBWT\n"; return 1; }
    const bool bidir = index.bidirectional();
    const size_t_ n_seqs  = index.sequences();
    const size_t_ n_paths = bidir ? n_seqs / 2 : n_seqs;
    std::cerr << "GBWT: " << n_seqs << " sequences, " << n_paths
              << " paths, bidirectional=" << bidir << "\n";
    if (n_seqs > std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("GBWT has >= 2^32 sequences; packed occurrence keys would overflow");

    // Sample of each path, for IBD/HBD classification and the output labels.
    const bool has_meta = index.hasMetadata() && index.metadata.hasPathNames();
    auto sampleOf = [&](size_t_ path_id) -> size_t_ {
        if (has_meta && path_id < index.metadata.paths())
            return index.metadata.path(path_id).sample;
        return path_id / 2;
    };
    // Hot pair classification reads this dense array instead of metadata.
    std::vector<uint32_t> sample_by_path(n_paths);
    for (size_t_ p = 0; p < n_paths; ++p)
        sample_by_path[p] = static_cast<uint32_t>(sampleOf(p));
    auto sampleName = [&](size_t_ path_id) -> std::string {
        const size_t_ sid = sampleOf(path_id);
        if (has_meta && index.metadata.hasSampleNames()
            && sid < index.metadata.samples())
            return index.metadata.sample(sid);
        return std::to_string(sid);
    };
    auto hapIndex = [&](size_t_ path_id) -> unsigned {
        if (has_meta && path_id < index.metadata.paths())
            return static_cast<unsigned>(index.metadata.path(path_id).phase) + 1;
        return static_cast<unsigned>(path_id % 2) + 1;
    };

    // Node-id range of the GBWT alphabet: per-node tables are dense over it.
    const size_t_ id_lo = gbwt::Node::id(index.firstNode());
    const size_t_ id_hi = gbwt::Node::id(index.sigma() - 1);
    NodeLens nlen;
    nlen.init(id_lo, id_hi);
    if (!opt.gfa_file.empty()) {
        nodelenFromGFA(opt.gfa_file, nlen);
        std::cerr << "Node lengths from GFA: " << nlen.size() << "\n";
    } else if (!opt.nlen_file.empty()) {
        nodelenFromTSV(opt.nlen_file, nlen);
        std::cerr << "Node lengths from TSV: " << nlen.size() << "\n";
    } else {
        nlen.finalize();
        std::cerr << "No node lengths; using 1 bp/node\n";
    }
    std::unique_ptr<gbwt::FastLocate> fl;   // built from the move structure below

    std::ofstream fout;
    std::ostream* out = &std::cout;
    if (!opt.output.empty() && opt.output != "-") {
        fout.open(opt.output);
        if (!fout) throw std::runtime_error("Cannot open: " + opt.output);
        out = &fout;
    }

    const bool ref_active = (opt.ref_path >= 0) || !opt.ref_prefix.empty();

    if (opt.map_file.empty()) {
        *out << "hap_a\thap_b\tinverted\tend_node\tlen_nodes\tlen_bp\ttype";
        *out << "\ta_start\ta_end\tb_start";
        if (ref_active)  *out << "\tinsertion\tdeletion";
        *out << '\n';
    }

    GeneticMap gmap;
    if (!opt.map_file.empty()) gmap.load(opt.map_file);

    static const double   SAFE_EXTEND_CM = 0.005;
    static const uint64_t SAFE_FLOOR_BP  = 100;
    const double eff_extend_cM = std::max(opt.min_extend_cM, SAFE_EXTEND_CM);
    uint64_t gmem_floor = opt.min_bp;
    if (gmap.loaded()) {
        const double rate = gmap.maxRatePerBp();
        gmem_floor = (rate > 0.0)
            ? static_cast<uint64_t>(std::floor(eff_extend_cM / rate))
            : SAFE_FLOOR_BP;
        if (gmem_floor < SAFE_FLOOR_BP) gmem_floor = SAFE_FLOOR_BP;
    }
    // Inverted matches are found below the forward floor: an inversion that one
    // haplotype of an IBD pair carries is a short inverted match between two
    // long forward ones, and it must be seen for the chain to cross it. One
    // shorter than max-gap-bp is not needed (the gap rule bridges it).
    uint64_t inv_floor = gmem_floor;
    if (gmap.loaded() && !opt.no_inversion)
        inv_floor = std::min<uint64_t>(gmem_floor,
                                       std::max<uint64_t>(SAFE_FLOOR_BP, static_cast<uint64_t>(opt.max_gap_bp)));

    RefCoords  refco;
    std::vector<char> is_ref_path(n_paths, 0);
    if (ref_active) {
        refco.init(id_lo, id_hi);
        std::vector<size_t_> ref_paths;
        if (!opt.ref_prefix.empty()) {
            if (!has_meta)
                throw std::runtime_error("--ref-prefix needs path-name metadata");
            for (size_t_ p = 0; p < index.metadata.paths(); ++p) {
                const gbwt::FullPathName fp = index.metadata.fullPath(p);
                const std::string name = fp.sample_name + "#" +
                    std::to_string(fp.haplotype) + "#" + fp.contig_name;
                if (name.rfind(opt.ref_prefix, 0) == 0) ref_paths.push_back(p);
            }
            if (ref_paths.empty())
                throw std::runtime_error("--ref-prefix '" + opt.ref_prefix +
                                         "' matched no path");
            std::cerr << "Reference prefix '" << opt.ref_prefix << "': "
                      << ref_paths.size() << " subpath fragment(s)\n";
        } else {
            if (static_cast<size_t_>(opt.ref_path) >= n_paths)
                throw std::runtime_error("--ref-path out of range");
            ref_paths.push_back(static_cast<size_t_>(opt.ref_path));
        }

        // Project each reference fragment onto chromosome bp. An explicit
        // --ref-start-bp applies only to a single reference path; otherwise use
        // metadata offsets when available.
        const bool single = ref_paths.size() == 1;
        for (size_t_ pidx : ref_paths) {
            const auto ref_walk = index.extract(bidir ? 2 * pidx : pidx);
            int64_t ref_off = 0;
            if (opt.ref_start_bp_set && single) {
                ref_off = opt.ref_start_bp;
            } else if (has_meta && pidx < index.metadata.paths()) {
                ref_off = static_cast<int64_t>(index.metadata.fullPath(pidx).offset);
            }
            refco.addWalk(ref_walk, nlen, ref_off);
            is_ref_path[pidx] = 1;
        }
        std::cerr << "Reference: " << refco.size() << " node(s) over "
                  << ref_paths.size() << " path(s), bp [" << refco.firstBp()
                  << "," << refco.lastBp() << "]\n";
        if (refco.revisits() > 0)
            std::cerr << "WARNING: reference path revisits "
                      << refco.revisits() << " node(s) (non-tandem cycle); "
                      << "their coordinates use the first occurrence.\n";
        if (refco.lastBp() > INT32_MAX)
            throw std::runtime_error("reference longer than 2^31 bp: merge pieces store int32 coordinates");
        if (gmap.loaded()) {
            std::cerr << "Reference coords: " << refco.size()
                      << " node(s), bp [" << refco.firstBp() << ","
                      << refco.lastBp() << "], cM ["
                      << gmap.cMAtBp(refco.firstBp()) << ","
                      << gmap.cMAtBp(refco.lastBp()) << "] (span "
                      << (gmap.cMAtBp(refco.lastBp()) - gmap.cMAtBp(refco.firstBp()))
                      << " cM)\n";
        }
    }

    // Without a map, gMEMs are written as the enumerator finds them, with the
    // insertion and deletion flags of their span. With a map, each gMEM that
    // can be part of a segment is kept as a merge piece (see CmFilter).
    const bool segments = gmap.loaded();
    const bool do_typing = !segments && ref_active && refco.size() > 0;
    SeedRowWriter writer;
    writer.out        = out;
    writer.ref_active = ref_active;
    writer.reserve();
    std::vector<SeedOut> barrier_seeds;   // inverted gMEMs dropped by --no-inversion
    SeedBufferSink barrier_sink{&barrier_seeds};

    PathView pv;
    std::vector<gbwt::vector_type> full_paths;
    CmFilter cm_filter{&pv, &refco, &nlen, &gmap, eff_extend_cM, opt.min_seed_cM};
    PieceSink piece_sink{&cm_filter, nullptr};
    if (segments) {
        full_paths.resize(n_paths);
        for (size_t_ p = 0; p < n_paths; ++p)
            full_paths[p] = index.extract(bidir ? 2 * p : p);
        std::cerr << "Extracted " << n_paths
                  << " full path(s)\n";
    }
    // Rare differences are masked before matching (mask_rare.h): matches,
    // pieces and gaps live on the masked paths, SV typing on the original ones.
    std::vector<gbwt::vector_type> orig_paths;
    MaskMap mask;
    const bool masked = segments && opt.min_mac > 1;
    if (masked) {
        std::vector<gbwt::vector_type> mpaths;
        const MaskStats ms = maskRare(index, full_paths, opt.min_mac, id_lo, id_hi, mpaths, mask);
        orig_paths.swap(full_paths);
        full_paths.swap(mpaths);
        index = buildBidirGBWT(full_paths, index.metadata, has_meta);
        std::cerr << "Rare masking (min-mac " << opt.min_mac << "): " << ms.rare_nodes << " rare node(s), "
                  << ms.rare_edges << " rare edge(s), " << ms.stretches << " rare stretch(es), "
                  << ms.rewritten << " rewritten to the common route\n";
    }
    // Marker sites (hap-ibd min-markers) on the masked walks; see MarkerSites.
    MarkerSites msites;
    if (segments && opt.min_markers > 0) {
        const size_t nid = id_hi - id_lo + 1;
        std::vector<uint32_t> cnt(nid, 0), last(nid, UINT32_MAX);
        uint32_t n_h = 0;
        for (size_t_ p = 0; p < n_paths; ++p) {
            if (is_ref_path[p] || full_paths[p].empty()) continue;
            ++n_h;
            for (auto w : full_paths[p]) {
                const size_t i = gbwt::Node::id(w) - id_lo;
                if (last[i] != p) { last[i] = static_cast<uint32_t>(p); ++cnt[i]; }
            }
        }
        std::vector<std::pair<int64_t, size_t>> poly;   // (bp start, node index)
        for (size_t i = 0; i < nid; ++i)
            if (cnt[i] > 0 && cnt[i] < n_h && refco.has(i + id_lo))
                poly.push_back({refco.bpStart(i + id_lo), i});
        std::sort(poly.begin(), poly.end());
        int64_t prev_end = INT64_MIN; uint32_t prev_cnt = 0;
        for (const auto& [bp, i] : poly) {
            if (!(bp == prev_end && cnt[i] == prev_cnt)) msites.bp.push_back(bp);
            prev_end = bp + static_cast<int64_t>(nlen(gbwt::Node::encode(i + id_lo, false)));
            prev_cnt = cnt[i];
        }
        msites.bp.erase(std::unique(msites.bp.begin(), msites.bp.end()), msites.bp.end());
        msites.seed_markers = opt.min_markers;
        msites.extend_markers = opt.min_seed_cM > 0.0
            ? static_cast<uint32_t>(std::floor(std::min(1.0, opt.min_extend_cM / opt.min_seed_cM) * opt.min_markers))
            : opt.min_markers;
        cm_filter.sites = &msites;
        std::cerr << "Marker sites: " << msites.bp.size() << " (seed >= " << msites.seed_markers
                  << ", extension >= " << msites.extend_markers << ")\n";
    }
    // The piece filter reads reference spans from a block summary of the
    // (masked) walks.
    RefBlocks ref_blocks;
    std::vector<std::vector<MergePiece>> kept_buckets;
    std::vector<Walk32> full32, orig32;   // the walks above, held as 32-bit node values
    if (segments) {
        toWalk32(full_paths, full32);
        toWalk32(orig_paths, orig32);
        pv.full = &full32;
        ref_blocks.build(full32, refco, nlen);
        pv.blk = &ref_blocks;
        kept_buckets.resize(kPieceBuckets);
        piece_sink.buckets = &kept_buckets;
    }
    // SV detection along a walk (scanWalk). The net reference gap versus
    // off-reference bridge is tested between consecutive reference anchors,
    // never summed across bubbles. An event is its reference interval [left
    // anchor end, right anchor start] and net length; two haplotypes share it
    // when the intervals touch and the lengths agree within 10% (at least
    // 10 bp). Node identity is not compared: a substitution inside inserted
    // sequence or on a flanking anchor copies the node but leaves the event
    // where it is. Anchors are reference nodes walked in the haplotype's own
    // direction (the orientation in which its whole walk covers the most
    // reference bp, so a contig stored reverse-complemented is read
    // backwards); a reference node walked the other way is bridge sequence,
    // so a balanced inversion nets to zero instead of reading as a deletion
    // of its length. A run of such nodes between two anchors is an inversion
    // allele when it lies between those anchors; its event is the reference
    // interval the run covers.
    struct SvEvent { int64_t lo, hi, len; };
    struct Stretch {
        std::vector<SvEvent> ins, del, inv;
        std::vector<std::pair<int64_t, int64_t>> anchors;   // [start, end) bp, when collected
        size_t  n_anchors = 0;                              // counted always
        int64_t an_lo = INT64_MAX, an_hi = INT64_MIN;       // least start, greatest end
        void clear() {
            ins.clear(); del.clear(); inv.clear(); anchors.clear();
            n_anchors = 0; an_lo = INT64_MAX; an_hi = INT64_MIN;
        }
    };
    // Positions [off, off + n) of a full walk, read in place.
    struct WalkSlice {
        const Walk32* w; size_t off, n;
        node_t operator[](size_t i) const { return (*w)[off + i]; }
        size_t size() const { return n; }
    };
    auto refRevDom = [&](const auto& walk) {
        int64_t fwd_bp = 0, rev_bp = 0;
        for (node_t w : walk)
            if (refco.has(gbwt::Node::id(w)))
                (gbwt::Node::is_reverse(w) ? rev_bp : fwd_bp) += static_cast<int64_t>(nlen(w));
        return rev_bp > fwd_bp;
    };
    auto scanWalk = [&](const auto& walk, bool rev_dom, bool inv_only, bool collect_anchors,
                        Stretch& st) {
        int64_t prev_ref_end = -1, bridge = 0;
        int64_t inv_lo = INT64_MAX, inv_hi = INT64_MIN, inv_bp = 0;
        const size_t n = walk.size();
        for (size_t si = 0; si < n; ++si) {
            const node_t w = rev_dom ? walk[n - 1 - si] : walk[si];
            const size_t_ id = gbwt::Node::id(w);
            const int64_t len = static_cast<int64_t>(nlen(w));
            if (!refco.has(id)) { bridge += len; continue; }
            const int64_t rs = refco.bpStart(id);
            if (gbwt::Node::is_reverse(w) != rev_dom) {
                bridge += len;
                inv_lo = std::min(inv_lo, rs); inv_hi = std::max(inv_hi, rs + len);
                inv_bp += len;
                continue;
            }
            if (prev_ref_end >= 0) {
                // In place only: an inverted copy of distant reference (an
                // inverted duplication) stays bridge sequence, not an allele.
                const int64_t slack = std::max<int64_t>(10, (inv_hi - inv_lo) / 10);
                if (inv_bp >= static_cast<int64_t>(opt.min_inv_bp)
                    && inv_lo >= std::min(prev_ref_end, rs) - slack
                    && inv_hi <= std::max(prev_ref_end, rs) + slack)
                    st.inv.push_back({inv_lo, inv_hi, inv_hi - inv_lo});
                const int64_t ref_gap = rs - prev_ref_end;
                if (!inv_only) {
                    if (ref_gap >= 0
                        && bridge - ref_gap >= static_cast<int64_t>(opt.min_ins_bp))
                        st.ins.push_back({prev_ref_end, rs, bridge - ref_gap});
                    if (ref_gap - bridge >= static_cast<int64_t>(opt.min_del_bp))
                        st.del.push_back({std::min(prev_ref_end, rs), std::max(prev_ref_end, rs),
                                          ref_gap - bridge});
                }
            }
            if (!inv_only) {
                ++st.n_anchors;
                st.an_lo = std::min(st.an_lo, rs); st.an_hi = std::max(st.an_hi, rs + len);
                if (collect_anchors) st.anchors.push_back({rs, rs + len});
            }
            inv_lo = INT64_MAX; inv_hi = INT64_MIN; inv_bp = 0;
            prev_ref_end = rs + len;
            bridge = 0;
        }
    };

    // Each haplotype's own direction, and the panel's inversion loci: the
    // reference intervals that at least two haplotypes walk inverted
    // (overlapping runs are one locus), from the original (unmasked) walks.
    std::vector<int8_t> path_rev;
    std::vector<std::pair<int64_t, int64_t>> inv_loci;
    auto computeLoci = [&](const std::vector<Walk32>& Wv) {
        path_rev.assign(Wv.size(), -1);
        struct Run { int64_t lo, hi; uint32_t p; };
        std::vector<Run> runs;
        Stretch st;
        for (size_t_ p = 0; p < Wv.size(); ++p) {
            const auto& W = Wv[p];
            if (W.empty()) continue;
            path_rev[p] = refRevDom(W) ? 1 : 0;
            st.inv.clear();
            scanWalk(W, path_rev[p] == 1, true, false, st);
            for (const SvEvent& e : st.inv) runs.push_back({e.lo, e.hi, static_cast<uint32_t>(p)});
        }
        std::sort(runs.begin(), runs.end(),
                  [](const Run& x, const Run& y) { return x.lo < y.lo; });
        std::vector<uint32_t> carriers;
        for (size_t i = 0; i < runs.size();) {
            int64_t lo = runs[i].lo, hi = runs[i].hi;
            carriers.clear();
            size_t j = i;
            while (j < runs.size() && runs[j].lo < hi) {
                hi = std::max(hi, runs[j].hi);
                carriers.push_back(runs[j].p);
                ++j;
            }
            std::sort(carriers.begin(), carriers.end());
            const size_t nc = std::unique(carriers.begin(), carriers.end()) - carriers.begin();
            if (nc >= 2) inv_loci.push_back({lo, hi});
            i = j;
        }
        std::cerr << "Inversion loci (>= 2 haplotypes): " << inv_loci.size() << "\n";
    };
    // Orientation-aware IBS (see the piece loop): the emission filter keeps
    // forward matches shorter than min-extend that end at an inversion locus.
    if (segments) {
        computeLoci(masked ? orig32 : full32);
        if (!opt.no_inversion) { cm_filter.loci = &inv_loci; cm_filter.tol = opt.max_gap_bp; }
        // A flank stub can be shorter than the forward floor, so with inversion
        // loci forward matches are found down to the inverted floor as well; the
        // emission filter keeps only those that end at a locus edge.
        if (cm_filter.loci && !inv_loci.empty() && inv_floor < gmem_floor) {
            gmem_floor = inv_floor;
            std::cerr << "Inversion loci present: forward floor lowered to " << gmem_floor
                      << " bp (flank stubs at locus edges)\n";
        }
    }
    // gMEM enumeration. Without a map the gMEMs are written and the run ends.
    MoveLF mv;
    MovePhi mp;
    {
        mv.build(index, MoveLF::ALPHA);   // balanced: every seek takes O(1) steps
        std::cerr << "Move structure: " << mv.rows.size() << " runs\n";
        // FastLocate (r-index) sampled by walking the paths over the move
        // structure instead of with per-step record decoding.
        fl = std::make_unique<gbwt::FastLocate>();
        sdsl::int_vector<0> head_row;   // Move-Phi interval of every run head
        std::vector<uint32_t> seq_len;  // path lengths, for capping the Move-Phi intervals
        buildFastLocate(index, mv, *fl, &head_row, &seq_len);
        std::cerr << "FastLocate (r-index): " << fl->size() << " logical runs\n";
        // Path lengths for position recovery. The endmarker row of each
        // sequence is its path start, whose distance to the path end is the
        // path length; one pass over the endmarker's rows gives them all
        // without walking any path.
        std::vector<uint32_t> dpl(n_paths, 0);
        for (gbwt::size_type x : fl->decompressSA(gbwt::ENDMARKER)) {
            const size_t_ seq = fl->seqId(x);
            if (!bidir || !gbwt::Path::is_reverse(seq))
                dpl[bidir ? gbwt::Path::id(seq) : seq] = static_cast<uint32_t>(fl->seqOffset(x));
        }
        // Move-Phi replaces FastLocate's Elias-Fano predecessor search, whose
        // structures are then no longer needed.
        mp.build(*fl, std::move(head_row), &seq_len);
        sdsl::util::clear(fl->last);
        sdsl::util::clear(fl->last_to_run);
        std::cerr << "Move-Phi structure: " << mp.rows.size() << " intervals\n";
        std::cerr << "Enumerating gMEMs (floor " << gmem_floor << " bp, inverted " << inv_floor << " bp)...\n";
        if (!segments) {
            enumerate_gmems(index, *fl, mv, mp, nlen, refco, do_typing,
                            opt.min_del_bp, opt.min_ins_bp, gmem_floor, inv_floor, opt.no_inversion,
                            sample_by_path, dpl, writer, barrier_sink);
            writer.flush();
            std::cerr << "Done. " << writer.n << " gMEM(s).\n";
            return 0;
        }
        enumerate_gmems(index, *fl, mv, mp, nlen, refco, do_typing,
                        opt.min_del_bp, opt.min_ins_bp, gmem_floor, inv_floor, opt.no_inversion,
                        sample_by_path, dpl, piece_sink, barrier_sink);
        size_t n_kept = 0;
        for (const auto& b : kept_buckets) n_kept += b.size();
        std::cerr << "Enumerator: " << n_kept << " gMEM(s) kept as merge pieces\n";
    }
    // IBD segments: chain the merge pieces of each pair, then type the SVs of
    // every reported segment.
    {
        using Piece = MergePiece;
        // Inverted matches suppressed by --no-inversion become merge barriers:
        // forward flanks cannot be chained across their reference intervals.
        auto barKey = [](uint32_t a, uint32_t b) {
            return (static_cast<uint64_t>(std::min(a, b)) << 32) | std::max(a, b);
        };
        std::unordered_map<uint64_t, std::vector<std::pair<int64_t,int64_t>>> barriers;
        for (const SeedOut& s : barrier_seeds) {
            if (!s.have_pos) continue;
            int64_t sb, eb;
            if (!refBpRange(pv, s, refco, nlen, sb, eb))
                continue;
            barriers[barKey(s.hap_a, s.hap_b)].push_back({sb, eb});
        }

        // The piece sink built the pieces during enumeration, bucketed by pair.
        std::vector<std::vector<Piece>> buckets;
        std::vector<Piece> pieces;
        buckets.swap(kept_buckets);

        // GIBD's IBD model, on the masked walks (rare alleles masked, see
        // mask_rare.h):
        //  (a) a forward exact match >= min-extend is a piece;
        //  (b) consecutive pieces of a pair at most max-gap apart chain (gap
        //      test in the piece loop below);
        //  (c) orientation: a match lying at least half inside a longer match
        //      of the other orientation is a copy (inverted repeat or
        //      duplication inside shared sequence) and is dropped. Inverted
        //      matches that fill a hole of the forward matching (from the forward
        //      piece before it to the one after it, every break <= --max-gap-bp)
        //      are pieces whatever their length when there are at most TWO of
        //      them, i.e. at most one break inside the inverted sequence: two
        //      haplotypes IBD across an inversion one of them carries are
        //      identical there in opposite orientation up to one derived
        //      mutation (the single discordance hap-ibd's max-gap tolerates).
        //      More breaks are identity by state between unrelated haplotypes,
        //      so there only inverted matches >= min-extend count. An inverted match
        //      >= min-extend is a piece anywhere (contigs on opposite strands);
        //  (d) a chain needs a piece >= min-seed and is reported at >= min-output.
        // Discordant bp of a pair over the reference interval [x, y) next to a
        // forward piece p: both walks are read from p's edge away from p (after
        // it when `after`, else before it) while they stay in [x, y) or off the
        // reference, and the two node lists are aligned greedily (identical
        // nodes match; each mismatched stretch costs the longer side's bp).
        // Unaligned within a 64-node window: "infinite".
        auto discordantBp = [&](const Piece& p, int dir_a, int dir_b, bool after,
                                int64_t x, int64_t y) -> int64_t {
            constexpr int64_t kInf = INT64_MAX / 4;
            auto collect = [&](uint32_t hap, uint32_t lo, uint32_t hi, int dir) {
                std::vector<node_t> out;
                const auto& w = (*pv.full)[hap];
                const int step = after ? dir : -dir;
                int64_t i = (step > 0) ? int64_t(hi) + 1 : int64_t(lo) - 1;
                for (; i >= 0 && i < int64_t(w.size()) && out.size() < 200000; i += step) {
                    const size_t_ id = gbwt::Node::id(w[i]);
                    if (refco.has(id)) {
                        const int64_t s = refco.bpStart(id), e = s + int64_t(nlen(w[i]));
                        if (e <= x || s >= y) break;
                    }
                    out.push_back(w[i]);
                }
                return out;
            };
            const auto A = collect(p.hap_a, p.a_lo, p.a_hi, dir_a);
            const auto B = collect(p.hap_b, p.b_lo, p.b_hi, dir_b);
            auto bp = [&](const std::vector<node_t>& v, size_t lo, size_t hi) {
                int64_t t = 0;
                for (size_t i = lo; i < hi; ++i) t += int64_t(nlen(v[i]));
                return t;
            };
            constexpr size_t K = 64;
            int64_t disc = 0;
            size_t i = 0, j = 0;
            while (i < A.size() && j < B.size()) {
                if (A[i] == B[j]) { ++i; ++j; continue; }
                bool synced = false;
                for (size_t d = 1; d <= 2 * K && !synced; ++d)
                    for (size_t di = 0; di <= d; ++di) {
                        const size_t dj = d - di;
                        if (i + di < A.size() && j + dj < B.size() && A[i + di] == B[j + dj]) {
                            disc += std::max(bp(A, i, i + di), bp(B, j, j + dj));
                            i += di; j += dj; synced = true;
                            break;
                        }
                    }
                if (!synced || disc > opt.max_gap_bp) return synced ? disc : kInf;
            }
            return disc + std::max(bp(A, i, A.size()), bp(B, j, B.size()));
        };

        // One bucket's pieces: sorted by pair and start, then (c) per pair.
        auto prepareBucket = [&](std::vector<Piece>& pieces) {
            std::sort(pieces.begin(), pieces.end(),
                      [](const Piece& a, const Piece& b) {
                          if (a.hap_a != b.hap_a) return a.hap_a < b.hap_a;
                          if (a.hap_b != b.hap_b) return a.hap_b < b.hap_b;
                          return a.start_bp < b.start_bp;
                      });
            // Kept pieces are written back into `pieces` (a group never keeps
            // more pieces than it has), each group read from a copy.
            size_t n_kept = 0;
            std::vector<Piece> grp;
            std::vector<const Piece*> fwd, inv;
            for (size_t g = 0; g < pieces.size();) {
                size_t h = g;
                while (h < pieces.size() && pieces[h].hap_a == pieces[g].hap_a
                       && pieces[h].hap_b == pieces[g].hap_b) ++h;
                grp.assign(pieces.begin() + g, pieces.begin() + h);
                auto copyOf = [&](const Piece& q) {
                    const int64_t lq = q.end_bp - q.start_bp;
                    for (size_t j = 0; j < grp.size(); ++j) {
                        const Piece& f = grp[j];
                        if (f.inverted == q.inverted || f.end_bp - f.start_bp <= lq) continue;
                        const int64_t ov = std::min(f.end_bp, q.end_bp) - std::max(f.start_bp, q.start_bp);
                        if (ov > 0 && 2 * ov >= lq) return true;
                    }
                    return false;
                };
                fwd.clear(); inv.clear();
                for (size_t i = 0; i < grp.size(); ++i) {
                    if (copyOf(grp[i])) continue;
                    (grp[i].inverted ? inv : fwd).push_back(&grp[i]);
                }
                for (const Piece* q : fwd) pieces[n_kept++] = *q;
                std::vector<char> use(inv.size(), 0);
                for (size_t k = 0; k < inv.size(); ++k)
                    if (inv[k]->ge_extend) use[k] = 1;
                // Holes of the forward matching: [end of the forward coverage so
                // far, start of the next forward piece]. fwd is in start order.
                // A break between the forward flank and the inverted matches is
                // measured as the sequence on which the pair is discordant there
                // (discordantBp), not as reference span: a short identical run
                // below min-extend between the inversion and a mutation is shared
                // sequence, not part of the gap.
                int64_t cov_end = INT64_MIN;
                const Piece* cov_p = nullptr;
                for (size_t f = 0; f < fwd.size(); ++f) {
                    const int64_t hole_lo = cov_end, hole_hi = fwd[f]->start_bp;
                    if (cov_end != INT64_MIN && hole_hi > hole_lo) {
                        // Inverted pieces overlapping the hole, in start order (inv is).
                        std::vector<size_t> in;
                        for (size_t k = 0; k < inv.size(); ++k)
                            if (inv[k]->end_bp > hole_lo && inv[k]->start_bp < hole_hi) in.push_back(k);
                        const int dir_a = cov_p->a_hi < fwd[f]->a_lo ? 1 : -1;
                        const int dir_b = cov_p->b_hi < fwd[f]->b_lo ? 1 : -1;
                        int64_t reach = hole_lo;
                        for (size_t k : in) {
                            const int64_t gap = inv[k]->start_bp - reach;
                            if (gap > opt.max_gap_bp
                                && (reach != hole_lo
                                    || discordantBp(*cov_p, dir_a, dir_b, true, hole_lo, inv[k]->start_bp)
                                           > opt.max_gap_bp))
                                break;
                            reach = std::max<int64_t>(reach, inv[k]->end_bp);
                        }
                        // At most one break inside the inverted sequence; with
                        // more, only matches >= min-extend (pieces already) count.
                        if (!in.empty() && in.size() <= 2
                            && (hole_hi - reach <= opt.max_gap_bp
                                || discordantBp(*fwd[f], dir_a, dir_b, false, reach, hole_hi)
                                       <= opt.max_gap_bp))
                            for (size_t k : in) use[k] = 1;
                    }
                    if (fwd[f]->end_bp > cov_end) { cov_end = fwd[f]->end_bp; cov_p = fwd[f]; }
                }
                for (size_t k = 0; k < inv.size(); ++k) if (use[k]) pieces[n_kept++] = *inv[k];
                g = h;
            }
            pieces.resize(n_kept);
            std::sort(pieces.begin(), pieces.end(),
                      [](const Piece& a, const Piece& b) {
                          if (a.hap_a != b.hap_a) return a.hap_a < b.hap_a;
                          if (a.hap_b != b.hap_b) return a.hap_b < b.hap_b;
                          return a.start_bp < b.start_bp;
                      });
            // Orientation-aware IBS runs. hap-ibd's seed and extension lengths
            // apply to IBS: an uninterrupted stretch of identical sequence. Two
            // haplotypes IBD across an inversion one of them carries are identical
            // there up to orientation, with no discordance at the breakpoints, so
            // consecutive pieces of opposite orientation that abut (gap <= max-gap-bp)
            // with the inverted one at an inversion locus are ONE IBS run: the run,
            // not each piece, must reach min-extend to be kept and min-seed (with the
            // marker floor) to seed. A discordance (a mutation, an error) still ends a
            // run, as in hap-ibd. Forward flank stubs (< min-extend, kept by the
            // emission filter only at locus edges) not in such a run are dropped.
            if (cm_filter.loci) {
                const int64_t tol = opt.max_gap_bp;
                auto atLocus = [&](const Piece& q) {
                    auto it = std::lower_bound(inv_loci.begin(), inv_loci.end(),
                                               std::make_pair(static_cast<int64_t>(q.end_bp), INT64_MIN));
                    if (it != inv_loci.begin()) --it;
                    for (; it != inv_loci.end() && it->first < q.end_bp; ++it)
                        if (it->second > q.start_bp) return true;
                    return false;
                };
                size_t w = 0;
                for (size_t g = 0; g < pieces.size();) {
                    size_t h = g;
                    while (h < pieces.size() && pieces[h].hap_a == pieces[g].hap_a
                           && pieces[h].hap_b == pieces[g].hap_b) ++h;
                    for (size_t r0 = g; r0 < h;) {
                        size_t r1 = r0 + 1;
                        int64_t re = pieces[r0].end_bp;
                        bool flip = false;
                        while (r1 < h) {
                            const Piece& p = pieces[r1 - 1];
                            const Piece& q = pieces[r1];
                            if (p.inverted == q.inverted || std::llabs(q.start_bp - re) > tol
                                || !atLocus(p.inverted ? p : q)) break;
                            flip = true;
                            re = std::max<int64_t>(re, q.end_bp);
                            ++r1;
                        }
                        if (flip) {
                            const int64_t rs = pieces[r0].start_bp;
                            const double rc = gmap.cMAtBp(re) - gmap.cMAtBp(rs);
                            const uint32_t mk = msites.count(rs, re);
                            const bool ext  = rc >= eff_extend_cM && mk >= msites.extend_markers;
                            const bool seed = rc >= opt.min_seed_cM && mk >= msites.seed_markers;
                            for (size_t k = r0; k < r1; ++k) {
                                if (ext)  pieces[k].ge_extend = true;
                                if (seed) pieces[k].ge_seed = true;
                            }
                        }
                        r0 = r1;
                    }
                    for (size_t k = g; k < h; ++k)
                        if (pieces[k].inverted || pieces[k].ge_extend) pieces[w++] = pieces[k];
                    g = h;
                }
                pieces.resize(w);
            }
        };

        // SV typing reads the original (unmasked) walks.
        PathView opv = pv;
        opv.blk = nullptr;   // the block summary describes the masked walks only
        if (masked) opv.full = &orig32;

        std::vector<node_t> sv_span;
        // lo, hi are masked-walk positions (= original ones when unmasked).
        auto typeStretch = [&](uint32_t hap, uint32_t lo, uint32_t hi, bool collect_anchors,
                               Stretch& st) {
            if (masked) { lo = mask.lo(hap, lo); hi = mask.hi(hap, hi); }
            if (hap < path_rev.size() && path_rev[hap] >= 0) {
                // The walk's direction is known: scan positions [lo, hi] in place.
                const auto& w = (*opv.full)[hap];
                const size_t n = (lo <= hi && lo < w.size())
                               ? std::min<size_t>(hi, w.size() - 1) - lo + 1 : 0;
                scanWalk(WalkSlice{&w, lo, n}, path_rev[hap] == 1, false, collect_anchors, st);
                return;
            }
            sv_span.clear();
            if (lo <= hi) opv.forRange(hap, lo, hi, [&](node_t w) { sv_span.push_back(w); });
            scanWalk(sv_span, refRevDom(sv_span), false, collect_anchors, st);
        };
        auto sameEvent = [](const SvEvent& x, const SvEvent& y) {
            if (x.hi < y.lo || y.hi < x.lo) return false;
            const int64_t d = std::llabs(x.len - y.len);
            return d <= std::max<int64_t>(10, std::max(x.len, y.len) / 10);
        };
        auto shared = [&](const std::vector<SvEvent>& a, const std::vector<SvEvent>& b) {
            for (const SvEvent& x : a)
                for (const SvEvent& y : b)
                    if (sameEvent(x, y)) return true;
            return false;
        };
        // Reference orientation at locus L: anchors flank L on both sides, no
        // inverted run of this haplotype touches L, and at least half of L's
        // reference bp is walked as anchors (not deleted).
        auto refAllele = [&](Stretch& st, const std::pair<int64_t, int64_t>& L,
                             const std::vector<int64_t>& starts,
                             const std::vector<int64_t>& cum) {
            if (starts.empty() || starts.front() > L.first || st.anchors.back().second < L.second)
                return false;
            for (const SvEvent& e : st.inv)
                if (e.lo < L.second && L.first < e.hi) return false;
            const size_t i = std::lower_bound(starts.begin(), starts.end(), L.first) - starts.begin();
            const size_t j = std::lower_bound(starts.begin(), starts.end(), L.second) - starts.begin();
            return 2 * (cum[j] - cum[i]) >= L.second - L.first;
        };
        auto sortAnchors = [](Stretch& st, std::vector<int64_t>& starts, std::vector<int64_t>& cum) {
            // Walk order is usually already reference order; equal anchors are
            // identical, so skipping the sort then changes nothing.
            if (!std::is_sorted(st.anchors.begin(), st.anchors.end()))
                std::sort(st.anchors.begin(), st.anchors.end());
            starts.clear(); cum.assign(1, 0);
            int64_t hi = INT64_MIN;
            for (auto& an : st.anchors) {
                starts.push_back(an.first);
                cum.push_back(cum.back() + (an.second - an.first));
                hi = std::max(hi, an.second);
            }
            if (!st.anchors.empty()) st.anchors.back().second = hi;   // max end, for refAllele
        };

        // Shared-SV typing for one forward stretch of a merged segment. Label
        // only variants carried by both haplotypes; for an inversion locus,
        // either allele (inverted or reference orientation) when both carry
        // it, reported with its reference interval. kind: 'I' both inverted,
        // 'R' both reference orientation, 'O' opposite (one of each).
        struct InvAllele { char kind; int64_t lo, hi; };
        auto kindRank = [](char k) { return k == 'I' ? 0 : (k == 'R' ? 1 : 2); };
        std::vector<int64_t> sa_, ca_, sb_, cb_;
        Stretch A, B;   // reused across calls
        auto typeSegment = [&](uint32_t ha, uint32_t hb,
                               uint32_t a_lo, uint32_t a_hi,
                               uint32_t b_lo, uint32_t b_hi,
                               bool& ins, bool& del, std::vector<InvAllele>& alleles) {
            A.clear(); B.clear();
            typeStretch(ha, a_lo, a_hi, false, A);
            typeStretch(hb, b_lo, b_hi, false, B);
            if (!ins && shared(A.ins, B.ins)) ins = true;
            if (!del && shared(A.del, B.del)) del = true;
            for (const SvEvent& x : A.inv)
                for (const SvEvent& y : B.inv)
                    if (sameEvent(x, y))
                        alleles.push_back({'I', std::min(x.lo, y.lo), std::max(x.hi, y.hi)});
            if (inv_loci.empty() || A.n_anchors == 0 || B.n_anchors == 0) return;
            const int64_t lo = A.an_lo, hi = A.an_hi;
            auto it = std::lower_bound(inv_loci.begin(), inv_loci.end(), std::make_pair(lo, INT64_MIN));
            if (it == inv_loci.end() || it->second > hi) return;
            // The stretch covers an inversion locus: rescan, collecting anchors.
            A.clear(); B.clear();
            typeStretch(ha, a_lo, a_hi, true, A);
            typeStretch(hb, b_lo, b_hi, true, B);
            sortAnchors(A, sa_, ca_);
            sortAnchors(B, sb_, cb_);
            for (; it != inv_loci.end() && it->second <= hi; ++it)
                if (refAllele(A, *it, sa_, ca_) && refAllele(B, *it, sb_, cb_))
                    alleles.push_back({'R', it->first, it->second});
        };
        // Opposite alleles at each inversion locus inside a chain that holds an
        // inverted match: one haplotype walks the locus inverted, the other in
        // reference orientation (the same tests as above, per haplotype). Both
        // walks are read over the whole chain, so the reference flanks of the
        // locus and the inverted run are seen together.
        auto invAt = [](const Stretch& st, const std::pair<int64_t, int64_t>& L) {
            for (const SvEvent& e : st.inv)
                if (e.lo < L.second && L.first < e.hi) return true;
            return false;
        };
        auto oppositeAlleles = [&](uint32_t ha, uint32_t hb,
                                   uint32_t a_lo, uint32_t a_hi, uint32_t b_lo, uint32_t b_hi,
                                   int64_t cs, int64_t ce, std::vector<InvAllele>& alleles) {
            auto it = std::lower_bound(inv_loci.begin(), inv_loci.end(), std::make_pair(cs, INT64_MIN));
            if (it == inv_loci.end() || it->second > ce) return;
            A.clear(); B.clear();
            typeStretch(ha, a_lo, a_hi, true, A);
            typeStretch(hb, b_lo, b_hi, true, B);
            sortAnchors(A, sa_, ca_);
            sortAnchors(B, sb_, cb_);
            for (; it != inv_loci.end() && it->second <= ce; ++it) {
                const bool ia = invAt(A, *it), ib = invAt(B, *it);
                if ((ia && !ib && refAllele(B, *it, sb_, cb_)) || (ib && !ia && refAllele(A, *it, sa_, ca_)))
                    alleles.push_back({'O', it->first, it->second});
            }
        };

        // hap-ibd-compatible columns plus graph-specific annotations.
        const std::string chrom = gmap.chrom().empty() ? "." : gmap.chrom();
        *out << "sample_a\thap_a\tsample_b\thap_b\tchrom\tstart_bp\tend_bp\tcM"
                "\ttype\tinverted\tinsertion\tdeletion\tinv_allele\n";

        // Buffer segments so a final pass can enforce non-overlap per haplotype pair.
        struct OutSeg {
            uint32_t ha, hb; int64_t cs, ce; double cM; SegType type;
            bool inv, ins, del;
            std::vector<InvAllele> alleles;
        };
        std::vector<OutSeg> outsegs;
        // Orientation runs of the current chain, typed only if the chain passes.
        struct Run4 { uint32_t a_lo, a_hi, b_lo, b_hi; };
        std::vector<Run4> runs;
        size_t n_pieces = 0;
        for (auto& bucket : buckets) {
            pieces.swap(bucket);
            std::vector<Piece>().swap(bucket);
            prepareBucket(pieces);
            n_pieces += pieces.size();
            for (size_t i = 0; i < pieces.size();) {
                int64_t  cs = pieces[i].start_bp, ce = pieces[i].end_bp;
                double   ce_cM = gmap.cMAtBp(ce);   // kept equal to cMAtBp(ce)
                runs.clear();
                // Inverted = the orientation flips WITHIN the chain (forward and inverted
                // pieces both present). A chain that is inverted throughout is a strand
                // difference between the two walks (a contig stored reverse-complemented),
                // not evidence of an inversion between the haplotypes.
                bool     has_fwd = !pieces[i].inverted, has_inv = pieces[i].inverted;
                bool     ins = false, del = false;
                std::vector<InvAllele> alleles;
                bool     has_seed = pieces[i].ge_seed;
                const SegType type = pieces[i].type;
                // SV typing runs over each run of same-orientation pieces of the chain:
                // an inverted piece's coordinates on b run against the forward ones,
                // so mixing them would corrupt a whole-chain range. An inverted piece
                // stores b in reverse-walk coordinates; they are flipped to b's own
                // walk, and typing reads each haplotype in its own direction, so a
                // pair whose contigs are stored on opposite strands is typed as well.
                bool     run = false, run_inv = false;
                uint32_t a_lo = 0, a_hi = 0, b_lo = 0, b_hi = 0;
                uint32_t ca_lo = UINT32_MAX, ca_hi = 0, cb_lo = UINT32_MAX, cb_hi = 0;   // whole chain
                auto addPiece = [&](const Piece& q) {
                    uint32_t qb_lo = q.b_lo, qb_hi = q.b_hi;
                    if (q.inverted) {
                        const uint32_t lb = pv.len(q.hap_b);
                        qb_lo = lb - 1 - q.b_hi; qb_hi = lb - 1 - q.b_lo;
                    }
                    ca_lo = std::min(ca_lo, q.a_lo); ca_hi = std::max(ca_hi, q.a_hi);
                    cb_lo = std::min(cb_lo, qb_lo);  cb_hi = std::max(cb_hi, qb_hi);
                    if (run && run_inv != q.inverted) {
                        runs.push_back({a_lo, a_hi, b_lo, b_hi});
                        run = false;
                    }
                    if (!run) { a_lo = q.a_lo; a_hi = q.a_hi; b_lo = qb_lo; b_hi = qb_hi; run = true; run_inv = q.inverted; return; }
                    a_lo = std::min(a_lo, q.a_lo); a_hi = std::max(a_hi, q.a_hi);
                    b_lo = std::min(b_lo, qb_lo); b_hi = std::max(b_hi, qb_hi);
                };
                addPiece(pieces[i]);
                auto gapOk = [&](const Piece& q) {
                    return opt.max_gap_cM_set ? (gmap.cMAtBp(q.start_bp) - ce_cM <= opt.max_gap_cM)
                                              : (q.start_bp - ce <= opt.max_gap_bp);
                };

                const uint64_t bkey = barKey(pieces[i].hap_a, pieces[i].hap_b);
                const auto barIt = barriers.find(bkey);
                auto gapBlocked = [&](int64_t lo, int64_t hi) {
                    if (barIt == barriers.end()) return false;
                    for (const auto& iv : barIt->second)
                        if (iv.second > lo && iv.first < hi) return true;
                    return false;
                };

                size_t k = i + 1;
                while (k < pieces.size()
                       && pieces[k].hap_a == pieces[i].hap_a
                       && pieces[k].hap_b == pieces[i].hap_b
                       && gapOk(pieces[k])
                       && !gapBlocked(ce, pieces[k].start_bp)) {
                    if (pieces[k].end_bp > ce) { ce = pieces[k].end_bp; ce_cM = gmap.cMAtBp(ce); }
                    has_fwd = has_fwd || !pieces[k].inverted;
                    has_inv = has_inv || pieces[k].inverted;
                    addPiece(pieces[k]);
                    has_seed = has_seed || pieces[k].ge_seed;
                    ++k;
                }

                if (run) runs.push_back({a_lo, a_hi, b_lo, b_hi});

                const bool   inv = has_fwd && has_inv;
                const double cM = ce_cM - gmap.cMAtBp(cs);
                const bool pass = has_seed && cM >= opt.min_output_cM;
                // SV labels are reported only for passing chains, so only they are typed.
                if (pass)
                    for (const Run4& r : runs)
                        typeSegment(pieces[i].hap_a, pieces[i].hap_b,
                                    r.a_lo, r.a_hi, r.b_lo, r.b_hi, ins, del, alleles);
                if (pass && has_inv && !inv_loci.empty())
                    oppositeAlleles(pieces[i].hap_a, pieces[i].hap_b,
                                    ca_lo, ca_hi, cb_lo, cb_hi, cs, ce, alleles);
                // One entry per allele: stretches of the same segment can repeat it.
                std::sort(alleles.begin(), alleles.end(), [&](const InvAllele& x, const InvAllele& y) {
                    return x.lo != y.lo ? x.lo < y.lo : kindRank(x.kind) < kindRank(y.kind);
                });
                size_t na = 0;
                for (const InvAllele& x : alleles) {
                    if (na && alleles[na - 1].kind == x.kind && x.lo <= alleles[na - 1].hi)
                        alleles[na - 1].hi = std::max(alleles[na - 1].hi, x.hi);
                    else alleles[na++] = x;
                }
                alleles.resize(na);
                if (pass)
                    outsegs.push_back({pieces[i].hap_a, pieces[i].hap_b, cs, ce, cM,
                                       type, inv, ins, del, std::move(alleles)});
                i = k;
            }
            std::vector<Piece>().swap(pieces);
        }
        // Buckets interleave pairs: restore pair order (chain order within a
        // pair) before the deterministic per-pair sort below.
        std::stable_sort(outsegs.begin(), outsegs.end(), [](const OutSeg& a, const OutSeg& b) {
            return a.ha != b.ha ? a.ha < b.ha : a.hb < b.hb;
        });

        // Hard no-overlap guarantee per haplotype pair. Sort by start and keep
        // the longest segment first when starts tie.
        auto pairKey = [](const OutSeg& s) {
            uint32_t lo = std::min(s.ha, s.hb), hi = std::max(s.ha, s.hb);
            return (static_cast<uint64_t>(lo) << 32) | hi;
        };
        std::sort(outsegs.begin(), outsegs.end(),
                  [&](const OutSeg& a, const OutSeg& b) {
                      uint64_t ka = pairKey(a), kb = pairKey(b);
                      if (ka != kb) return ka < kb;
                      if (a.cs != b.cs) return a.cs < b.cs;
                      return (a.ce - a.cs) > (b.ce - b.cs);
                  });
        size_t n_out = 0, n_dropped = 0;
        size_t g = 0;
        while (g < outsegs.size()) {
            size_t h = g;
            while (h < outsegs.size() && pairKey(outsegs[h]) == pairKey(outsegs[g])) ++h;
            int64_t kept_end = INT64_MIN;
            for (size_t j = g; j < h; ++j) {
                const OutSeg& s = outsegs[j];
                if (s.cs < kept_end) { ++n_dropped; continue; }
                kept_end = s.ce;
                *out << sampleName(s.ha) << '\t' << hapIndex(s.ha) << '\t'
                     << sampleName(s.hb) << '\t' << hapIndex(s.hb) << '\t'
                     << chrom << '\t' << s.cs << '\t' << s.ce << '\t' << s.cM << '\t'
                     << (s.type == SegType::HBD ? "HBD" : "IBD") << '\t'
                     << (s.inv ? 1 : 0) << '\t' << (s.ins ? 1 : 0) << '\t'
                     << (s.del ? 1 : 0) << '\t';
                if (s.alleles.empty()) *out << '.';
                for (size_t q = 0; q < s.alleles.size(); ++q)
                    *out << (q ? "," : "")
                         << (s.alleles[q].kind == 'I' ? "INV:" : (s.alleles[q].kind == 'R' ? "REF:" : "OPP:"))
                         << s.alleles[q].lo << '-' << s.alleles[q].hi;
                *out << '\n';
                ++n_out;
            }
            g = h;
        }
        std::cerr << "Done. " << n_pieces << " piece(s) -> "
                  << n_out << " merged segment(s)";
        if (n_dropped) std::cerr << " (" << n_dropped
                                 << " overlapping/contained dropped)";
        std::cerr << ".\n";
    }
    return 0;
}
