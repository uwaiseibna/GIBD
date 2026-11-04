/*
 * Shared types of GIBD: gMEM records and the sinks that receive them, node
 * lengths, reference coordinates, the genetic map, the options, and the
 * reference-span helpers used to turn gMEMs into merge pieces.
 */
#pragma once

#include <gbwt/gbwt.h>
#include <gbwt/dynamic_gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwt/internal.h>

#include "move_lf.h"
#include "move_locate.h"
#include "move_phi.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <charconv>
#include <tuple>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using node_t   = gbwt::node_type;
// A haplotype walk held as 32-bit GBWT node values (id * 2 + orientation): half
// the memory of gbwt::vector_type. toWalk32 releases each source walk as it goes.
using Walk32 = std::vector<uint32_t>;
inline void toWalk32(std::vector<gbwt::vector_type>& in, std::vector<Walk32>& out) {
    out.assign(in.size(), {});
    for (size_t p = 0; p < in.size(); ++p) {
        auto& w = in[p];
        out[p].reserve(w.size());
        for (node_t x : w) {
            if (x > UINT32_MAX) throw std::runtime_error("GBWT node value exceeds 32 bits");
            out[p].push_back(static_cast<uint32_t>(x));
        }
        gbwt::vector_type().swap(w);
    }
    std::vector<gbwt::vector_type>().swap(in);
}
using size_t_  = gbwt::size_type;

enum class SegType : uint8_t { IBD, HBD };

// One gMEM: haplotypes (path ids) a < b, orientation, the walk positions of
// the match on a ([a_start, a_end]) and its start on b, and its length.
struct SeedOut {
    uint32_t hap_a, hap_b;
    uint32_t len_nodes;
    uint32_t a_start = 0, a_end = 0, b_start = 0;
    node_t   end_node;
    uint64_t len_bp;
    SegType  type;
    bool     inverted;
    bool     have_pos = false;
    bool     has_insertion = false;
    bool     has_deletion  = false;
};

// Sinks for the enumerator. Its two instantiations (enumerate.cpp) are
// SeedRowWriter for gMEM output and PieceSink for IBD segments; inverted gMEMs
// dropped by --no-inversion go to a SeedBufferSink in both.

// A gMEM in reference coordinates, as the segment builder chains it.
// 36 bytes: billions are held at once on large panels. Reference bp fit in
// int32 (checked when the reference is loaded); the piece's cM length is only
// ever compared with min-seed and min-extend, so the two results are stored.
struct MergePiece {
    uint32_t hap_a, hap_b;
    int32_t  start_bp, end_bp;
    uint32_t a_lo, a_hi;
    uint32_t b_lo, b_hi;
    bool     inverted;
    SegType  type;
    bool     ge_seed, ge_extend;   // cM >= min-seed, cM >= min-extend
};
static_assert(sizeof(MergePiece) == 36, "MergePiece layout");

// Kept pieces go to one of kPieceBuckets buckets by ordered pair, so a pair's
// pieces stay together in emission order and no single buffer has to grow (and
// be copied) to hold them all; buckets are then chained one at a time.
constexpr size_t kPieceBuckets = 64;
inline size_t pieceBucket(uint32_t a, uint32_t b) {
    return (static_cast<uint64_t>(a) * 1000003u + b) % kPieceBuckets;
}

// Stores every gMEM it receives.
struct SeedBufferSink {
    std::vector<SeedOut>* seeds;
    void operator()(const SeedOut& s) const { seeds->push_back(s); }
};

// Keeps the gMEMs that can be part of an IBD segment (CmFilter) as merge
// pieces; the filter has already computed each one's reference span.
struct CmFilter;
struct PieceSink {
    const CmFilter* keep;
    std::vector<std::vector<MergePiece>>* buckets;
    inline void operator()(const SeedOut& s) const;
};

// Streaming TSV writer for gMEM output. It writes rows into a reused buffer
// with to_chars, which keeps formatting out of std::string allocation and lets
// the enumerator emit directly without building a gMEM vector.
struct SeedRowWriter {
    std::ostream* out = nullptr;
    bool      ref_active = false;
    uint64_t  n          = 0;      // rows written

    // Enough room for the longest supported row, so each row needs one capacity
    // check rather than one check per field.
    static constexpr size_t kRowMax  = 256;
    static constexpr size_t kFlushAt = 1u << 17;
    std::vector<char> buf;
    char* p = nullptr;

    void reserve() { buf.resize(kFlushAt + kRowMax); p = buf.data(); }
    void flush() {
        const size_t nb = (size_t)(p - buf.data());
        if (nb) { out->write(buf.data(), (std::streamsize)nb); p = buf.data(); }
    }
    void operator()(const SeedOut& s) {
        char* q = p;
        auto U = [&q](uint64_t x) { q = std::to_chars(q, q + 20, x).ptr; };
        U(s.hap_a); *q++ = '\t'; U(s.hap_b); *q++ = '\t';
        *q++ = (s.inverted ? '1' : '0'); *q++ = '\t';
        U(s.end_node); *q++ = '\t';
        U(s.len_nodes); *q++ = '\t'; U(s.len_bp); *q++ = '\t';
        if (s.type == SegType::HBD) { *q++='H'; *q++='B'; *q++='D'; }
        else                        { *q++='I'; *q++='B'; *q++='D'; }
        if (s.have_pos) {
            *q++ = '\t'; U(s.a_start); *q++ = '\t'; U(s.a_end);
            *q++ = '\t'; U(s.b_start);
        } else {
            *q++='\t'; *q++='.'; *q++='\t'; *q++='.'; *q++='\t'; *q++='.';
        }
        if (ref_active) {
            *q++ = '\t'; *q++ = (s.has_insertion ? '1' : '0');
            *q++ = '\t'; *q++ = (s.has_deletion  ? '1' : '0');
        }
        *q++ = '\n';
        p = q;
        if ((size_t)(p - buf.data()) > kFlushAt) flush();
        ++n;
    }
};

// Node lengths, stored densely over the GBWT's node-id range [lo, hi]. Ids
// outside the range never occur on an indexed path and are ignored; a node with
// no recorded length counts as 1 bp.
class NodeLens {
public:
    void init(size_t_ lo, size_t_ hi) { lo_ = lo; len_.assign(hi - lo + 1, kUnset); n_ = 0; }
    void set(size_t_ id, uint64_t len) {
        if (id < lo_ || id - lo_ >= len_.size()) return;
        uint32_t& slot = len_[id - lo_];
        if (slot == kUnset) ++n_;
        slot = static_cast<uint32_t>(len);
    }
    // Replace unset slots by the 1 bp default once loading is done.
    void finalize() { for (uint32_t& x : len_) if (x == kUnset) x = 1; }
    uint64_t operator()(node_t v) const { return len_[gbwt::Node::id(v) - lo_]; }
    const uint32_t* data() const { return len_.data(); }
    size_t_ lo() const { return lo_; }
    size_t  size() const { return n_; }   // ids given a length
private:
    static constexpr uint32_t kUnset = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> len_;
    size_t_ lo_ = 0;
    size_t  n_  = 0;
};

inline void nodelenFromGFA(const std::string& path, NodeLens& m) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open GFA: " + path);
    std::string line, tag, name, seq;
    while (std::getline(in, line)) {
        if (line.size() < 2 || line[0] != 'S') continue;
        std::istringstream ss(line);
        ss >> tag >> name >> seq;
        if (!ss) continue;
        try {
            size_t_ nid = std::stoull(name);
            if (seq != "*") {
                m.set(nid, seq.size());
            } else {
                std::string field;
                while (ss >> field) {
                    if (field.size() > 5 && field.substr(0,5) == "LN:i:") {
                        m.set(nid, std::stoull(field.substr(5)));
                        break;
                    }
                }
            }
        } catch (...) {}
    }
    m.finalize();
}

inline void nodelenFromTSV(const std::string& path, NodeLens& m) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open TSV: " + path);
    size_t_ nid; uint64_t len;
    while (in >> nid >> len) m.set(nid, len);
    m.finalize();
}

// PLINK-style genetic map: chrom marker cM bp. Positions are interpolated
// linearly between markers and extrapolated at the ends, as hap-ibd does.
class GeneticMap {
public:
    void load(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Cannot open map: " + path);

        std::string line, chrom_seen;
        size_t lno = 0;
        int64_t prev_bp = std::numeric_limits<int64_t>::min();
        double  prev_cM = -std::numeric_limits<double>::infinity();
        while (std::getline(in, line)) {
            ++lno;
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            std::string chrom, marker;
            double cM; int64_t bp;
            if (!(iss >> chrom >> marker >> cM >> bp))
                throw std::runtime_error("Bad map line " + std::to_string(lno));
            if (chrom_seen.empty()) { chrom_seen = chrom; chrom_ = chrom; }
            else if (chrom != chrom_seen)
                throw std::runtime_error("Map has multiple chromosomes (saw '"
                                         + chrom_seen + "' then '" + chrom
                                         + "' at line " + std::to_string(lno) + ")");
            if (bp < prev_bp)
                throw std::runtime_error("Map bp decreases at line "
                                         + std::to_string(lno));
            if (cM < prev_cM)
                throw std::runtime_error("Map cM decreases at line "
                                         + std::to_string(lno));
            prev_bp = bp; prev_cM = cM;
            bp_.push_back(bp);
            cM_.push_back(cM);
        }
        if (bp_.size() < 2) throw std::runtime_error("Map needs >= 2 entries");
        buildBuckets();

        const double total_cM = cM_.back() - cM_.front();
        if (total_cM < 5.0)
            std::cerr << "WARNING: map total cM span is " << total_cM
                      << ", implausibly small for a human autosome; check units "
                         "(cM not Morgans) and chromosome.\n";
        std::cerr << "Loaded genetic map: " << bp_.size() << " entries, bp ["
                  << bp_.front() << "," << bp_.back() << "], cM ["
                  << cM_.front() << "," << cM_.back() << "]\n";
    }

    bool loaded() const { return bp_.size() >= 2; }
    const std::string& chrom() const { return chrom_; }

    // Interpolate between bracketing markers; extrapolate outside the map.
    double cMAtBp(int64_t bp) const {
        if (bp_.empty()) return 0.0;
        if (bp <= bp_.front()) {
            const int64_t dbp = bp_[1] - bp_[0];
            const double sl = dbp > 0 ? (cM_[1] - cM_[0]) / (double)dbp : 0.0;
            return cM_.front() + sl * (double)(bp - bp_.front());
        }
        if (bp >= bp_.back()) {
            const size_t n  = bp_.size();
            const int64_t dbp = bp_[n-1] - bp_[n-2];
            const double sl = dbp > 0 ? (cM_[n-1] - cM_[n-2]) / (double)dbp : 0.0;
            return cM_.back() + sl * (double)(bp - bp_.back());
        }
        // The bucket of bp bounds the bracketing marker, so the search reads
        // only that bucket's markers.
        const size_t k = static_cast<size_t>((bp - bp_.front()) / kBucketBp);
        auto it  = std::upper_bound(bp_.begin() + bkt_lo_[k], bp_.begin() + bkt_hi_[k], bp);
        const size_t hi = it - bp_.begin(), lo = hi - 1;
        const int64_t dbp = bp_[hi] - bp_[lo];
        const double fr = dbp > 0 ? (double)(bp - bp_[lo]) / (double)dbp : 0.0;
        return cM_[lo] + fr * (cM_[hi] - cM_[lo]);
    }

    // Steepest cM-per-bp rate. This converts a cM extension floor into a safe
    // bp floor for the enumeration.
    double maxRatePerBp() const {
        double r = 0.0;
        for (size_t i = 1; i < bp_.size(); ++i) {
            const int64_t dbp = bp_[i] - bp_[i-1];
            if (dbp > 0) r = std::max(r, (cM_[i] - cM_[i-1]) / (double)dbp);
        }
        return r;
    }

private:
    // Bucket k covers bp [front + k*W, front + (k+1)*W). For bp in it, the
    // first marker > bp lies in [first marker >= bucket start, first marker >
    // bucket end], so upper_bound over that index range finds the same marker.
    static constexpr int64_t kBucketBp = 4096;
    void buildBuckets() {
        const size_t nb = static_cast<size_t>((bp_.back() - bp_.front()) / kBucketBp) + 1;
        bkt_lo_.resize(nb); bkt_hi_.resize(nb);
        for (size_t k = 0; k < nb; ++k) {
            const int64_t b0 = bp_.front() + static_cast<int64_t>(k) * kBucketBp;
            bkt_lo_[k] = std::lower_bound(bp_.begin(), bp_.end(), b0) - bp_.begin();
            bkt_hi_[k] = std::upper_bound(bp_.begin(), bp_.end(), b0 + kBucketBp - 1) - bp_.begin();
        }
    }
    std::vector<int64_t> bp_;
    std::vector<double>  cM_;
    std::vector<size_t>  bkt_lo_, bkt_hi_;
    std::string          chrom_;
};

// Reference-coordinate lookup, dense over the same node-id range as NodeLens.
// Reference walks give node ids chromosome bp starts; off-reference nodes have
// no entry (bpStart -1). If a reference path revisits a node, the first
// coordinate is kept and the revisit count is reported.
class RefCoords {
public:
    void init(size_t_ lo, size_t_ hi) { lo_ = lo; bp_.assign(hi - lo + 1, kUnset); }

    // Add one reference fragment starting at chromosome bp `off`. With
    // --ref-prefix this is called once per fragment and keeps the bp span of
    // the union.
    void addWalk(const gbwt::vector_type& walk, const NodeLens& nlen, int64_t off) {
        int64_t bp = off;
        if (!any_) { first_bp_ = off; last_bp_ = off; any_ = true; }
        else       { first_bp_ = std::min(first_bp_, off); }
        for (node_t w : walk) {
            int64_t& slot = bp_[gbwt::Node::id(w) - lo_];
            if (slot != kUnset) ++revisits_;
            else { slot = bp; ++n_; }
            bp += static_cast<int64_t>(nlen(w));
        }
        last_bp_ = std::max(last_bp_, bp);
    }

    bool    has(size_t_ id)     const { return bp_[id - lo_] != kUnset; }
    int64_t bpStart(size_t_ id) const { const int64_t b = bp_[id - lo_]; return b == kUnset ? -1 : b; }
    int64_t firstBp()  const { return first_bp_; }
    int64_t lastBp()   const { return last_bp_;  }
    size_t  size()     const { return n_; }
    size_t  revisits() const { return revisits_; }
    // Raw stored start for the enumerator's hot path: any negative value means
    // off-reference (unset slots hold INT64_MIN).
    int64_t raw(size_t_ id) const { return bp_[id - lo_]; }

private:
    static constexpr int64_t kUnset = std::numeric_limits<int64_t>::min();
    std::vector<int64_t> bp_;
    size_t_ lo_ = 0;
    int64_t first_bp_ = 0, last_bp_ = 0;
    size_t  n_ = 0;
    size_t  revisits_ = 0;
    bool    any_ = false;
};

// Command-line options (parsed in main.cpp).
struct Options {
    std::string gbwt_file;
    std::string gfa_file;
    std::string nlen_file;
    std::string output;
    std::string map_file;
    int64_t     ref_start_bp = 0;
    bool        ref_start_bp_set = false;
    double      min_seed_cM   = 2.0;
    double      min_extend_cM = 1.0;
    double      min_output_cM = 2.0;
    int64_t     max_gap_bp    = 1000;
    double      max_gap_cM    = 0.05;
    bool        max_gap_cM_set = false;   // chain on the cM gap only when --max-gap-cM is given
    uint32_t    min_markers   = 100;   // variant sites in a starting gMEM (an extending one: scaled)
    uint32_t    min_mac       = 2;     // rarer nodes and edges are masked before matching
    uint64_t    min_bp    = 1000;
    int64_t     ref_path  = -1;
    std::string ref_prefix;
    bool        no_inversion = false;
    uint64_t    min_ins_bp = 50;
    uint64_t    min_del_bp = 50;
    uint64_t    min_inv_bp = 50;
    double      sv_len_tol = 0.10;   // shared SV: lengths differ by <= max(sv_len_tol_bp, sv_len_tol * longer)
    int64_t     sv_len_tol_bp = 10;
};

// Reference span summary of full walks: for each block of kB consecutive walk
// positions, the least start and greatest end bp of its on-reference nodes
// (INT64_MAX / INT64_MIN when it has none). O(n / kB) space, O(n) to build.
struct RefBlocks {
    static constexpr uint64_t kB = 64;
    std::vector<std::vector<std::pair<int64_t, int64_t>>> b;
    void build(const std::vector<Walk32>& paths, const RefCoords& refco,
               const NodeLens& nlen) {
        b.assign(paths.size(), {});
        for (size_t p = 0; p < paths.size(); ++p) {
            const auto& w = paths[p];
            auto& bp = b[p];
            bp.assign((w.size() + kB - 1) / kB, {INT64_MAX, INT64_MIN});
            for (size_t i = 0; i < w.size(); ++i) {
                const size_t_ id = gbwt::Node::id(w[i]);
                if (!refco.has(id)) continue;
                const int64_t s = refco.bpStart(id), e = s + static_cast<int64_t>(nlen(w[i]));
                auto& x = bp[i / kB];
                x.first = std::min(x.first, s); x.second = std::max(x.second, e);
            }
        }
    }
};

// Read access to the full haplotype walks (and their block summary).
struct PathView {
    const std::vector<Walk32>* full = nullptr;
    const RefBlocks*           blk  = nullptr;   // summary of *full

    uint32_t len(size_t_ h) const { return static_cast<uint32_t>((*full)[h].size()); }
    template <class F>
    void forRange(size_t_ h, uint32_t lo, uint32_t hi, F f) const {
        const auto& w = (*full)[h];
        for (uint32_t i = lo; i <= hi && i < w.size(); ++i) f(w[i]);
    }
};

// Reference bp span of walk positions [a_start, a_end] of `hap`: the least
// start and greatest end of its on-reference nodes, reading whole blocks from
// the summary and only the partial blocks at the two ends node by node.
inline bool refBpRange(
        const PathView& pv, size_t_ hap,
        uint32_t a_start, uint32_t a_end,
        const RefCoords& refco, const NodeLens& nlen,
        int64_t& start_bp, int64_t& end_bp)
{
    const auto& w  = (*pv.full)[hap];
    const auto& bk = pv.blk->b[hap];
    int64_t lo = INT64_MAX, hi = INT64_MIN;
    auto node = [&](uint64_t i) {
        const size_t_ id = gbwt::Node::id(w[i]);
        if (!refco.has(id)) return;
        const int64_t s = refco.bpStart(id);
        lo = std::min(lo, s); hi = std::max(hi, s + static_cast<int64_t>(nlen(w[i])));
    };
    if (!w.empty() && a_start <= a_end && a_start < w.size()) {
        const uint64_t last = std::min<uint64_t>(a_end, w.size() - 1);
        uint64_t i = a_start;
        for (; i <= last && i % RefBlocks::kB != 0; ++i) node(i);
        for (; i + RefBlocks::kB - 1 <= last; i += RefBlocks::kB) {
            const auto& x = bk[i / RefBlocks::kB];
            lo = std::min(lo, x.first); hi = std::max(hi, x.second);
        }
        for (; i <= last; ++i) node(i);
    }
    if (lo == INT64_MAX) return false;
    start_bp = lo; end_bp = hi;
    return true;
}

// First reference node of `hap`'s walk past the match [i_start, i_end] when the
// match's edge node is off-reference: the reference position where the
// off-reference stretch at that edge (the bubble the edge allele sits in, as
// vg deconstruct positions an alternate allele) rejoins the reference. The
// walk may carry further off-reference nodes first (e.g. its own mutated copy
// of the junction). `after` = the edge after i_end. Returns false when the
// edge node is on-reference or the walk never rejoins the reference.
inline bool anchorPastEdge(const PathView& pv, size_t_ hap,
                           uint32_t i_start, uint32_t i_end, bool after,
                           const RefCoords& refco, const NodeLens& nlen,
                           int64_t& s, int64_t& e)
{
    const auto& w = (*pv.full)[hap];
    if (i_end >= w.size()) return false;
    const uint32_t edge = after ? i_end : i_start;
    if (refco.has(gbwt::Node::id(w[edge]))) return false;
    if (after) {
        for (size_t i = size_t(i_end) + 1; i < w.size(); ++i)
            if (refco.has(gbwt::Node::id(w[i]))) {
                s = refco.bpStart(gbwt::Node::id(w[i]));
                e = s + static_cast<int64_t>(nlen(w[i]));
                return true;
            }
    } else {
        for (size_t i = i_start; i-- > 0; )
            if (refco.has(gbwt::Node::id(w[i]))) {
                s = refco.bpStart(gbwt::Node::id(w[i]));
                e = s + static_cast<int64_t>(nlen(w[i]));
                return true;
            }
    }
    return false;
}

// Reference span of a match. A match that ends (starts) in off-reference
// sequence covers the reference interval that sequence replaces, up to the
// walk's next (previous) reference node: both haplotypes carry the same
// allele there (a shared deletion bridge or substitution), so the span reaches
// that node, even when the walks differ in off-reference sequence before it.
// Otherwise a shared deletion would count as a gap between the pieces
// flanking it. Either haplotype's walk can show the anchor (the other
// may continue into its own mutated copy of the anchor node); when both do,
// the nearer one is taken. Only forward matches read b's walk (an inverted
// match stores b in reverse-walk coordinates).
inline bool refBpRange(const PathView& pv, const SeedOut& sd,
                       const RefCoords& refco, const NodeLens& nlen,
                       int64_t& start_bp, int64_t& end_bp)
{
    int64_t lo, hi;
    if (!refBpRange(pv, sd.hap_a, sd.a_start, sd.a_end, refco, nlen, lo, hi)) return false;
    const uint32_t b_end = sd.b_start + sd.len_nodes - 1;
    int64_t new_hi = INT64_MAX, new_lo = INT64_MIN;
    auto offer = [&](size_t_ hap, uint32_t i0, uint32_t i1, bool after) {
        int64_t s, e;
        if (!anchorPastEdge(pv, hap, i0, i1, after, refco, nlen, s, e)) return;
        if (s >= hi)      new_hi = std::min(new_hi, s);   // walk runs with the reference
        else if (e <= lo) new_lo = std::max(new_lo, e);   // walk runs against it
    };
    for (bool after : {true, false}) {
        offer(sd.hap_a, sd.a_start, sd.a_end, after);
        if (!sd.inverted) offer(sd.hap_b, sd.b_start, b_end, after);
    }
    if (new_hi != INT64_MAX) hi = new_hi;
    if (new_lo != INT64_MIN) lo = new_lo;
    start_bp = lo; end_bp = hi;
    return true;
}

// hap-ibd's min-markers floor on the graph. A marker is a variant site: a
// reference node that some but not all haplotypes walk (after rare masking, so
// a site whose minor allele is below min-mac is not a marker, as hap-ibd drops
// it). A run of such nodes abutting on the reference with the same carrier
// count is one site (a long allele chopped into several nodes). A seed needs
// seed_markers sites, an extension floor(min-extend / min-seed * min-markers)
// (hap-ibd PbwtIbd: minExtendMarkersM1).
struct MarkerSites {
    std::vector<int64_t> bp;   // sorted reference starts of the sites
    uint32_t seed_markers = 0, extend_markers = 0;
    uint32_t count(int64_t sb, int64_t eb) const {
        return static_cast<uint32_t>(std::lower_bound(bp.begin(), bp.end(), eb)
                                     - std::lower_bound(bp.begin(), bp.end(), sb));
    }
};

// The segment builder's piece test, applied at emission so that only gMEMs
// that can be part of a segment are kept: a forward gMEM needs an on-reference
// span of at least min-extend (cM and markers), an inverted one is kept at any
// length.
struct CmFilter {
    const PathView*   pv;
    const RefCoords*  refco;
    const NodeLens*   nlen;
    const GeneticMap* gmap;
    double            min_cM;    // min-extend (effective)
    double            seed_cM;   // min-seed
    const MarkerSites* sites = nullptr;   // marker floor (null = off)
    // Inversion loci (sorted, disjoint); a forward match shorter than min-extend
    // that ends within tol bp of a locus edge is kept as a flank stub: across an
    // inversion one haplotype carries, two IBD haplotypes are identical up to
    // orientation, so the stub can be part of one orientation-aware IBS run.
    const std::vector<std::pair<int64_t, int64_t>>* loci = nullptr;
    int64_t tol = 0;
    bool atLocusEdge(int64_t sb, int64_t eb) const {
        if (!loci || loci->empty()) return false;
        // Loci are sorted and disjoint: an edge within tol of x belongs to the
        // first locus starting at or after x - tol, a later one starting by
        // x + tol, or the locus before them (its end).
        auto near = [&](int64_t x) {
            auto it = std::lower_bound(loci->begin(), loci->end(), std::make_pair(x - tol, INT64_MIN));
            for (auto j = (it == loci->begin() ? it : it - 1); j != loci->end() && j->first <= x + tol; ++j)
                if (std::llabs(j->first - x) <= tol || std::llabs(j->second - x) <= tol) return true;
            return false;
        };
        return near(eb) || near(sb);
    }
    // Fills the gMEM's reference span and cM length when it is kept.
    bool operator()(const SeedOut& s, int64_t& sb, int64_t& eb, double& cM) const {
        if (!s.have_pos) return false;
        if (!refBpRange(*pv, s, *refco, *nlen, sb, eb))
            return false;
        cM = gmap->cMAtBp(eb) - gmap->cMAtBp(sb);
        // Inverted matches are kept at any length: in the piece loop they can
        // fill a hole of the forward matching (see main.cpp).
        if (s.inverted) return true;
        if (cM >= min_cM && (!sites || sites->count(sb, eb) >= sites->extend_markers)) return true;
        return atLocusEdge(sb, eb);   // flank stub (ge_extend false)
    }
    bool isSeed(int64_t sb, int64_t eb, double cM) const {
        return cM >= seed_cM && (!sites || sites->count(sb, eb) >= sites->seed_markers);
    }
    bool isExtend(int64_t sb, int64_t eb, double cM) const {
        return cM >= min_cM && (!sites || sites->count(sb, eb) >= sites->extend_markers);
    }
};

inline void PieceSink::operator()(const SeedOut& s) const {
    int64_t sb, eb; double cM;
    if (!(*keep)(s, sb, eb, cM)) return;
    const uint32_t b_end = s.b_start + s.len_nodes - 1;
    (*buckets)[pieceBucket(s.hap_a, s.hap_b)].push_back({s.hap_a, s.hap_b, static_cast<int32_t>(sb), static_cast<int32_t>(eb),
                       s.a_start, s.a_end, s.b_start, b_end, s.inverted, s.type,
                       keep->isSeed(sb, eb, cM), keep->isExtend(sb, eb, cM)});
}

// Enumerator entry point (enumerate.cpp). It is templated on its sinks so the
// hot emit path can inline; sinks are mutable references because
// SeedRowWriter owns a buffer and a row count.
template <class Emit, class EmitBarrier>
void enumerate_gmems(const gbwt::GBWT& index, const gbwt::FastLocate& fl,
                      const MoveLF& mv, const MovePhi& mp, const NodeLens& nlen, const RefCoords& refco, bool do_typing,
                      uint64_t min_del_bp, uint64_t min_ins_bp, uint64_t floor_bp,
                      uint64_t inv_floor_bp, bool no_inversion, const std::vector<uint32_t>& sample_by_path,
                      const std::vector<uint32_t>& path_len,
                      Emit& emit, EmitBarrier& emit_barrier);

extern template void enumerate_gmems<SeedRowWriter, SeedBufferSink>(
    const gbwt::GBWT&, const gbwt::FastLocate&, const MoveLF&, const MovePhi&, const NodeLens&,
    const RefCoords&, bool, uint64_t, uint64_t, uint64_t, uint64_t, bool,
    const std::vector<uint32_t>&, const std::vector<uint32_t>&,
    SeedRowWriter&, SeedBufferSink&);
extern template void enumerate_gmems<PieceSink, SeedBufferSink>(
    const gbwt::GBWT&, const gbwt::FastLocate&, const MoveLF&, const MovePhi&, const NodeLens&,
    const RefCoords&, bool, uint64_t, uint64_t, uint64_t, uint64_t, bool,
    const std::vector<uint32_t>&, const std::vector<uint32_t>&,
    PieceSink&, SeedBufferSink&);
