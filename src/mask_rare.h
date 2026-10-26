/*
 * Rare-difference masking: the graph form of hap-ibd's min-mac.
 *
 * hap-ibd drops markers whose minor allele is carried by fewer than min-mac
 * haplotypes before it matches anything, so a private difference (assembly or
 * sequencing error, private mutation) never breaks a match. Here a node or an
 * edge walked by fewer than min-mac haplotypes is rare. In each path, a rare
 * stretch is the walk between two consecutive common nodes u and v that holds
 * a rare node or a rare edge; it is rewritten to the route from u to v that the
 * most haplotypes walk (substitution, insertion and deletion alike). Matching
 * runs on the rewritten paths; SV typing reads the original walks through
 * MaskMap. A stretch with no common route (no other haplotype walks from u to v
 * nearby, e.g. the breakpoints of a private inversion) is left as it is.
 */
#pragma once
#include "sweep_common.h"
#include <functional>

// One rewritten stretch: masked walk [m_lo, m_lo + m_len) stands for original
// walk [o_lo, o_lo + o_len), between the same two common nodes.
struct MaskRewrite { uint32_t m_lo, m_len, o_lo, o_len; };

// Masked walk index -> original walk index, per path.
struct MaskMap {
    std::vector<std::vector<MaskRewrite>> rw;   // per path, ascending m_lo

    // First original position that masked position i stands for.
    uint32_t lo(size_t_ h, uint32_t i) const { return map(h, i, false); }
    // Last original position that masked position i stands for.
    uint32_t hi(size_t_ h, uint32_t i) const { return map(h, i, true); }

private:
    uint32_t map(size_t_ h, uint32_t i, bool last) const {
        if (h >= rw.size() || rw[h].empty()) return i;
        const auto& v = rw[h];
        auto it = std::upper_bound(v.begin(), v.end(), i,
                                   [](uint32_t x, const MaskRewrite& r) { return x < r.m_lo; });
        if (it == v.begin()) return i;
        const MaskRewrite& r = *(it - 1);
        if (i < r.m_lo + r.m_len)   // inside a rewritten route
            return last ? r.o_lo + r.o_len - 1 : r.o_lo;
        return i - (r.m_lo + r.m_len) + (r.o_lo + r.o_len);
    }
};

struct MaskStats { size_t rare_nodes = 0, rare_edges = 0, stretches = 0, rewritten = 0; };

// Rewrites `paths` into `masked` (same path order) and fills `map`.
inline MaskStats maskRare(const gbwt::GBWT& index,
                          const std::vector<gbwt::vector_type>& paths,
                          uint32_t min_mac, size_t_ id_lo, size_t_ id_hi,
                          std::vector<gbwt::vector_type>& masked, MaskMap& map)
{
    MaskStats st;
    // Haplotype counts, capped at min_mac. A path counts once per node / edge;
    // an edge and its reverse are the same edge.
    const size_t_ n_ids = id_hi - id_lo + 1;
    std::vector<uint32_t> node_mac(n_ids, 0), node_last(n_ids, UINT32_MAX);
    auto pack = [](node_t x, node_t y) { return (static_cast<uint64_t>(x) << 32) | static_cast<uint64_t>(y); };
    auto edgeKey = [&](node_t x, node_t y) {
        return std::min(pack(x, y), pack(gbwt::Node::reverse(y), gbwt::Node::reverse(x)));
    };
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> edge_mac;   // key -> (last path, count)
    for (size_t_ p = 0; p < paths.size(); ++p) {
        const auto& W = paths[p];
        for (size_t k = 0; k < W.size(); ++k) {
            if (W[k] > UINT32_MAX) throw std::runtime_error("mask_rare: node handle exceeds 32 bits");
            const size_t_ id = gbwt::Node::id(W[k]) - id_lo;
            if (node_last[id] != p) { node_last[id] = static_cast<uint32_t>(p); if (node_mac[id] < min_mac) ++node_mac[id]; }
            if (k == 0) continue;
            auto& e = edge_mac.try_emplace(edgeKey(W[k - 1], W[k]), UINT32_MAX, 0).first->second;
            if (e.first != p) { e.first = static_cast<uint32_t>(p); if (e.second < min_mac) ++e.second; }
        }
    }
    node_last.clear(); node_last.shrink_to_fit();
    for (uint32_t c : node_mac) st.rare_nodes += (c > 0 && c < min_mac);
    for (const auto& kv : edge_mac) st.rare_edges += (kv.second.second < min_mac);
    auto commonNode = [&](node_t w) { return node_mac[gbwt::Node::id(w) - id_lo] >= min_mac; };
    auto rareEdge = [&](node_t x, node_t y) {
        const auto it = edge_mac.find(edgeKey(x, y));
        return it == edge_mac.end() || it->second.second < min_mac;
    };

    // Route from u to v walked by the most haplotypes, at most `cap` nodes
    // between them; empty optional if none reaches min_mac. Depth-first, larger
    // branches first, pruned by the best count so far (counts only shrink).
    struct Cached { size_t cap; bool found; gbwt::vector_type r; };
    std::unordered_map<uint64_t, Cached> cache;
    gbwt::vector_type cur;
    auto route = [&](node_t u, node_t v, size_t cap, gbwt::vector_type& best) -> bool {
        const uint64_t key = pack(u, v);
        const auto ci = cache.find(key);
        if (ci != cache.end() && ci->second.cap == cap) { best = ci->second.r; return ci->second.found; }
        size_t best_n = min_mac - 1;
        bool found = false;
        cur.clear();
        std::function<void(gbwt::SearchState)> dfs = [&](gbwt::SearchState s) {
            std::vector<std::pair<size_t, gbwt::SearchState>> next;
            for (const gbwt::edge_type& e : index.edges(s.node)) {
                if (e.first == gbwt::ENDMARKER) continue;
                const gbwt::SearchState ns = index.extend(s, e.first);
                if (ns.size() > best_n) next.push_back({ns.size(), ns});
            }
            std::sort(next.begin(), next.end(),
                      [](const auto& a, const auto& b) { return a.first > b.first; });
            for (const auto& [n, ns] : next) {
                if (n <= best_n) break;
                if (ns.node == v) { best_n = n; best = cur; found = true; continue; }
                if (cur.size() >= cap) continue;
                cur.push_back(ns.node);
                dfs(ns);
                cur.pop_back();
            }
        };
        dfs(index.find(u));
        cache[key] = {cap, found, found ? best : gbwt::vector_type{}};
        return found;
    };

    masked.assign(paths.size(), {});
    map.rw.assign(paths.size(), {});
    gbwt::vector_type alt;
    for (size_t_ p = 0; p < paths.size(); ++p) {
        const auto& W = paths[p];
        auto& M = masked[p];
        M.reserve(W.size());
        size_t prev = SIZE_MAX;   // last common node's position in W
        for (size_t k = 0; k < W.size(); ++k) {
            if (!commonNode(W[k])) {
                if (prev == SIZE_MAX) M.push_back(W[k]);   // before the first common node: kept
                continue;
            }
            if (prev != SIZE_MAX) {
                const size_t o_len = k - prev - 1;
                if (o_len > 0 || rareEdge(W[prev], W[k])) {
                    ++st.stretches;
                    if (route(W[prev], W[k], o_len + 8, alt)) {
                        ++st.rewritten;
                        map.rw[p].push_back({static_cast<uint32_t>(M.size()), static_cast<uint32_t>(alt.size()),
                                             static_cast<uint32_t>(prev + 1), static_cast<uint32_t>(o_len)});
                        M.insert(M.end(), alt.begin(), alt.end());
                    } else {
                        M.insert(M.end(), W.begin() + prev + 1, W.begin() + k);
                    }
                }
            }
            M.push_back(W[k]);
            prev = k;
        }
        if (prev != SIZE_MAX)   // after the last common node: kept
            M.insert(M.end(), W.begin() + prev + 1, W.end());
    }
    return st;
}

// Bidirectional GBWT of `paths` (forward and reverse of each, in path order),
// carrying `meta` so sample names and path ids are those of the source index.
inline gbwt::GBWT buildBidirGBWT(const std::vector<gbwt::vector_type>& paths,
                                 const gbwt::Metadata& meta, bool has_meta)
{
    gbwt::vector_type text;
    size_t total = 0;
    for (const auto& w : paths) total += 2 * (w.size() + 1);
    text.reserve(total);
    gbwt::vector_type rev;
    for (const auto& w : paths) {
        text.insert(text.end(), w.begin(), w.end());
        text.push_back(gbwt::ENDMARKER);
        rev.clear();
        gbwt::reversePath(w, rev);
        text.insert(text.end(), rev.begin(), rev.end());
        text.push_back(gbwt::ENDMARKER);
    }
    gbwt::DynamicGBWT dynamic;
    dynamic.insert(text, /*has_both_orientations=*/true);
    if (has_meta) {
        dynamic.metadata = meta;
        dynamic.header.set(gbwt::GBWTHeader::FLAG_METADATA);
    }
    return gbwt::GBWT(dynamic);
}
