/*
 * gMEM enumeration over the bidirectional GBWT.
 *
 * Each haplotype is a string of oriented nodes, so the enumeration needs no
 * order on the graph and handles cycles and inverted repeats directly.
 *
 * enumerate_gmems() is intentionally self-contained: the traversal owns one
 * DFS stack and a set of reusable scratch buffers, with lambdas declared in the
 * order the algorithm uses them:
 *   leftExt             enumerate left extensions in one range-distinct pass
 *   build_root_runs     index the fixed forward record for one DFS root
 *   rightMaxByRuns      test right-maximality from that root index
 *   cellKeys            locate the rows belonging to one context cell
 *   build_cells         overlay left and right contexts
 *   report              pair compatible cells, coalesce, and emit
 *   type_advance        update insertion/deletion typing along the DFS path, O(1) per node
 *
 * Output order is part of the contract. Cell block order, stable sorting, and
 * the last-write-wins drain together decide which witness survives for an equal
 * canonical pair key.
 */
#include "common.h"

// Enumerate the pairwise-maximal matches (gMEMs) of all haplotype pairs.
// gMEMs go to `emit`; inverted gMEMs suppressed by --no-inversion go to
// `emit_barrier`, so that segments cannot be extended across them.
template <class Emit, class EmitBarrier>
void enumerate_gmems(
        const gbwt::GBWT& index,
        const gbwt::FastLocate& fl,
        const MoveLF& mv,
        const MovePhi& mp,
        const NodeLens& nlen,
        const RefCoords& refco,
        bool do_typing,
        uint64_t min_del_bp, uint64_t min_ins_bp,
        uint64_t floor_bp,
        uint64_t inv_floor_bp,   // <= floor_bp: inverted matches are reported from here
        bool no_inversion,
        const std::vector<uint32_t>& sample_by_path,
        const std::vector<uint32_t>& path_len,
        Emit& emit,
        EmitBarrier& emit_barrier)
{
    using gbwt::BidirectionalState;
    // Dense id-indexed node length and reference start (hot path).
    auto NL = [&](node_t w) -> uint64_t { return nlen(w); };
    auto RS = [&](node_t w) -> int64_t  { return refco.raw(gbwt::Node::id(w)); };
    // The forward record is fixed for an entire DFS root subtree because left
    // extension does not change forward.node. The backward record is different
    // at each state, so it stays transient rather than growing a graph-wide
    // cache.
    const gbwt::CompressedRecord* fwd_rec = nullptr;

    // Distinct non-empty backward extensions of a bidirectional interval, read
    // from the move structure. `brun` and `erun` are the move rows holding the
    // first and last rows of st.backward.range. A range of at most 2B runs is
    // read in one pass; a longer one scans at most B runs from each end and
    // takes every other successor's first and last run from the two range
    // samples, so a state costs O(B + outdegree) however many runs it spans.
    // A successor's rows map to consecutive rows of its record, so its count is
    // the distance between the LF images of its first and last row in range.
    // Each child's backward range is that image; its forward range is the
    // parent's, shifted by the rows whose successor x has reverse(x) <
    // reverse(to) -- the bidirectional synchronisation of CompressedRecord::bdLF.
    struct Child { node_t first; BidirectionalState second; uint32_t run, erun; };
    using Ext = std::vector<Child>;
    struct Sym { uint32_t succ; uint32_t cnt; uint32_t first_pos; uint32_t first_run; uint32_t last_run; };
    std::vector<Sym> syms;             // successors present in the range, reused
    auto symOf = [&](uint32_t succ) -> Sym* {
        for (Sym& x : syms) if (x.succ == succ) return &x;
        return nullptr;
    };
    // One pass over runs [brun, erun]: first and last run and row count of each successor.
    auto scanRange = [&](uint64_t s, uint64_t e, uint32_t brun, uint32_t erun) {
        for (uint32_t j = brun; j <= erun; ++j) {
            const MoveLF::Row& R = mv.rows[j];
            const uint64_t lo = std::max<uint64_t>(s, R.start);
            const uint64_t hi = std::min<uint64_t>(e, (uint64_t)R.start + R.len - 1);
            Sym* y = symOf(R.succ);
            if (!y) { syms.push_back({R.succ, 0, (uint32_t)lo, j, j}); y = &syms.back(); }
            y->cnt += (uint32_t)(hi - lo + 1);
            y->last_run = j;
        }
    };
    // Long range: the same result from O(B) runs at each end plus two samples.
    auto sampleRange = [&](uint64_t s, uint64_t e, uint32_t brun, uint32_t erun) {
        constexpr uint32_t B = MoveLF::B, NONE = MoveLF::NONE;
        const uint32_t t1 = (brun + B - 1) / B * B, t2 = erun / B * B;   // brun <= t1 < t2 <= erun
        for (uint32_t j = brun; j < t1; ++j) {          // first runs before the left sample
            const MoveLF::Row& R = mv.rows[j];
            if (!symOf(R.succ)) syms.push_back({R.succ, 0, (uint32_t)std::max<uint64_t>(s, R.start), j, NONE});
        }
        for (uint32_t p = mv.samp_off[t1 / B]; p < mv.samp_off[t1 / B + 1]; ++p) {
            const uint32_t a = mv.samp_nxt[p];
            if (a == NONE || a > erun) continue;
            const MoveLF::Row& R = mv.rows[a];
            if (!symOf(R.succ)) syms.push_back({R.succ, 0, (uint32_t)std::max<uint64_t>(s, R.start), a, NONE});
        }
        for (uint32_t j = erun; j > t2; --j) {          // last runs after the right sample
            Sym* y = symOf(mv.rows[j].succ);
            if (y->last_run == NONE) y->last_run = j;
        }
        for (uint32_t p = mv.samp_off[t2 / B]; p < mv.samp_off[t2 / B + 1]; ++p) {
            const uint32_t b = mv.samp_prv[p];
            if (b == NONE || b < brun) continue;
            Sym* y = symOf(mv.rows[b].succ);
            if (y && y->last_run == NONE) y->last_run = b;
        }
        uint64_t real = 0;                               // rows with a real successor
        for (Sym& y : syms) {
            if (y.succ == gbwt::ENDMARKER) continue;
            const MoveLF::Row& A = mv.rows[y.first_run];
            const MoveLF::Row& Z = mv.rows[y.last_run];
            const uint64_t cs = (uint64_t)A.dest + (y.first_pos - A.start);
            const uint64_t ce = (uint64_t)Z.dest + (std::min<uint64_t>(e, (uint64_t)Z.start + Z.len - 1) - Z.start);
            y.cnt = (uint32_t)(ce - cs + 1);
            real += y.cnt;
        }
        if (Sym* z = symOf(gbwt::ENDMARKER)) z->cnt = (uint32_t)((e - s + 1) - real);
    };
    auto leftExt = [&](const BidirectionalState& st, uint32_t brun, uint32_t erun, Ext& out) {
        out.clear(); syms.clear();
        const uint64_t s = st.backward.range.first, e = st.backward.range.second;
        if (erun - brun + 1 <= 2 * MoveLF::B) scanRange(s, e, brun, erun);
        else {
            sampleRange(s, e, brun, erun);
        }
        if (syms.size() > 1)
            std::sort(syms.begin(), syms.end(),
                      [](const Sym& a, const Sym& b) { return a.succ < b.succ; });
        for (const Sym& y : syms) {
            if (y.succ == gbwt::ENDMARKER) continue;
            const node_t rto = gbwt::Node::reverse(y.succ);
            gbwt::size_type ro = 0;
            for (const Sym& x : syms)
                if (x.succ != y.succ && gbwt::Node::reverse(x.succ) < rto) ro += x.cnt;
            const MoveLF::Row& R = mv.rows[y.first_run];
            const gbwt::size_type cs = (gbwt::size_type)R.dest + (y.first_pos - R.start);
            BidirectionalState ns;
            ns.forward.node   = st.forward.node;
            ns.forward.range  = {st.forward.range.first + ro, st.forward.range.first + ro + y.cnt - 1};
            ns.backward.node  = y.succ;
            ns.backward.range = {cs, cs + y.cnt - 1};
            // Last row's run: from the LF image of the successor's last run in range.
            const MoveLF::Row& Z = mv.rows[y.last_run];
            out.push_back({rto, ns, mv.seek(R.dest_run, cs), mv.seek(Z.dest_run, cs + y.cnt - 1)});
        }
    };
    Ext ext_buf;                       // reused across states (see leftExt above)

    // Per-root run index for the fixed forward record. It stores all runs in
    // order plus a CSR view by outgoing rank, and two lookup tables that answer
    // the right-maximality tests and context-cell queries in O(1): the run
    // holding each row, and for every run j and rank k the first run of rank k
    // at or after j. Building them costs O(len + nsym * runs) per root, which
    // is O(n) over all roots when node degree is bounded.
    struct RootRuns {
        std::vector<uint32_t> starts, syms;      // run start offset / outgoing rank, ascending
        std::vector<uint32_t> sym_off, sym_runs; // CSR: runs of rank k = sym_runs[sym_off[k]..)
        std::vector<uint32_t> logical;           // FastLocate logical run index of each run head
        std::vector<uint32_t> run_of;            // row -> run holding it
        std::vector<uint32_t> nxt;               // (run j, rank k) -> CSR position of the first
                                                 // run of rank k at or after j; row j*nsym
        gbwt::size_type len = 0;                 // record length (one past the last run)
        size_t nsym = 0;
        // Runs partition [0, len); return the run containing offset i.
        size_t runAt(gbwt::size_type i) const { return run_of[i]; }
        gbwt::size_type runEnd(size_t ri) const {          // one past run ri
            return (ri + 1 < starts.size()) ? starts[ri + 1] : len;
        }
        // First run of rank k that can intersect a range starting at lo. Runs
        // before the one holding lo end at or before lo; that run and all later
        // ones end after it.
        size_t firstRunOfSym(uint32_t k, gbwt::size_type lo) const {
            return nxt[runAt(lo) * nsym + k];
        }
        bool symInRange(uint32_t k, gbwt::range_type range) const {
            const size_t p = firstRunOfSym(k, range.first);
            return p < sym_off[k + 1] && starts[sym_runs[p]] <= range.second;
        }
    };
    RootRuns rr;
    std::vector<uint32_t> sym_cursor;
    auto build_root_runs = [&](const gbwt::CompressedRecord& rec) {
        rr.nsym = rec.outdegree();
        rr.starts.clear(); rr.syms.clear();
        gbwt::CompressedRecordIterator it(rec);
        gbwt::size_type end = 0;
        while (!it.end()) {
            end = it.offset();
            rr.starts.push_back((uint32_t)(end - it->second));
            rr.syms.push_back((uint32_t)it->first);
            ++it;
        }
        rr.len = end;
        // Logical runs number FastLocate's samples: an endmarker-destined run
        // is one logical run per row, any other run is one logical run.
        rr.logical.resize(rr.starts.size());
        uint32_t lg = 0;
        for (size_t i = 0; i < rr.starts.size(); ++i) {
            rr.logical[i] = lg;
            lg += (rec.successor(rr.syms[i]) == gbwt::ENDMARKER)
                      ? (uint32_t)(rr.runEnd(i) - rr.starts[i]) : 1u;
        }
        rr.sym_off.assign(rr.nsym + 1, 0);
        for (uint32_t s : rr.syms) rr.sym_off[s + 1]++;
        for (size_t k = 0; k < rr.nsym; ++k) rr.sym_off[k + 1] += rr.sym_off[k];
        rr.sym_runs.resize(rr.syms.size());
        sym_cursor.assign(rr.sym_off.begin(), rr.sym_off.end());   // reused across roots
        for (size_t i = 0; i < rr.syms.size(); ++i)
            rr.sym_runs[sym_cursor[rr.syms[i]]++] = (uint32_t)i;
        // O(1) lookup tables. sym_cursor[k] now equals sym_off[k + 1], so a
        // backward pass hands out each rank's CSR positions from the end.
        const size_t R = rr.starts.size(), K = rr.nsym;
        rr.run_of.resize(rr.len);
        for (size_t i = 0; i < R; ++i)
            std::fill(rr.run_of.begin() + rr.starts[i], rr.run_of.begin() + rr.runEnd(i), (uint32_t)i);
        rr.nxt.resize((R + 1) * K);
        for (size_t k = 0; k < K; ++k) rr.nxt[R * K + k] = rr.sym_off[k + 1];
        for (size_t i = R; i-- > 0; ) {
            std::copy(rr.nxt.begin() + (i + 1) * K, rr.nxt.begin() + (i + 2) * K, rr.nxt.begin() + i * K);
            rr.nxt[i * K + rr.syms[i]] = --sym_cursor[rr.syms[i]];
        }
    };
    // A range of at least two rows is right-maximal if it spans two successor
    // runs, or if its only successor class is the endmarker.
    auto rightMaxByRuns = [&](const gbwt::CompressedRecord& rec, gbwt::range_type range) -> bool {
        if (range.first > range.second || range.second - range.first + 1 < 2) return false;
        const size_t ri = rr.runAt(range.first);
        if (rr.runEnd(ri) <= range.second) return true;      // a second run intersects
        return rec.successor(rr.syms[ri]) == gbwt::ENDMARKER;
    };
    auto right_maximal_state = [&](const gbwt::CompressedRecord& rec,
                                   const BidirectionalState& S) {
        return !S.empty() && rightMaxByRuns(rec, S.forward.range);
    };

    // Per-state candidate pairs. The packed key orders by
    // (hap_a, hap_b, inverted, type); stable sorting then groups equivalent
    // records while preserving the last witness in each group.
    const uint32_t n_paths = (uint32_t)path_len.size();
    unsigned path_bits = 1;
    while (path_bits < 31 && ((uint64_t)1 << path_bits) < n_paths) path_bits++;
    if (n_paths > ((uint64_t)1 << path_bits))
        throw std::runtime_error("2^31 or more paths; packed pair key would overflow");
    const unsigned key_bits = 2 * path_bits + 2;
    const unsigned key_bytes = (key_bits + 7) / 8;
    using LKey = uint64_t;
    auto pack_key = [&](uint32_t pa, uint32_t pb, unsigned inv, unsigned type) -> LKey {
        return ((LKey)pa << (path_bits + 2)) | ((LKey)pb << 2) | ((LKey)inv << 1) | (LKey)type;
    };
    std::vector<std::pair<LKey, std::pair<uint32_t,uint32_t>>> localbuf, radix_buf;
    std::vector<uint32_t> radix_cnt;

    // A context cell C_W(l,r) is represented lazily as a row interval and an
    // outgoing rank. Occurrences are located only if the cell has a compatible
    // partner.
    struct Cell {
        node_t l, r, node;
        gbwt::range_type row;
        uint32_t col;                 // outgoing rank of r in the root record
    };
    std::vector<Cell> cells;
    // Overlay scratch, reused across states.
    std::vector<uint32_t> cols;
    std::vector<gbwt::range_type> rows, left_term;
    std::vector<gbwt::range_type> lr_scratch;    // build_cells left-interval sort scratch
    // report() scratch. Inner vectors keep their capacity across states, which
    // avoids repeated allocation during cell materialization.
    std::vector<char> needed_buf;
    std::vector<std::vector<uint64_t>> seqs_buf;

    // Locate one lazy cell. Every cell of a root's subtree is a row interval
    // of the same root record, so located rows are kept per root: each run of
    // the root record owns an SA cursor that is extended by Phi only as far as
    // the deepest row any cell has asked for. A row is therefore Phi-located
    // at most once per root, never re-walked from its run head, and its run is
    // found through the root run index instead of a record scan. A cursor
    // starts in its head's recorded Move-Phi interval, so no step searches.
    std::vector<uint64_t> sa_row;      // SA value of each root-record row, valid below sa_filled
    std::vector<uint32_t> sa_filled;   // rows of each root run located so far
    std::vector<uint32_t> sa_k;        // Move-Phi interval of each run's last located row
    gbwt::size_type root_run_base = 0; // FastLocate run id of the root record's first run
    auto cellKeys = [&](const Cell& c, std::vector<uint64_t>& out) {
        out.clear();
        for (size_t p = rr.firstRunOfSym(c.col, c.row.first); p < rr.sym_off[c.col + 1]; ++p) {
            const size_t ri = rr.sym_runs[p];
            if (rr.starts[ri] > c.row.second) break;
            const gbwt::size_type rs = rr.starts[ri];
            const gbwt::size_type a = std::max<gbwt::size_type>(rs, c.row.first);
            const gbwt::size_type b = std::min<gbwt::size_type>(rr.runEnd(ri) - 1, c.row.second);
            uint32_t& f = sa_filled[ri];
            const uint32_t need = (uint32_t)(b - rs + 1);
            if (f < need) {
                if (f == 0) {
                    const size_t run = root_run_base + rr.logical[ri];
                    sa_row[rs] = fl.samples[run]; sa_k[ri] = (uint32_t)mp.head[run]; f = 1;
                }
                uint32_t k = sa_k[ri];
                for (; f < need; ++f) sa_row[rs + f] = mp.next(sa_row[rs + f - 1], k);
                sa_k[ri] = k;
            }
            for (gbwt::size_type i = a; i <= b; ++i) {
                const uint64_t q = sa_row[i];
                out.push_back((fl.seqId(q) << 32) | (uint64_t)fl.seqOffset(q));
            }
        }
    };

    // Build the context-cell grid for one state. Cell block order is part of
    // the output contract because equal canonical pairs keep the last witness.
    auto build_cells = [&](const BidirectionalState& S, const Ext& L) {
        const node_t TERMINAL = gbwt::ENDMARKER;
        cells.clear();
        const node_t seed = S.forward.node;
        const gbwt::CompressedRecord& rec = *fwd_rec;

        // Columns are the successor ranks present in the state's interval, kept
        // in record order for deterministic cell construction.
        cols.clear();
        for (uint32_t k = 0; k < rr.nsym; ++k) {
            if (rr.symInRange(k, S.forward.range)) cols.push_back(k);
        }
        if (cols.empty()) return;
        const size_t ncols = cols.size();

        // Rows are ordinary left-extension ranges followed by terminal-left
        // residual ranges. Together they partition S.forward.range.
        rows.clear();
        for (auto& le : L) rows.push_back(le.second.forward.range);
        const size_t nreal = rows.size();
        left_term.clear();
        {
            std::vector<gbwt::range_type>& lr = lr_scratch;   // reused across calls
            lr.assign(rows.begin(), rows.end());
            std::sort(lr.begin(), lr.end());
            gbwt::size_type cur = S.forward.range.first;
            for (auto& rg : lr) {
                if (rg.first > cur) left_term.push_back({cur, rg.first - 1});
                if (rg.second + 1 > cur) cur = rg.second + 1;
            }
            if (cur <= S.forward.range.second) left_term.push_back({cur, S.forward.range.second});
        }
        for (auto& rg : left_term) rows.push_back(rg);

        // A cell exists iff at least one run of `col` intersects `row`.
        auto push_cell = [&](node_t l, node_t r, size_t row, size_t col) {
            const gbwt::range_type rg = rows[row];
            if (rg.first > rg.second) return;
            if (!rr.symInRange(cols[col], rg)) return;
            cells.push_back({l, r, seed, rg, cols[col]});
        };
        size_t ecol = ncols;   // the endmarker column, if present
        for (size_t c = 0; c < ncols; ++c)
            if (rec.successor(cols[c]) == gbwt::ENDMARKER) { ecol = c; break; }
        // 1. ordinary (a,b): right-major, then left in extension order.
        for (size_t c = 0; c < ncols; ++c) {
            if (c == ecol) continue;
            const node_t b = rec.successor(cols[c]);
            for (size_t i = 0; i < nreal; ++i) push_cell(L[i].first, b, i, c);
        }
        if (ecol < ncols) {
            // 2. (a, TERMINAL) in L order.
            for (size_t i = 0; i < nreal; ++i) push_cell(L[i].first, TERMINAL, i, ecol);
            // 3. (TERMINAL, TERMINAL) per residual interval, ascending.
            for (size_t j = 0; j < left_term.size(); ++j)
                push_cell(TERMINAL, TERMINAL, nreal + j, ecol);
        }
        // 4. (TERMINAL, b): right-major, then residual interval ascending.
        for (size_t c = 0; c < ncols; ++c) {
            if (c == ecol) continue;
            const node_t b = rec.successor(cols[c]);
            for (size_t j = 0; j < left_term.size(); ++j) push_cell(TERMINAL, b, nreal + j, c);
        }
    };

    // Pair compatible cells, deduplicate the state-local candidate list, and
    // emit SeedOut records. S is already known to be right-maximal.
    auto report = [&](const BidirectionalState& S, uint32_t dn, uint64_t dbp,
                      bool span_ins, bool span_del) {
        const node_t end_node = S.forward.node;
        if (dbp < inv_floor_bp) return;
        if (S.size() < 2) return;
        const bool fwd_ok = dbp >= floor_bp;   // forward pairs need the full floor
        const node_t TERMINAL = gbwt::ENDMARKER;

        localbuf.clear();
        auto emit_pair = [&](uint64_t ka, uint64_t kb){
            uint64_t seq_a = ka >> 32, seq_b = kb >> 32;
            uint32_t off_a = (uint32_t)ka, off_b = (uint32_t)kb;
            uint64_t pa = gbwt::Path::id(seq_a), pb = gbwt::Path::id(seq_b);
            if (pa == pb) return;
            bool inv = gbwt::Path::is_reverse(seq_a) != gbwt::Path::is_reverse(seq_b);
            if (!inv && !fwd_ok) return;
            if (pa > pb) { std::swap(pa, pb); std::swap(seq_a, seq_b); std::swap(off_a, off_b); }
            // Keep one of the two mirror enumerations: the one in which the
            // smaller-path-id sequence is forward-oriented. A non-inverted pair
            // is enumerated once with both sequences forward and once with both
            // reverse, an inverted pair once with each orientation on either
            // side, so exactly one enumeration qualifies, and its offsets are
            // forward-sequence offsets. (The mirror is enumerated at the reverse
            // of the subpath's first node, so the rule cannot key on the ending
            // node.)
            if (gbwt::Path::is_reverse(seq_a)) return;
            int type = (sample_by_path[pa] == sample_by_path[pb]) ? 1 : 0;
            localbuf.push_back({pack_key((uint32_t)pa, (uint32_t)pb, inv ? 1u : 0u,
                                        (unsigned)type), {off_a, off_b}});
        };
        auto side_ok = [&](node_t a, node_t b){ return a != b || a == TERMINAL; };

        // Materialize only cells that have at least one compatible partner.

        std::vector<char>& needed = needed_buf;
        needed_buf.assign(cells.size(), 0);
        for (size_t a = 0; a < cells.size(); ++a)
            for (size_t b = a; b < cells.size(); ++b) {
                if (!side_ok(cells[a].r, cells[b].r) || !side_ok(cells[a].l, cells[b].l)) continue;
                needed[a] = needed[b] = 1;
            }
        if (seqs_buf.size() < cells.size()) seqs_buf.resize(cells.size());
        std::vector<std::vector<uint64_t>>& seqs = seqs_buf;
        // Clear all visible slots; cells rejected by needed[] must look empty.
        for (size_t i = 0; i < cells.size(); ++i) seqs[i].clear();
        for (size_t i = 0; i < cells.size(); ++i) {
            if (!needed[i]) continue;
            cellKeys(cells[i], seqs[i]);
        }

        for (size_t a = 0; a < cells.size(); ++a)
            for (size_t b = a; b < cells.size(); ++b) {
                if (!side_ok(cells[a].r, cells[b].r) || !side_ok(cells[a].l, cells[b].l)) continue;
                if (a == b) {
                    const auto& g = seqs[a];
                    for (size_t i = 0; i < g.size(); ++i)
                        for (size_t j = i + 1; j < g.size(); ++j) emit_pair(g[i], g[j]);
                } else {
                    for (uint64_t x : seqs[a]) for (uint64_t y : seqs[b]) emit_pair(x, y);
                }
            }

        // Stable sort by packed key. Large buffers use LSD radix passes; small
        // ones use stable_sort with the same ordering.
        if (localbuf.size() < 128) {
            std::stable_sort(localbuf.begin(), localbuf.end(),
                             [](const auto& a, const auto& b){ return a.first < b.first; });
        } else {
            radix_buf.resize(localbuf.size());
            radix_cnt.resize(256);
            for (unsigned pass = 0; pass < key_bytes; ++pass) {
                const unsigned shift = pass * 8;
                std::fill(radix_cnt.begin(), radix_cnt.end(), 0u);
                for (const auto& e : localbuf) radix_cnt[(e.first >> shift) & 0xFF]++;
                uint32_t sum = 0;
                for (unsigned d = 0; d < 256; ++d) { const uint32_t c = radix_cnt[d];
                                                     radix_cnt[d] = sum; sum += c; }
                for (const auto& e : localbuf) radix_buf[radix_cnt[(e.first >> shift) & 0xFF]++] = e;
                localbuf.swap(radix_buf);
            }
        }
        for (size_t bi = 0; bi < localbuf.size(); ) {
            size_t bj = bi;                                     // keep the last equal-key witness
            while (bj + 1 < localbuf.size() &&
                   localbuf[bj + 1].first == localbuf[bi].first) ++bj;
            const auto& t   = localbuf[bj].first;
            const auto& val = localbuf[bj].second;
            uint32_t pa = (uint32_t)(t >> (path_bits + 2));
            uint32_t pb = (uint32_t)((t >> 2) & (((LKey)1 << path_bits) - 1));
            bool inv = (t >> 1) & 1;
            SeedOut so{};
            so.hap_a = pa; so.hap_b = pb; so.inverted = inv;
            so.end_node = end_node; so.len_nodes = dn; so.len_bp = dbp;
            so.type = (t & 1) ? SegType::HBD : SegType::IBD;
            so.a_end   = path_len[pa] - 1 - val.first;
            so.a_start = so.a_end - dn + 1;
            so.b_start = (path_len[pb] - 1 - val.second) - dn + 1;
            so.have_pos = true;
            if (do_typing && !inv) { so.has_insertion = span_ins; so.has_deletion = span_del; }
            if (inv && no_inversion) emit_barrier(so);
            else                     emit(so);
            bi = bj + 1;
        }
        // Do not keep a rare high-multiplicity buffer for the rest of the run.
        constexpr size_t kArenaCap = 1u << 14;   // ~16K entries
        if (localbuf.capacity() > kArenaCap) {
            std::vector<std::pair<LKey, std::pair<uint32_t,uint32_t>>>().swap(localbuf);
            std::vector<std::pair<LKey, std::pair<uint32_t,uint32_t>>>().swap(radix_buf);
        }
    };

    // Incremental SV span typing along the DFS path. The verdict is an OR over
    // consecutive on-reference node pairs, scored in the reference direction
    // of the whole span, and a prepended node adds exactly one pair. Each pair
    // is scored under BOTH directions, so a change of direction only selects
    // the other accumulated verdict: O(1) per node, no path rescan.
    struct TypeState {
        bool ins_asc = false, del_asc = false, ins_desc = false, del_desc = false;
        int64_t first_rs = -1, first_end = -1, last_rs = -1, bridge = 0;
        bool asc = true; uint32_t non = 0;
        bool ins() const { return asc ? ins_asc : ins_desc; }
        bool del() const { return asc ? del_asc : del_desc; }
    };
    // Extend the carried typing state by one prepended node.
    auto type_advance = [&](TypeState& ts, node_t w) {
        const int64_t rs = RS(w), nl = (int64_t)NL(w);
        if (rs < 0) { if (ts.non >= 1) ts.bridge += nl; return; }   // off-ref: accrue between-bridge
        if (ts.non == 0) { ts.first_rs = rs; ts.first_end = rs + nl; ts.last_rs = rs;
                           ts.bridge = 0; ts.asc = true; ts.non = 1; return; }
        // w pairs with the previous leftmost on-reference node.
        const int64_t gap_asc  = ts.first_rs - (rs + nl);   // ascending: w is the new scan start
        const int64_t gap_desc = rs - ts.first_end;          // descending: w is the new scan end
        if (gap_asc - ts.bridge >= (int64_t)min_del_bp)  ts.del_asc  = true;
        if (ts.bridge - gap_asc >= (int64_t)min_ins_bp)  ts.ins_asc  = true;
        if (gap_desc - ts.bridge >= (int64_t)min_del_bp) ts.del_desc = true;
        if (ts.bridge - gap_desc >= (int64_t)min_ins_bp) ts.ins_desc = true;
        ts.asc = !(ts.last_rs < rs);   // span direction: leftmost vs rightmost on-ref node
        ts.first_rs = rs; ts.first_end = rs + nl; ts.bridge = 0; ts.non += 1;
    };

    // DFS frames extend matches to the left while the end node stays fixed.
    // `brun` / `erun` are the move rows holding the first / last backward row.
    // `ptype` is the parent's typing state before this frame's node is added.
    struct Frame { BidirectionalState st; uint32_t brun, erun; uint32_t dn; uint64_t dbp; node_t lnode; TypeState ptype; };
    std::vector<Frame> stk;
    // Roots in descending node order; each root owns one forward record index
    // for its whole subtree.
    for (node_t v = index.sigma(); v-- > index.firstNode(); ) {
        if (!index.contains(v) || index.empty(v)) continue;
        BidirectionalState st = index.bdFind(v);
        gbwt::CompressedRecord root_rec = index.record(v);
        // Decode this root's forward runs before the root right-maximality test.
        build_root_runs(root_rec);
        if (!right_maximal_state(root_rec, st)) continue;
        fwd_rec = &root_rec;
        root_run_base = fl.comp_to_run[index.toComp(v)];
        sa_filled.assign(rr.starts.size(), 0);
        if (sa_k.size() < rr.starts.size()) sa_k.resize(rr.starts.size());
        if (sa_row.size() < rr.len) sa_row.resize(rr.len);
        // The root spans its whole backward record: one pass over its runs,
        // O(r) over all roots.
        const uint32_t b0 = mv.firstRun(st.backward.node);
        stk.push_back({st, b0, mv.seek(b0, st.backward.range.second), 1, NL(v), v, TypeState{}});
        while (!stk.empty()) {
            Frame f = std::move(stk.back()); stk.pop_back();
            TypeState ts;
            // Conduit transitions in O(1). When the backward range lies inside
            // one move row whose successor is a real node, every occurrence
            // has the same predecessor: the state emits nothing and its only
            // child is the row's LF image with the same forward range. That
            // child would be pushed and popped next anyway, so advance in place.
            // The last row's run is needed only where the chain ends; it lies in
            // the image of the chain's final move row, so one seek finds it.
            bool chained = false;
            for (;;) {
                ts = f.ptype;
                if (do_typing) type_advance(ts, f.lnode);
                const MoveLF::Row& R = mv.rows[f.brun];
                const gbwt::size_type bs = f.st.backward.range.first, be = f.st.backward.range.second;
                if (R.succ == gbwt::ENDMARKER || be >= (gbwt::size_type)R.start + R.len) break;
                const gbwt::size_type cs = (gbwt::size_type)R.dest + (bs - R.start);
                const node_t ln = gbwt::Node::reverse(R.succ);
                f.st.backward.node  = R.succ;
                f.st.backward.range = {cs, cs + (be - bs)};
                f.brun  = mv.seek(R.dest_run, cs);
                chained = true;
                f.dn   += 1;
                f.dbp  += NL(ln);
                f.lnode = ln;
                f.ptype = ts;
            }
            if (chained) {
                f.erun = mv.seek(f.brun, f.st.backward.range.second);
            }
            Ext& L = ext_buf;
            leftExt(f.st, f.brun, f.erun, L);
            // A conduit has one ordinary left extension covering the whole
            // interval, so every cell it could build carries the same
            // non-terminal left context and no cell pair passes side_ok. It
            // therefore emits nothing, but its child may still become
            // left-maximal, so the DFS continues through it.
            const bool conduit = (L.size() == 1 &&
                                  L[0].second.forward.range == f.st.forward.range);
            if (f.dbp >= inv_floor_bp && f.st.size() >= 2 && !conduit) {
                build_cells(f.st, L);
                report(f.st, f.dn, f.dbp, ts.ins(), ts.del());
            }
            if (conduit) {
                // The conduit child has the same rows and successor runs, so
                // it inherits right-maximality from the current state.
                stk.push_back({L[0].second, L[0].run, L[0].erun, f.dn + 1, f.dbp + NL(L[0].first), L[0].first, ts});
            } else {
                // Each child is tested on its own sub-range of the same record.
                for (auto& c : L) {
                    if (rightMaxByRuns(root_rec, c.second.forward.range)) {
                        stk.push_back({c.second, c.run, c.erun, f.dn + 1, f.dbp + NL(c.first), c.first, ts});
                    }
                }
            }
        }
        fwd_rec = nullptr;
    }
}

// Explicit instantiations for the sink pairs used by main.cpp.
template void enumerate_gmems<SeedRowWriter, SeedBufferSink>(
    const gbwt::GBWT&, const gbwt::FastLocate&, const MoveLF&, const MovePhi&, const NodeLens&,
    const RefCoords&, bool, uint64_t, uint64_t, uint64_t, uint64_t, bool,
    const std::vector<uint32_t>&, const std::vector<uint32_t>&,
    SeedRowWriter&, SeedBufferSink&);
template void enumerate_gmems<PieceSink, SeedBufferSink>(
    const gbwt::GBWT&, const gbwt::FastLocate&, const MoveLF&, const MovePhi&, const NodeLens&,
    const RefCoords&, bool, uint64_t, uint64_t, uint64_t, uint64_t, bool,
    const std::vector<uint32_t>&, const std::vector<uint32_t>&,
    PieceSink&, SeedBufferSink&);
