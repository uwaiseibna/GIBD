/*
 * Move structure for LF over the GBWT (Nishimoto & Tabei, ICALP 2021).
 *
 * One row per stored run of every node record. A row keeps the run's local
 * start and length, its successor node, the LF image of its first row
 * (`dest`, a local offset in the successor's record) and the index of the run
 * containing `dest` (`dest_run`). Because rows of one run share a successor,
 * LF(i) = dest + (i - start) for every row i of the run, and the run holding
 * LF(i) is found by stepping forward from dest_run. A search state that knows
 * the runs holding its first and last rows therefore extends with no record
 * lookup, no header decode and no scan from the record start.
 *
 * Range samples. Every B-th run (global index t = i * B) stores, for each
 * outgoing rank of its record, the first run of that rank at or after t and
 * the last run at or before t, within the record. The distinct successors of a
 * run range of any length then cost O(B + outdegree): scan at most B runs from
 * each end and take the remaining successors from the two samples.
 *
 * Space: 20 bytes per run plus one run offset per record, plus
 * 8 * outdegree bytes per sample.
 */
#pragma once

#include <gbwt/gbwt.h>
#include <gbwt/internal.h>

#include "move_balance.h"

#ifdef __GLIBC__
#include <malloc.h>
#endif

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

struct MoveLF {
    static constexpr uint32_t B    = 32;             // range-sample spacing in runs
    static constexpr uint32_t ALPHA = 8;             // balancing: < 2 * ALPHA steps per seek
    static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max();
    struct Row {
        uint32_t start, len;   // local row interval [start, start + len) in its record
        uint32_t dest;         // LF(start): local offset in record(succ)
        uint32_t dest_run;     // global index of the run containing dest
        uint32_t succ;         // successor node (gbwt::ENDMARKER for path ends)
    };
    std::vector<Row>      rows;        // records in comp order, runs in row order
    std::vector<uint32_t> first_run;   // comp -> index of its first row; size effective + 1
    std::vector<uint32_t> samp_off;    // sample i -> its slots [samp_off[i], samp_off[i + 1])
    std::vector<uint32_t> samp_nxt;    // slot k: first run of rank k at or after the sample, or NONE
    std::vector<uint32_t> samp_prv;    // slot k: last run of rank k at or before the sample, or NONE
    std::vector<bool>     cont;        // row continues the previous row's run (balancing split)
    const gbwt::GBWT*     index = nullptr;

    uint32_t firstRun(gbwt::node_type node) const { return first_run[index->toComp(node)]; }

    // Step from run `k` (holding some row <= pos of the same record) to the run holding pos.
    uint32_t seek(uint32_t k, uint64_t pos) const {
        while ((uint64_t)rows[k].start + rows[k].len <= pos) ++k;
        return k;
    }

    // alpha >= 2 balances the structure (MoveBalancerT, move_balance.h): every output interval
    // then holds fewer than 2 alpha run starts, so every seek takes fewer than
    // 2 alpha steps. Split pieces keep their run's successor and rank; `cont`
    // marks a piece that continues the previous row's run.
    void build(const gbwt::GBWT& gbwt, uint32_t alpha = 0) {
        index = &gbwt;
        const gbwt::size_type eff = gbwt.effective();
        if (gbwt.runs().first >= std::numeric_limits<uint32_t>::max())
            throw std::runtime_error("MoveLF: 2^32 or more runs");
        rows.clear();
        rows.reserve(gbwt.runs().first);
        first_run.assign(eff + 1, 0);
        std::vector<uint32_t> rank;          // outgoing rank of each row
        rank.reserve(gbwt.runs().first);
        std::vector<gbwt::size_type> seen;   // rows of each rank before the current run
        gbwt.bwt.forEach([&](gbwt::size_type c, const gbwt::CompressedRecord& rec) {
            first_run[c] = (uint32_t)rows.size();
            if (rec.outdegree() == 0) return;   // empty record: no runs
            seen.assign(rec.outdegree(), 0);
            gbwt::CompressedRecordIterator it(rec);
            while (!it.end()) {
                const gbwt::rank_type k = it->first;
                const gbwt::size_type len = it->second, end = it.offset();
                const gbwt::node_type to = rec.successor(k);
                Row row{(uint32_t)(end - len), (uint32_t)len, 0, 0, (uint32_t)to};
                if (to != gbwt::ENDMARKER) row.dest = (uint32_t)(rec.offset(k) + seen[k]);
                seen[k] += len;
                rows.push_back(row);
                rank.push_back((uint32_t)k);
                ++it;
            }
        });
        first_run[eff] = (uint32_t)rows.size();
        cont.assign(rows.size(), false);
        if (alpha >= 2) balance(gbwt, alpha, rank);
        // dest_run. LF is order-preserving: over the records in comp order, the
        // rows entering one successor record arrive at increasing offsets. One
        // forward cursor per successor record therefore finds every dest_run
        // in a single pass, O(r) in total.
        std::vector<uint32_t> cur(first_run.begin(), first_run.end() - 1);
        for (Row& row : rows) {
            if (row.succ == gbwt::ENDMARKER) continue;
            uint32_t& k = cur[gbwt.toComp(row.succ)];
            while ((uint64_t)rows[k].start + rows[k].len <= row.dest) ++k;
            row.dest_run = k;
        }
        // Range samples, per record: O(rows + samples * outdegree).
        samp_off.clear(); samp_nxt.clear(); samp_prv.clear();
        std::vector<uint32_t> slot;
        for (size_t c = 0; c < eff; ++c) {
            const size_t f = first_run[c], l = first_run[c + 1];
            if (f == l) continue;
            size_t deg = 0;
            for (size_t j = f; j < l; ++j) deg = std::max<size_t>(deg, rank[j] + 1);
            const size_t i0 = (f + B - 1) / B, i1 = (l + B - 1) / B;
            const size_t base = samp_nxt.size();
            if (base + (i1 - i0) * deg >= NONE) throw std::runtime_error("MoveLF: too many range-sample slots");
            for (size_t i = i0; i < i1; ++i) samp_off.push_back((uint32_t)(base + (i - i0) * deg));
            samp_nxt.resize(base + (i1 - i0) * deg); samp_prv.resize(samp_nxt.size());
            slot.assign(deg, NONE);
            for (size_t j = l; j-- > f; ) {            // first run of each rank at or after j
                slot[rank[j]] = (uint32_t)j;
                if (j % B == 0) std::copy(slot.begin(), slot.end(), samp_nxt.begin() + samp_off[j / B]);
            }
            slot.assign(deg, NONE);
            for (size_t j = f; j < l; ++j) {           // last run of each rank at or before j
                slot[rank[j]] = (uint32_t)j;
                if (j % B == 0) std::copy(slot.begin(), slot.end(), samp_prv.begin() + samp_off[j / B]);
            }
        }
        samp_off.push_back((uint32_t)samp_nxt.size());
        if (samp_off.size() != (rows.size() + B - 1) / B + 1)
            throw std::runtime_error("MoveLF: range samples out of step with the runs");
    }

private:
    // LF as a permutation of all rows: records concatenated in comp order. Rows
    // with a real successor map to their LF image; the endmarker-destined rows,
    // whose images are not stored, are mapped to the endmarker record in row
    // order, which completes the bijection without changing any real image.
    // Peak extra memory is the balancer's two lists, about 32 bytes per row
    // when positions fit in 32 bits: the rows themselves are released while
    // they are balanced, keeping only each row's successor.
    void balance(const gbwt::GBWT& gbwt, uint32_t alpha, std::vector<uint32_t>& rank) {
        const size_t eff = first_run.size() - 1;
        uint64_t N = 0;
        for (size_t c = 0; c < eff; ++c) {
            const size_t l = first_run[c + 1];
            N += (l > first_run[c]) ? (uint64_t)rows[l - 1].start + rows[l - 1].len : 0;
        }
        if (N <= std::numeric_limits<uint32_t>::max()) balanceAs<uint32_t>(gbwt, alpha, rank, N);
        else                                            balanceAs<uint64_t>(gbwt, alpha, rank, N);
    }
    template <class Pos>
    void balanceAs(const gbwt::GBWT& gbwt, uint32_t alpha, std::vector<uint32_t>& rank, uint64_t N) {
        const size_t eff = first_run.size() - 1, R = rows.size();
        std::vector<Pos> base(eff + 1, 0);                 // global start of each record
        for (size_t c = 0; c < eff; ++c) {
            const size_t l = first_run[c + 1];
            base[c + 1] = base[c] + (l > first_run[c] ? (Pos)(rows[l - 1].start + rows[l - 1].len) : 0);
        }
        auto target = [&](size_t j) -> size_t {
            return rows[j].succ == gbwt::ENDMARKER ? 0 : (size_t)gbwt.toComp(rows[j].succ);
        };
        MoveBalancerT<Pos> b;
        b.reserve(R, R / 8);
        for (size_t c = 0; c < eff; ++c)
            for (size_t j = first_run[c]; j < first_run[c + 1]; ++j) b.addInput(base[c] + rows[j].start);
        {
            // Outputs in output order: grouped by target record (counting sort),
            // within a target in row order, which is output order since LF is
            // order-preserving. The endmarker record takes the endmarker rows.
            std::vector<uint32_t> qorder(R), cnt(eff + 1, 0);
            for (size_t j = 0; j < R; ++j) cnt[target(j) + 1]++;
            for (size_t c = 0; c < eff; ++c) cnt[c + 1] += cnt[c];
            for (size_t j = 0; j < R; ++j) qorder[cnt[target(j)]++] = (uint32_t)j;
            Pos ends = 0;
            for (uint32_t j : qorder) {
                if (rows[j].succ == gbwt::ENDMARKER) { b.addOutput(base[0] + ends, j); ends += rows[j].len; }
                else b.addOutput(base[gbwt.toComp(rows[j].succ)] + rows[j].dest, j);
            }
            if (ends != base[1] - base[0]) throw std::runtime_error("MoveLF: endmarker rows do not match the endmarker record");
        }
        std::vector<uint32_t> succ(R), comp(R);            // what the rebuild needs of each row
        for (size_t c = 0; c < eff; ++c)
            for (size_t j = first_run[c]; j < first_run[c + 1]; ++j) { succ[j] = rows[j].succ; comp[j] = (uint32_t)c; }
        std::vector<Row>().swap(rows);
        b.balance((Pos)N, alpha);
        // Pieces, compactly; then the lists go.
        const size_t R2 = b.size();
        if (R2 >= std::numeric_limits<uint32_t>::max()) throw std::runtime_error("MoveLF: 2^32 or more rows");
        std::vector<Pos> pin(R2), pout(R2);
        std::vector<uint32_t> porig(R2);
        { size_t i = 0; b.forEach([&](Pos in, Pos o, uint32_t j) { pin[i] = in; pout[i] = o; porig[i] = j; ++i; }); }
        b.clear();
        // Rebuild the rows from the pieces.
        rows.resize(R2);
        std::vector<uint32_t> rank2(R2);
        cont.assign(R2, false);
        size_t c = 0;
        first_run[0] = 0;
        for (size_t i = 0; i < R2; ++i) {
            const uint32_t j = porig[i];
            const size_t cj = comp[j];
            while (c < cj) first_run[++c] = (uint32_t)i;
            const Pos st = pin[i] - base[cj];
            const Pos next = (i + 1 < R2) ? pin[i + 1] : (Pos)N;
            Row row{(uint32_t)st, (uint32_t)(next - pin[i]), 0, 0, succ[j]};
            if (succ[j] != gbwt::ENDMARKER) row.dest = (uint32_t)(pout[i] - base[gbwt.toComp(succ[j])]);
            cont[i] = (i > 0 && porig[i - 1] == j);
            rows[i] = row; rank2[i] = rank[j];
        }
        while (c < eff) first_run[++c] = (uint32_t)R2;
        rank.swap(rank2);
#ifdef __GLIBC__
        malloc_trim(0);   // hand the balancer's freed lists back to the OS
#endif
    }
public:
};
