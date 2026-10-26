/*
 * Move structure for Phi^{-1} over FastLocate's samples (SA[i] -> SA[i+1]).
 *
 * FastLocate answers locateNext(p) as samples[last_to_run[k] + 1] + (p - t_k),
 * where t_k is the largest run-tail text position <= p, found by a predecessor
 * search in an Elias-Fano vector. The tails cut the text into intervals
 * [t_k, t_{k+1}), and Phi^{-1} shifts each interval as a block. Storing, per
 * interval, its start, its image `dest` and the interval holding `dest` turns a
 * walk of consecutive rows into one row read plus a short forward seek per row
 * (Nishimoto & Tabei; Brown, Gagie & Rossi). A walk starts at a run head,
 * and the interval holding every run head is recorded while FastLocate is
 * sampled (buildFastLocate), so neither the build nor a walk searches.
 */
#pragma once

#include <gbwt/fast_locate.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#ifndef MOVE_FF   // defined by move_lf.h (fast-forward counters)
#define MOVE_FF(site, n) ((void)0)
#endif

struct MovePhi {
    struct __attribute__((packed)) Row {
        uint64_t start;      // run-tail text position t_k (packed seq/offset)
        uint64_t dest;       // Phi^{-1}(t_k) = head sample of the next run
        uint32_t dest_run;   // interval holding dest
    };
    std::vector<Row>    rows;
    sdsl::int_vector<0> head;   // interval holding samples[j], for every logical run j

    // `head_row` is the head-interval array produced by buildFastLocate. The
    // row of tail k moves to the head of run last_to_run[k] + 1, so its
    // dest_run is that run's head interval: O(r), no search.
    // `seq_len` (path lengths from buildFastLocate), when given, caps the
    // intervals; see cap().
    void build(const gbwt::FastLocate& fl, sdsl::int_vector<0>&& head_row,
               const std::vector<uint32_t>* seq_len = nullptr) {
        const size_t R = fl.last_to_run.size();
        if (R >= std::numeric_limits<uint32_t>::max())
            throw std::runtime_error("MovePhi: 2^32 or more runs");
        if (head_row.size() != fl.samples.size())
            throw std::runtime_error("MovePhi: head intervals do not match the samples");
        head = std::move(head_row);
        if (seq_len) rows.reserve(cappedSize(fl, *seq_len));   // capping then expands in place
        rows.resize(R);
        size_t k = 0;
        for (auto it = fl.last.one_begin(); it != fl.last.one_end(); ++it, ++k) {
            const size_t run = fl.last_to_run[k];
            rows[k].start = it->second;
            if (run + 1 < fl.samples.size()) {
                rows[k].dest     = fl.samples[run + 1];
                rows[k].dest_run = (uint32_t)head[run + 1];
            } else {
                rows[k].dest = 0; rows[k].dest_run = 0;   // the interval holding 0
            }
        }
        if (seq_len) cap(fl, *seq_len);
    }

    // Length capping (Brown & Langmead 2026). Every interval is cut into pieces
    // of at most L = CAP * ceil(n / R) text positions of its sequence; Phi^{-1}
    // stays linear on each piece. The sweep takes every Phi step from a distinct
    // row (each root-record row is located once), so a query's fast-forward is
    // charged to the run starts inside its piece's image: at most L per start,
    // (R + n / L) * L <= (CAP + 1) n + R in total, against R times the longest
    // interval uncapped. At most R / CAP pieces are added. O(R + pieces) time;
    // dest_run and head move to the pieces.
    static constexpr uint64_t CAP = 8;
    static uint64_t capLength(const std::vector<uint32_t>& seq_len, size_t R) {
        uint64_t n = 0;
        for (uint32_t l : seq_len) n += l;
        return CAP * std::max<uint64_t>(1, (n + R - 1) / R);
    }
    // Number of intervals after capping, from the tails alone.
    static size_t cappedSize(const gbwt::FastLocate& fl, const std::vector<uint32_t>& seq_len) {
        const size_t R = fl.last_to_run.size();
        if (R == 0) return 0;
        const uint64_t L = capLength(seq_len, R);
        size_t total = 0;
        uint64_t prev = 0; bool have = false;
        auto add = [&](uint64_t t, uint64_t next) {
            const uint64_t seq_end = fl.pack(fl.seqId(t), seq_len[fl.seqId(t)]);
            const uint64_t ve = std::min<uint64_t>(next, seq_end);
            total += std::max<uint64_t>(1, (ve - t + L - 1) / L);
        };
        for (auto it = fl.last.one_begin(); it != fl.last.one_end(); ++it) {
            if (have) add(prev, it->second);
            prev = it->second; have = true;
        }
        add(prev, std::numeric_limits<uint64_t>::max());
        return total;
    }
    void cap(const gbwt::FastLocate& fl, const std::vector<uint32_t>& seq_len) {
        const size_t R = rows.size();
        if (R == 0) return;
        const uint64_t L = capLength(seq_len, R);
        std::vector<uint32_t> first(R + 1, 0);           // index of each interval's first piece
        for (size_t k = 0; k < R; ++k) {
            const uint64_t t = rows[k].start;
            const uint64_t seq_end = fl.pack(fl.seqId(t), seq_len[fl.seqId(t)]);
            const uint64_t ve = (k + 1 < R) ? std::min<uint64_t>(rows[k + 1].start, seq_end) : seq_end;
            const uint64_t pieces = std::max<uint64_t>(1, (ve - t + L - 1) / L);
            if (first[k] + pieces >= std::numeric_limits<uint32_t>::max())
                throw std::runtime_error("MovePhi: 2^32 or more pieces");
            first[k + 1] = (uint32_t)(first[k] + pieces);
        }
        const size_t R2 = first[R];
        if (R2 == R) return;
        // Expand in place from the back; dest_run keeps the old interval for now.
        rows.resize(R2);
        for (size_t k = R; k-- > 0; ) {
            const Row o = rows[k];
            for (uint32_t i = 0; first[k] + i < first[k + 1]; ++i)
                rows[first[k] + i] = Row{o.start + i * L, o.dest + i * L, o.dest_run};
        }
        // dest_run of each piece: from the first piece of the old interval that
        // held the unsplit dest, forward. Over one interval's pieces the scan is
        // monotone, and the rows it passes lie in that interval's image.
        for (size_t k = 0; k < R; ++k) {
            uint32_t j = first[rows[first[k]].dest_run];
            for (uint32_t i = first[k]; i < first[k + 1]; ++i) {
                while (j + 1 < R2 && rows[j + 1].start <= rows[i].dest) ++j;
                rows[i].dest_run = j;
            }
        }
        // Head intervals: the piece of the old interval holding each head sample
        // (in place unless the wider indices need more bits).
        if (sdsl::bits::length(R2 - 1) > head.width()) {
            sdsl::int_vector<0> h(head.size(), 0, sdsl::bits::length(R2 - 1));
            for (size_t j = 0; j < head.size(); ++j) h[j] = head[j];
            head.swap(h);
        }
        for (size_t j = 0; j < head.size(); ++j) {
            const size_t m = head[j], f = first[m], np = first[m + 1] - f;
            head[j] = f + std::min<uint64_t>((fl.samples[j] - rows[f].start) / L, np - 1);
        }
    }

    // Interval holding text position p (predecessor tail), by binary search.
    // Reference only: a binary-search check for `head` and every dest_run.
    uint32_t find(uint64_t p) const {
        auto it = std::upper_bound(rows.begin(), rows.end(), p,
                                   [](uint64_t x, const Row& r) { return x < r.start; });
        return it == rows.begin() ? 0 : (uint32_t)((it - rows.begin()) - 1);
    }

    // SA[i] -> SA[i+1] for p = SA[i] in interval k; k is updated to the new interval.
    uint64_t next(uint64_t p, uint32_t& k) const {
        const Row& r = rows[k];
        const uint64_t q = r.dest + (p - r.start);
        uint32_t j = r.dest_run;
        while (j + 1 < rows.size() && rows[j + 1].start <= q) ++j;
        MOVE_FF(PHI_NEXT, j - r.dest_run);
        k = j;
        return q;
    }
};
