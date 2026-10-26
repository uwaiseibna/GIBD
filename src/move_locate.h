/*
 * FastLocate construction over the move structure.
 *
 * FastLocate stores the text position of every logical run head (samples) and
 * run tail (last, last_to_run). gbwt::FastLocate(const GBWT&) finds them by
 * walking every path with CompressedRecord::LF, which locates the record,
 * decodes its header and scans its runs from the record start at each of the n
 * steps: O(n * r_max). The same walk over MoveLF reads one row and seeks
 * forward from the row's dest_run per step, and it yields the samples already
 * in their stored order, so the library's two sorts are not needed either.
 * The result is identical to the library's structure.
 *
 * Optionally it also returns, for every logical run, the index of the run tail
 * that precedes the run head in text order: the Move-Phi interval holding the
 * head sample (MovePhi::head). The walk meets heads and tails of a sequence in
 * walk order, so this costs O(1) per head and no predecessor search.
 */
#pragma once

#include "move_lf.h"

#include <gbwt/fast_locate.h>
#include <sdsl/sd_vector.hpp>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

inline void buildFastLocate(const gbwt::GBWT& index, const MoveLF& mv, gbwt::FastLocate& fl,
                            sdsl::int_vector<0>* head_row = nullptr,
                            std::vector<uint32_t>* seq_len_out = nullptr) {
    using gbwt::size_type;
    fl.setGBWT(index);
    if (index.empty()) return;
    const size_type eff = index.effective(), nseq = index.sequences();

    // Logical runs number the samples. An endmarker-destined run is one
    // logical run per row, any other run is one logical run.
    std::vector<uint32_t> glog(mv.rows.size());   // global logical id of each row's first row
    fl.comp_to_run.resize(eff);
    size_type total = 0;
    for (size_type c = 0; c < eff; ++c) {
        fl.comp_to_run[c] = total;
        for (uint32_t k = mv.first_run[c]; k < mv.first_run[c + 1]; ++k) {
            if (total >= std::numeric_limits<uint32_t>::max())
                throw std::runtime_error("buildFastLocate: 2^32 or more logical runs");
            if (mv.cont[k] && mv.rows[k].succ != gbwt::ENDMARKER) {
                glog[k] = glog[k - 1];          // a balancing piece continues its run
                continue;
            }
            glog[k] = (uint32_t)total;
            total += (mv.rows[k].succ == gbwt::ENDMARKER) ? mv.rows[k].len : 1;
        }
    }
    sdsl::util::bit_compress(fl.comp_to_run);

    // Row i of the endmarker record starts sequence i. One pass over its move
    // rows gives every start and the move row holding it.
    std::vector<gbwt::edge_type> start(nseq, gbwt::invalid_edge());
    std::vector<uint32_t> start_row(nseq, 0);
    for (uint32_t k = mv.first_run[0]; k < mv.first_run[1]; ++k) {
        const MoveLF::Row& R = mv.rows[k];
        uint32_t kk = R.dest_run;
        for (size_type i = R.start; i < (size_type)R.start + R.len; ++i) {
            if (R.succ == gbwt::ENDMARKER) { start[i] = {gbwt::ENDMARKER, 0}; continue; }
            start[i] = {R.succ, R.dest + (i - R.start)};
            kk = mv.seek(kk, start[i].second);
            start_row[i] = kk;
        }
    }
    // Logical runs of the endmarker record, numbered as the library does.
    std::vector<uint32_t> erun(nseq, 0);
    {
        uint32_t id = 0;
        gbwt::node_type prev = start[0].first;
        for (size_type i = 1; i < nseq; ++i) {
            const gbwt::node_type cur = start[i].first;
            if (cur == gbwt::ENDMARKER || cur != prev) { id++; prev = cur; }
            erun[i] = id;
        }
    }

    // Walk every sequence. Offsets are counted from the path start here and
    // turned into distances to the path end once the length is known. For the
    // head intervals, each head also notes the first tail at or after it in
    // walk order: in text order, the largest tail <= the head on its sequence.
    const bool want_head = (head_row != nullptr);
    constexpr uint32_t NO_TAIL = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> head_seq(total), head_walk(total), seq_len(nseq);
    std::vector<uint32_t> head_tl(want_head ? total : 0);   // walk-order index into tails[seq]
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> tails(nseq);   // (walk offset, run)
    size_type max_length = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(max : max_length)
    for (size_type i = 0; i < nseq; ++i) {
        std::vector<std::pair<uint32_t, uint32_t>>& tl = tails[i];
        std::vector<uint32_t> pend;   // heads of this sequence waiting for a tail
        auto head = [&](uint32_t run, size_type w) {
            head_seq[run] = (uint32_t)i; head_walk[run] = (uint32_t)w;
            if (want_head) pend.push_back(run);
        };
        auto tail = [&](uint32_t run, size_type w) {
            tl.push_back({(uint32_t)w, run});
            if (want_head) {
                for (uint32_t h : pend) head_tl[h] = (uint32_t)(tl.size() - 1);
                pend.clear();
            }
        };
        const uint32_t er = erun[i];   // the endmarker's logical runs come first
        if (i == 0 || er != erun[i - 1]) head(er, 0);
        if (i + 1 >= nseq || er != erun[i + 1]) tail(er, 0);
        gbwt::edge_type curr = start[i];
        uint32_t k = start_row[i];
        size_type walk = 1;
        while (curr.first != gbwt::ENDMARKER) {
            const MoveLF::Row& R = mv.rows[k];
            const size_type p = curr.second;
            const bool to_end = (R.succ == gbwt::ENDMARKER);
            const uint32_t run = glog[k] + (to_end ? (uint32_t)(p - R.start) : 0u);
            // Heads and tails of the run, not of a balancing piece.
            if (to_end || (p == R.start && !mv.cont[k])) head(run, walk);
            if (to_end || (p + 1 == (size_type)R.start + R.len
                           && !(k + 1 < mv.rows.size() && mv.cont[k + 1]))) tail(run, walk);
            ++walk;
            if (to_end) break;
            curr = {R.succ, R.dest + (p - R.start)};
            k = mv.seek(R.dest_run, curr.second);
            MOVE_FF(FL_WALK, k - R.dest_run);
        }
        for (uint32_t h : pend) head_tl[h] = NO_TAIL;   // its predecessor is in an earlier sequence
        seq_len[i] = (uint32_t)walk;
        max_length = std::max(max_length, walk);
    }
    if (max_length >= std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("buildFastLocate: sequence of 2^32 or more nodes");
    fl.header.max_length = max_length;
    if (seq_len_out) *seq_len_out = seq_len;   // text positions of sequence s: offsets [0, seq_len[s])

    // Head samples, by run id.
    fl.samples.width(sdsl::bits::length(fl.pack(nseq - 1, max_length - 1)));
    fl.samples.resize(total);
    for (size_type r = 0; r < total; ++r) {
        const uint32_t s = head_seq[r];
        fl.samples[r] = fl.pack(s, seq_len[s] - 1 - head_walk[r]);
    }
    std::vector<uint32_t>().swap(head_walk);

    // Head intervals. Tails are numbered in text order, so the tail at
    // walk-order index j of sequence s is tail base[s] + |tails[s]| - 1 - j.
    // A head with no tail after it on its sequence follows the last tail of
    // the earlier sequences (or interval 0 if there is none).
    if (want_head) {
        std::vector<size_type> base(nseq + 1, 0);
        for (size_type i = 0; i < nseq; ++i) base[i + 1] = base[i] + tails[i].size();
        head_row->width(std::max<uint8_t>(1, sdsl::bits::length(total - 1)));
        head_row->resize(total);
        for (size_type r = 0; r < total; ++r) {
            const uint32_t s = head_seq[r], j = head_tl[r];
            (*head_row)[r] = (j != NO_TAIL) ? base[s] + tails[s].size() - 1 - j
                                            : (base[s] ? base[s] - 1 : 0);
        }
        std::vector<uint32_t>().swap(head_tl);
    }
    std::vector<uint32_t>().swap(head_seq);

    // Tail samples, by text position: sequences in order, and within one
    // sequence the reverse of walk order.
    sdsl::sd_vector_builder builder(nseq * max_length, total);
    fl.last_to_run.width(sdsl::bits::length(total - 1));
    fl.last_to_run.resize(total);
    size_type t = 0;
    for (size_type i = 0; i < nseq; ++i) {
        for (auto it = tails[i].rbegin(); it != tails[i].rend(); ++it) {
            builder.set_unsafe(fl.pack(i, seq_len[i] - 1 - it->first));
            fl.last_to_run[t++] = it->second;
        }
        std::vector<std::pair<uint32_t, uint32_t>>().swap(tails[i]);
    }
    fl.last = sdsl::sd_vector<>(builder);
}
