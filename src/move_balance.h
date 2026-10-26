/*
 * Balancing a move structure in O(r) time (Brown, Sanaullah, Zhang & Langmead,
 * "Optimal-Time Move Structure Construction", 2026; the balancing criterion is
 * Nishimoto & Tabei's, ICALP 2021).
 *
 * A permutation pi of [0, N) is given as r input intervals I_j = [P[j], P[j+1])
 * (P sorted, P[r] = N), each mapped contiguously to its output interval
 * [Ppi[j], Ppi[j] + |I_j|). The structure is alpha-balanced when every output
 * interval holds fewer than 2 alpha input starts strictly inside it, so a move
 * query (read one row, step forward from the row holding its image start) takes
 * O(alpha) steps. Here pi^{-1} is balanced as well, which is what lets every
 * split find its place in the other list without a search.
 *
 * Two linked lists hold the input (L_P) and output (L_Q) boundaries in position
 * order, with tau / tau^{-1} links between matching intervals and cross-list
 * predecessor pointers (the interval of the other list holding a node's start).
 * A sweep position t moves left to right. Every interval starting before t has
 * fewer than 2 alpha boundaries of the other list inside it up to t, and every
 * node at or before t has a correct predecessor pointer. Each step takes the
 * interval holding t with the smaller start (the longer on a tie), splits it at
 * its (alpha + 1)-th inside boundary if it is heavy, and mirrors the split into
 * the other list through tau. A mirrored split that lands before t is repaired
 * on the spot, in O(alpha) because that region is balanced, and may cascade;
 * one that lands after t is left for the sweep. Every step inserts one node per
 * list, so the work is O(alpha) per node and the result has at most
 * (alpha + 1) r / (alpha - 1) intervals.
 */
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

template <class Pos = uint64_t>
struct MoveBalancerT {
    static constexpr uint32_t NIL = std::numeric_limits<uint32_t>::max();
    struct List {
        std::vector<Pos>      idx;    // boundary position
        std::vector<uint32_t> next;   // next node in position order, NIL at the end
        std::vector<uint32_t> tau;    // matching interval in the other list
        std::vector<uint32_t> pred;   // interval of the other list holding idx (valid up to t)
        uint32_t add(Pos i, uint32_t nx) {
            idx.push_back(i); next.push_back(nx); tau.push_back(NIL); pred.push_back(NIL);
            return (uint32_t)(idx.size() - 1);
        }
        void reserve(size_t n) { idx.reserve(n); next.reserve(n); tau.reserve(n); pred.reserve(n); }
    };
    List L[2];                        // L[0] = inputs (P), L[1] = outputs (Q)
    Pos N = 0;
    uint32_t alpha = 2;
    size_t r = 0;                     // original intervals
    uint64_t steps = 0;               // list nodes visited, for checking linearity

    Pos end(int s, uint32_t v) const { return L[s].next[v] == NIL ? N : L[s].idx[L[s].next[v]]; }

    // Loading: the r input starts in order, then the r output starts in output
    // order, each with the input interval it belongs to.
    void reserve(size_t n, size_t extra) {
        for (int s = 0; s < 2; ++s) { L[s] = List{}; L[s].reserve(n + extra); }
        r = 0;
    }
    void addInput(Pos start) { L[0].add(start, NIL); ++r; }
    void addOutput(Pos start, uint32_t j) {
        const uint32_t k = L[1].add(start, NIL);
        L[1].tau[k] = j; L[0].tau[j] = k;
    }

    // Boundaries of list o = 1 - s strictly inside node v of list s, below `lim`.
    // Returns their number, capped at 2 alpha, and the (alpha + 1)-th of them.
    uint32_t weight(int s, uint32_t v, Pos lim, uint32_t& kth) {
        const int o = 1 - s;
        const Pos a = L[s].idx[v], b = std::min(end(s, v), lim);
        uint32_t w = 0; kth = NIL;
        for (uint32_t u = L[s].pred[v]; u != NIL && L[o].idx[u] < b; u = L[o].next[u]) {
            ++steps;
            if (L[o].idx[u] <= a) continue;
            if (++w == alpha + 1) kth = u;
            if (w >= 2 * alpha) break;
        }
        return w;
    }

    // Split node v of list s at the boundary u of the other list (strictly inside
    // v), and mirror the split into the other list through tau. t is the sweep.
    void split(int s, uint32_t v, uint32_t u, Pos t) {
        const int o = 1 - s;
        const Pos x = L[o].idx[u];
        if (L[s].idx.size() >= NIL - 1) throw std::runtime_error("MoveBalancer: too many intervals");
        const uint32_t vx = L[s].add(x, L[s].next[v]);
        L[s].next[v] = vx;
        L[s].pred[vx] = u;                               // x is a boundary of list o
        const uint32_t m = L[s].tau[v];                  // matching interval in list o
        const Pos xm = L[o].idx[m] + (x - L[s].idx[v]);
        const uint32_t mx = L[o].add(xm, L[o].next[m]);
        L[o].next[m] = mx;
        L[s].tau[vx] = mx; L[o].tau[mx] = vx;
        if (x < t) {
            // Nodes of list o inside [x, t) now lie in vx.
            for (uint32_t w = L[s].pred[vx]; w != NIL && L[o].idx[w] < std::min(end(s, vx), t); w = L[o].next[w]) {
                ++steps;
                if (L[o].idx[w] >= x) L[o].pred[w] = vx;
            }
        }
        if (xm >= t) return;                             // repaired later by the sweep
        // xm lies in the balanced region: find the interval of list s holding it.
        uint32_t h = L[o].pred[m];
        while (L[s].next[h] != NIL && L[s].idx[L[s].next[h]] <= xm) { h = L[s].next[h]; ++steps; }
        L[o].pred[mx] = h;
        // Nodes of list s inside [xm, t) that lay in m now lie in mx.
        for (uint32_t w = h; w != NIL && L[s].idx[w] < std::min(end(o, mx), t); w = L[s].next[w]) {
            ++steps;
            if (L[s].idx[w] >= xm) L[s].pred[w] = mx;
        }
        // xm is a new boundary strictly inside h unless it starts h; h may now be heavy.
        if (L[s].idx[h] == xm) return;
        uint32_t kth;
        if (weight(s, h, t, kth) >= 2 * alpha) split(s, h, kth, t);
    }

    // Balance the loaded intervals of a permutation of [0, n).
    void balance(Pos n, uint32_t a) {
        N = n; alpha = a;
        if (r == 0) return;
        if (L[1].idx.size() != r) throw std::runtime_error("MoveBalancer: inputs and outputs differ in number");
        for (int s = 0; s < 2; ++s)
            for (size_t i = 0; i + 1 < r; ++i) L[s].next[i] = (uint32_t)(i + 1);
        // Initial predecessor pointers: two-finger walk.
        {
            uint32_t q = 0;
            for (uint32_t p = 0; p < r; ++p) {
                while (L[1].next[q] != NIL && L[1].idx[L[1].next[q]] <= L[0].idx[p]) q = L[1].next[q];
                L[0].pred[p] = q;
            }
            uint32_t p = 0;
            for (uint32_t k = 0; k < r; ++k) {
                while (L[0].next[p] != NIL && L[0].idx[L[0].next[p]] <= L[1].idx[k]) p = L[0].next[p];
                L[1].pred[k] = p;
            }
        }
        // Sweep.
        Pos t = 0;
        uint32_t cur[2] = {0, 0};                        // intervals holding t
        while (t < N) {
            const Pos a0 = L[0].idx[cur[0]], a1 = L[1].idx[cur[1]];
            int s;
            if (a0 != a1) s = (a0 < a1) ? 0 : 1;
            else s = (end(0, cur[0]) >= end(1, cur[1])) ? 0 : 1;
            const uint32_t v = cur[s];
            uint32_t kth;
            Pos nt;
            if (a0 == a1 && end(0, cur[0]) == end(1, cur[1])) nt = end(0, cur[0]);   // same interval
            else if (weight(s, v, N, kth) >= 2 * alpha) {
                nt = L[1 - s].idx[kth];
                if (nt <= t) throw std::logic_error("MoveBalancer: sweep did not advance");
                split(s, v, kth, t);
            } else nt = end(s, v);
            t = nt;
            // Advance both fingers to the intervals holding t, fixing predecessor
            // pointers of every node passed (ties: both lists step together).
            for (;;) {
                const uint32_t n0 = L[0].next[cur[0]], n1 = L[1].next[cur[1]];
                const bool m0 = n0 != NIL && L[0].idx[n0] <= t, m1 = n1 != NIL && L[1].idx[n1] <= t;
                if (!m0 && !m1) break;
                ++steps;
                if (m0 && m1 && L[0].idx[n0] == L[1].idx[n1]) {
                    cur[0] = n0; cur[1] = n1; L[0].pred[n0] = n1; L[1].pred[n1] = n0;
                } else if (m0 && (!m1 || L[0].idx[n0] < L[1].idx[n1])) {
                    cur[0] = n0; L[0].pred[n0] = cur[1];
                } else {
                    cur[1] = n1; L[1].pred[n1] = cur[0];
                }
            }
        }
    }

    size_t size() const { return L[0].idx.size(); }

    // The balanced intervals in input order: f(input start, output start,
    // original interval). Nodes below r are the original intervals, in order;
    // every later node is a piece of the original interval before it.
    template <class F> void forEach(F f) const {
        uint32_t j = 0;
        for (uint32_t p = 0; p != NIL; p = L[0].next[p]) {
            if (p < r) j = p;
            f(L[0].idx[p], L[1].idx[L[0].tau[p]], j);
        }
    }
    void clear() { for (int s = 0; s < 2; ++s) L[s] = List{}; r = 0; }
};
