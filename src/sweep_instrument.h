/*
 * sweep_instrument.h  —  optional measurement scaffolding for the dag-free engine.
 *
 * Three independent, flag-gated frameworks. All of them compile to nothing unless their flag is
 * set, so the shipped binary is unaffected:
 *   -DSWEEP_STATS     traversal counters   (SSTAT)
 *   -DSWEEP_ALLOC     allocation ownership (SALLOC_INIT / SALLOC_SCOPE / SALLOC_CAP)
 *   -DSWEEP_PROFILE   phase timers         (PPHASE / PPROF_INIT)
 *
 * INCLUDE THIS FROM EXACTLY ONE TRANSLATION UNIT. Under -DSWEEP_ALLOC it defines global
 * operator new/delete; a second includer would be a duplicate-symbol link error.
 *
 * It lives here rather than at the top of sweep_dagfree.cpp so that file opens with the
 * algorithm. Keep it: these are the reproducibility record behind the measurements in
 * docs/RUN_CELLS_PLAN.md, and instrumentation that is never run rots.
 */
#ifndef SWEEP_INSTRUMENT_H
#define SWEEP_INSTRUMENT_H

// No includes of its own: this header is included immediately after sweep_common.h, which
// already provides <cstdint>, <cstdlib>, <algorithm>, <iostream> and <new>. Pulling them in
// again here changes nothing semantically but does change libstdc++ inlining decisions, and
// therefore codegen — see docs/RUN_CELLS_PLAN.md section 9.

// Optional traversal counters. These are compiled out unless SWEEP_STATS is set.
#ifdef SWEEP_STATS
namespace sweepstats {
struct Counters {
    uint64_t states = 0;          // report() calls on right-maximal states
    uint64_t bwd_ext = 0;         // bdExtendBackward calls (leftExt)
    uint64_t left_children = 0;   // sum of d_W = nonempty left extensions over visited states
    uint64_t rmax_tests = 0;      // right-maximality tests
    uint64_t run_scans = 0;       // runs decoded when building a root's run index
    uint64_t probes = 0;          // binary-search probes into that index (columns + cells)
    uint64_t cells = 0;           // cells constructed
    uint64_t cells_mat = 0;       // cells materialized (needed[] gate passed)
    uint64_t pieces = 0;          // run pieces in materialized cells
    uint64_t keys = 0;            // rows located (Phi walk output)
    uint64_t pairs = 0;           // candidate pairs pushed to localbuf
    uint64_t seeds = 0;           // seeds emitted after last-write-wins dedup
    uint64_t rec_ctor = 0;        // CompressedRecord constructions (all transient)
    uint64_t rec_edges = 0;       // sum of outdegree over those constructions
    uint64_t bwd_scans = 0;       // runs touched by leftExt's single range-distinct pass
    uint64_t bwd_scans_1x1 = 0;   // counterfactual cost of one bdLF scan per edge
    uint64_t bwd_scans_cond = 0;  // range-distinct work spent at conduit states
    uint64_t conduits = 0;        // states that were conduits
    uint64_t chains = 0, chain_len = 0, chain_max = 0;  // maximal conduit chains: count, total
    // Output that depends on treating stored endmarker rows as private path boundaries.
    uint64_t b_states = 0;        // emitting states with NO real left extension at all
    uint64_t b_pairs = 0;         // candidate pairs whose BOTH left contexts are the endmarker
    // Typing direction flips, and the path depth a from-scratch retype would rescan at each.
    uint64_t type_flips = 0, type_flip_depth = 0;
    uint64_t range_sampled = 0;   // left extensions that took the sampled long-range path
    ~Counters() {
        std::cerr << "[stats] states=" << states << " rmax_tests=" << rmax_tests
                  << " bwd_ext=" << bwd_ext << " left_children=" << left_children
                  << " run_scans=" << run_scans << " probes=" << probes << "\n"
                  << "[stats] cells=" << cells << " materialized=" << cells_mat
                  << " pieces=" << pieces << " keys=" << keys
                  << " pairs=" << pairs << " seeds=" << seeds << "\n"
                  << "[stats] rec_ctor=" << rec_ctor << " rec_edges=" << rec_edges
                  << " bwd_scans=" << bwd_scans << " bwd_scans_per_edge=" << bwd_scans_1x1 << "\n"
                  << "[stats] conduits=" << conduits << " bwd_scans_at_conduits=" << bwd_scans_cond
                  << " chains=" << chains << " chain_len_total=" << chain_len
                  << " chain_max=" << chain_max << "\n"
                  << "[stats] b_states=" << b_states << " b_pairs=" << b_pairs << "\n"
                  << "[stats] type_flips=" << type_flips << " type_flip_depth=" << type_flip_depth
                  << " range_sampled=" << range_sampled
                  << std::endl;
    }
};
inline Counters& c() { static Counters instance; return instance; }
}
#define SSTAT(field, n) (sweepstats::c().field += (n))
#else
#define SSTAT(field, n) ((void)0)
#endif

// Optional allocation attribution. The scoped buckets identify which local
// owner caused an operator-new allocation; direct malloc users remain in OTHER.
#ifdef SWEEP_ALLOC
#include <cstdlib>
namespace sweepalloc {
enum Bucket { OTHER = 0, BWD_REC, EXT_VEC, CELL_SCRATCH, PIECE_KEYS, LOCALBUF, FRAMES,
              ROOT_REC, NEEDED_VEC, SEQS_VEC, BD_EXT, NBUCKETS };
inline const char* name(int b) {
    static const char* n[NBUCKETS] = {"other/unattributed", "backward CompressedRecord",
        "leftExt result vector", "build_cells scratch", "pieceKeys occurrence vector",
        "localbuf (pairing)", "DFS frame stack", "root CompressedRecord",
        "report needed[] vector", "report seqs[] vector-of-vectors",
        "bdExtendBackward internals"};
    return n[b];
}
struct Counter { uint64_t calls = 0, bytes = 0, maxreq = 0; };
inline Counter* stats() { static Counter c[NBUCKETS]; return c; }
inline int& cur() { static int c = OTHER; return c; }
struct Scope {
    int prev;
    explicit Scope(int b) : prev(cur()) { cur() = b; }
    ~Scope() { cur() = prev; }
};
// Structure-level capacity peaks (bytes), updated once per state — not on the allocation path.
struct Caps { size_t cells = 0, root_runs = 0, localbuf = 0, stk = 0; };
inline Caps& caps() { static Caps c; return c; }
struct Report {
    ~Report() {
        uint64_t tc = 0, tb = 0;
        for (int i = 0; i < NBUCKETS; ++i) { tc += stats()[i].calls; tb += stats()[i].bytes; }
        std::cerr << "[alloc] operator-new attribution (calls / bytes / largest / % of calls)\n";
        for (int i = 0; i < NBUCKETS; ++i) {
            const auto& s = stats()[i];
            if (!s.calls) continue;
            std::cerr << "[alloc]   " << name(i) << " = " << s.calls << " / " << s.bytes
                      << " / " << s.maxreq << " / "
                      << (tc ? 100.0 * (double)s.calls / (double)tc : 0.0) << "%\n";
        }
        std::cerr << "[alloc]   TOTAL = " << tc << " calls / " << tb << " bytes\n";
        std::cerr << "[alloc] peak capacity (bytes): cells=" << caps().cells
                  << " root_runs=" << caps().root_runs << " localbuf=" << caps().localbuf
                  << " stk=" << caps().stk << std::endl;
    }
};
inline Report& report_at_exit() { static Report r; return r; }
}
#define SALLOC_INIT() (void)sweepalloc::report_at_exit()
#define SALLOC_SCOPE(b) sweepalloc::Scope _as_(sweepalloc::b)
#define SALLOC_CAP(field, v) \
    (sweepalloc::caps().field = std::max(sweepalloc::caps().field, (v)))
void* operator new(size_t n) {
    auto& c = sweepalloc::stats()[sweepalloc::cur()];
    ++c.calls; c.bytes += n; if (n > c.maxreq) c.maxreq = n;
    void* p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    return p;
}
void* operator new[](size_t n) { return operator new(n); }
void  operator delete(void* p) noexcept { std::free(p); }
void  operator delete[](void* p) noexcept { std::free(p); }
void  operator delete(void* p, size_t) noexcept { std::free(p); }
void  operator delete[](void* p, size_t) noexcept { std::free(p); }
#else
#define SALLOC_INIT() ((void)0)
#define SALLOC_SCOPE(b) ((void)0)
#define SALLOC_CAP(field, v) ((void)0)
#endif

// Optional phase timers. Use them for relative attribution only; each probe has
// nonzero cost and the code is compiled out of normal builds.
#ifdef SWEEP_PROFILE
#include <x86intrin.h>
namespace sweepprof {
enum Phase { TYPING = 0, LEFT_EXT, LEFT_REC, LEFT_BDLF, BUILD_CELLS, NEEDED_GATE, MATERIALIZE_PHI, PAIRING,
             STABLE_SORT, DRAIN_EMIT, CHILD_RMAX, TOTAL, NPHASES };
inline const char* name(int p) {
    static const char* n[NPHASES] = {"typing", "leftExt", "  ^ record ctor", "  ^ bdLF scans",
        "build_cells", "needed_gate",
        "materialize_phi", "pairing", "stable_sort", "drain_emit", "child_rmax", "TOTAL"};
    return n[p];
}
inline uint64_t* acc() { static uint64_t a[NPHASES] = {0}; return a; }
struct Timer {
    int phase; uint64_t t0;
    explicit Timer(int p) : phase(p), t0(__rdtsc()) {}
    ~Timer() { acc()[phase] += __rdtsc() - t0; }
};
struct Report {
    ~Report() {
        const double tot = (double)acc()[TOTAL];
        uint64_t sum = 0;
        std::cerr << "[prof] phase breakdown (rdtsc cycles, % of TOTAL)\n";
        for (int i = 0; i < TOTAL; ++i) {
            if (i != LEFT_REC && i != LEFT_BDLF) sum += acc()[i];   // nested inside LEFT_EXT
            std::cerr << "[prof]   " << name(i) << " = " << acc()[i] << "  ("
                      << (tot ? 100.0 * (double)acc()[i] / tot : 0.0) << "%)\n";
        }
        std::cerr << "[prof]   TOTAL = " << acc()[TOTAL] << "  (100%)\n"
                  << "[prof]   accounted = " << sum << "  ("
                  << (tot ? 100.0 * (double)sum / tot : 0.0)
                  << "% of TOTAL; the remainder is loop overhead + the probes themselves)"
                  << std::endl;
    }
};
inline Report& report_at_exit() { static Report r; return r; }
}
#define PPHASE(p) sweepprof::Timer _pp_(sweepprof::p)
#define PPROF_INIT() (void)sweepprof::report_at_exit()
#else
#define PPHASE(p) ((void)0)
#define PPROF_INIT() ((void)0)
#endif
#endif  // SWEEP_INSTRUMENT_H
