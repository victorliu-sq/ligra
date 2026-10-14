// X-Blossom-Pro (maximum matching in general graphs) written as a Ligra app.
//
// The algorithm is XB-Pro's (X-Blossom repo, src/xblsm_cpu/xb_cpu.cu with
// XB_REUSE=1): search phases grow a forest of alternating trees from the
// exposed nodes, each phase repeats
//   1. augment  -- an edge between even nodes of two different trees is an
//                  augmenting path; claim both trees and record the path
//   2. expand   -- an even node's unmatched-tree neighbour w joins the tree as
//                  an odd node and its mate M[w] as an even node
//   3. blossom  -- an edge between two even nodes of the same tree closes a
//                  blossom; its odd nodes become even
// until paths are found (flip them, next phase) or the forest stops growing
// (maximum). Trees whose root is still exposed are kept across phases.
//
// What Ligra does here: the frontier of each step is a vertexSubset and the
// three steps are edgeMap functors (sparse edgeMap, no_dense: the pull
// direction has no meaning for alternating trees), so the framework decides
// how frontier edges are spread over workers. What stays X-Blossom's own: the
// per-node path table, tree membership and the per-tree / per-match /
// per-blossom claims. One semantic difference from the hand-written XB-Pro: a
// functor sees one edge at a time, so "this node's tree is taken, skip its
// remaining edges" becomes "skip this edge".
//
// Output (the GACGE measurement format):
//   XBConfig: ...                 once
//   XBRound: index runtime_s matching_size valid     per timed round
//   XBResult: status matching_size valid rounds mean_runtime_s
// Each Ligra round is one full matching from scratch; the first -warmup rounds
// are untimed warm-ups, so -rounds is warm-up + timed rounds.
//   CILK_NWORKERS=<n> ./XBlossomPro_cilk -s -rounds 11 -warmup 1 -dataset Amazon <adj file>
// (built with OpenCilk: make XBlossomPro_cilk; without it the app runs serially)
#include "ligra.h"

// Compile-time switches (as X-Blossom's XB_REUSE / XB_COUNTERS):
//   XBP_REUSE    1 = keep trees whose root is still exposed (XB-Pro-Ligra),
//                0 = rebuild every tree each search phase (XB-Ligra)
//   XBP_COUNTERS 1 = count reused / reset tree nodes per round (XBCounters lines)
//   XBP_OPT      1 = (default) no serial or contended parts beyond the
//                algorithm's own: the two trees of an augmenting path are
//                claimed in id order (no mutual back-off, so no serial rescan
//                of the edges given up), augmenting paths go to one flat buffer per
//                worker, path and blossom walks reuse a per-worker scratch
//                buffer, a blossom is collected front to back and reversed
//                once, and the per-node path table reservation and the frontier
//                copies run in parallel (1.6x faster at 48 workers, GM over the
//                ten SNAP graphs, measure/xb-cpu-opt);
//                0 = XB-Pro's own design (XBlossomPro_xbpro_cilk): a mutex
//                around one shared path list, vectors allocated per path and
//                per blossom edge, serial copies
#ifndef XBP_REUSE
#define XBP_REUSE 1
#endif
#ifndef XBP_COUNTERS
#define XBP_COUNTERS 0
#endif
#ifndef XBP_OPT
#define XBP_OPT 1
#endif
// XBP_TRIPLET 1 = (default with XBP_OPT, since 2026-10-13) XB-Gunrock's
//              compact path table and cycle-free blossom handling
//              on the CPU: every node keeps one record of three node ids
//              (12 B) instead of a vector of its path towards the root
//              (XB-Pro reserves 100 ids per node):
//                root or odd node                 (-1, -1, -1)
//                x even by expand through v       (w, v, -1)    w = x's mate
//                u even by a blossom with closing (s, o, base)
//                edge (s, o) and base
//              A blossom is never written out: its base is found by walking the
//              records' anchors (v, or base), the odd nodes on both chains are
//              claimed and become even after the step with record (s, o, base);
//              an augmenting edge is kept as (v, w) and its paths are flipped by
//              unfolding the records (as XB-Gunrock's TuplePathTable and
//              ClaimOddNodesOnChain, src/gunrock/.../xblossom)
//              0 = a vector of path nodes per node, 100 reserved, and blossoms
//              written out (XB-Pro's path table; XBlossomPro_vector_cilk, and
//              the XB-Pro design build)
#ifndef XBP_TRIPLET
#define XBP_TRIPLET XBP_OPT
#endif
#if XBP_TRIPLET && !XBP_OPT
#error "XBP_TRIPLET builds on XBP_OPT"
#endif

#if defined(OPENCILK)
static inline int WorkerId() { return __cilkrts_get_worker_number(); }
#else
static inline int WorkerId() { return 0; }
#endif

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// GACGE_XB_LOAD_FIRST=1: read a claim flag before its compare-and-swap and
// skip the CAS when the flag is already taken (test-and-test-and-set). Every
// claim has the same outcome; the A/B question (2026-10-14) is whether the
// failed seq_cst CASes on flags other workers own are what makes web graphs
// (2-17 billion blossom closing edges) slow on the CPU.
static const bool kLoadFirst = [] {
  const char *e = std::getenv("GACGE_XB_LOAD_FIRST");
  return e != nullptr && e[0] == '1';
}();
static inline bool TryClaim(std::atomic<int> &flag) {
  if (kLoadFirst && flag.load(std::memory_order_relaxed) != 0) return false;
  int expected = 0;
  return flag.compare_exchange_strong(expected, 1);
}

// Cilkscale (GACGE_CILKSCALE): the work and span of a part of one Compute(),
// one CSV row tagged <name><call> (init, match, check) next to ligra.h's run<r>.
#ifdef GACGE_CILKSCALE
#define XBP_WS_BEGIN(v) wsp_t v = wsp_getworkspan()
#define XBP_WS_END(v, name, call)                                         \
  do {                                                                    \
    char t_[32];                                                          \
    std::snprintf(t_, sizeof t_, "%s%ld", name, static_cast<long>(call)); \
    wsp_dump(wsp_sub(wsp_getworkspan(), v), t_);                          \
  } while (0)
#else
#define XBP_WS_BEGIN(v) (void)0
#define XBP_WS_END(v, name, call) (void)0
#endif

namespace xbp {

// Vectors whose elements start uninitialized (XBP_OPT fills them in parallel):
// std::vector zero-fills new elements serially, on construction and on resize.
template <class T>
struct NoInitAlloc : std::allocator<T> {
  template <class U>
  struct rebind {
    using other = NoInitAlloc<U>;
  };
  NoInitAlloc() = default;
  template <class U>
  NoInitAlloc(const NoInitAlloc<U> &) noexcept {}
  template <class U>
  void construct(U *p) noexcept(std::is_nothrow_default_constructible<U>::value) {
    ::new (static_cast<void *>(p)) U;
  }
  template <class U, class... A>
  void construct(U *p, A &&...a) {
    ::new (static_cast<void *>(p)) U(std::forward<A>(a)...);
  }
};
using Frontier = std::vector<uintE, NoInitAlloc<uintE>>;  // a frontier, filled by parallel copies
using Ints = std::vector<int, NoInitAlloc<int>>;          // a per-node field
using Flags = std::vector<std::atomic<int>, NoInitAlloc<std::atomic<int>>>;

struct State {
  long n = 0;
  Ints M;         // mate, -1 = exposed
  Ints is_even;   // 1 = even node of a tree
  Ints belongs;   // tree (root) of a node, -1 = none
  Ints blossom_to_base;
  std::vector<std::vector<int>> path_table;  // node -> nodes towards its root
  Flags select_tree, select_match, select_blossom;

  std::vector<std::vector<int>> paths;  // augmenting paths of this phase
  std::mutex paths_mutex;
  // XBP_OPT: the augmenting paths per worker, path i of worker t being
  // path_nodes[t][path_start[t][i] .. path_start[t][i + 1] or the end); and a
  // per-worker scratch buffer for one root path or one blossom.
  std::vector<std::vector<int>> path_nodes, path_start, scratch;
  // XBP_TRIPLET: the records, three node ids per node; per worker, the
  // augmenting edges of this phase (v, w pairs) and the odd nodes claimed by
  // this step's blossoms (u, s, o, base).
  int *rec = nullptr;
  std::vector<std::vector<int>> aug_pairs, tasks;
  // Frontier nodes whose scan gave up a cross-tree edge under contention, per worker.
  std::vector<std::vector<uintE>> dropped;

  // Nodes produced by expand / blossom, one buffer per worker (as XB-Pro's
  // thread-local vectors: one shared counter would make every emit contend).
  std::vector<std::vector<uintE>> out;

  // Tree nodes reused / reset at the start of each search phase (XBP_COUNTERS).
  std::atomic<unsigned long long> reused_tree_nodes{0}, reset_tree_nodes{0};

  explicit State(long n_) : n(n_), M(n_), is_even(n_), belongs(n_), blossom_to_base(n_),
                            path_table(XBP_TRIPLET ? 0 : n_), select_tree(n_), select_match(n_), select_blossom(n_),
                            dropped(std::max(1, getWorkers())), out(std::max(1, getWorkers())),
                            path_nodes(XBP_OPT ? std::max(1, getWorkers()) : 0),
                            path_start(XBP_OPT ? std::max(1, getWorkers()) : 0),
                            scratch(XBP_OPT ? std::max(1, getWorkers()) : 0),
                            aug_pairs(XBP_TRIPLET ? std::max(1, getWorkers()) : 0),
                            tasks(XBP_TRIPLET ? std::max(1, getWorkers()) : 0) {
    // Every per-node field in one pass, parallel under XBP_OPT (XB-Pro's
    // design: serial, as std::vector's fill constructor).
    auto init = [this](long x) {
      M[x] = -1;
      is_even[x] = 0;
      belongs[x] = -1;
      blossom_to_base[x] = -1;
      select_tree[x].store(0, std::memory_order_relaxed);
      select_match[x].store(0, std::memory_order_relaxed);
      select_blossom[x].store(0, std::memory_order_relaxed);
    };
    if (XBP_OPT) {
      parallel_for(long x = 0; x < n_; x++) init(x);
    } else {
      for (long x = 0; x < n_; x++) init(x);
    }
    if (XBP_TRIPLET) {
      rec = newA(int, 3 * n_);
      parallel_for(long i = 0; i < 3 * n_; i++) rec[i] = -1;
      return;
    }
    // as XB-Pro: limits reallocation under concurrent readers
    if (XBP_OPT) {
      parallel_for(long i = 0; i < n_; i++) path_table[i].reserve(100);
    } else {
      for (auto &p : path_table) p.reserve(100);
    }
  }

  ~State() { free(rec); }
  State(const State &) = delete;
  State &operator=(const State &) = delete;

  bool have_paths() const {
    if (XBP_TRIPLET) {
      for (auto &p : aug_pairs)
        if (!p.empty()) return true;
      return false;
    }
    if (!XBP_OPT) return !paths.empty();
    for (auto &p : path_start)
      if (!p.empty()) return true;
    return false;
  }

  void emit(int v) { out[WorkerId()].push_back(v); }
};

// ---- XBP_TRIPLET: the records -----------------------------------------------
static inline void SetRec(State &S, int x, int a, int b, int c) {
  int *r = S.rec + 3 * static_cast<size_t>(x);
  r[0] = a;
  r[1] = b;
  r[2] = c;
}
static inline bool RecEmpty(const State &S, int x) { return S.rec[3 * static_cast<size_t>(x)] == -1; }
// the next even node towards the root: v of an expand record, base of a blossom record
static inline int Anchor(const State &S, int x) {
  const int *r = S.rec + 3 * static_cast<size_t>(x);
  return r[2] == -1 ? r[1] : r[2];
}

// The base of the blossom closed by (v, w): the first anchor the two anchor
// chains share (each walks to its root, then switches to the other side's
// start). -1 if the chains never meet (a forest read mid-update). *steps: the
// walk's iterations.
static int FindBaseByAnchor(const State &S, int v, int w, long *steps) {
  int cv = v, cw = w;
  long it = 0;
  while (cv != cw) {
    cv = RecEmpty(S, cv) ? w : Anchor(S, cv);
    cw = RecEmpty(S, cw) ? v : Anchor(S, cw);
    if (++it > 2 * S.n + 2) {
      *steps = it;
      return -1;
    }
  }
  *steps = it;
  return cv;
}

// Claim every node on s's anchor chain down to base that was odd when it was
// recorded (the w of each expand record) and is still odd: one task
// (u, s, o, base) per claim, applied after the step. Returns the records visited.
static long ClaimOddOnChain(State &S, int s, int o, int base) {
  std::vector<int> &tk = S.tasks[WorkerId()];
  int cur = s;
  long visited = 0;
  while (cur != base) {
    const int *r = S.rec + 3 * static_cast<size_t>(cur);
    visited++;
    if (r[0] == -1) return visited;  // the root without meeting the base
    if (r[2] == -1) {
      const int u = r[0];
      if (!S.is_even[u] && TryClaim(S.select_blossom[u])) {
        tk.push_back(u);
        tk.push_back(s);
        tk.push_back(o);
        tk.push_back(base);
      }
    }
    cur = r[2] == -1 ? r[1] : r[2];
  }
  return visited;
}

// Flip the matched pairs of route(x) up to its root. A blossom record
// (s, o, base) met at u stands for (s, o), route(s) before u, and route(o),
// which reaches base and goes on as the rest of this walk: route(s) is
// walked first, route(o) kept on a stack (one entry per nesting level).
static void FlipRoute(State &S, int x) {
  struct Walk {
    int x, stop;  // stop >= 0: end at the expand record whose w is stop (an odd node); -1: the root
  };
  Walk stack[256];
  int depth = 0;
  Walk cur{x, -1};
  unsigned long visited = 0;  // records (GACGE_REGIONS: the walk's length)
  while (true) {
    while (true) {
      const int *r = S.rec + 3 * static_cast<size_t>(cur.x);
      visited++;
      if (r[0] == -1) break;  // the root
      if (r[2] == -1) {
        if (cur.stop >= 0 && r[0] == cur.stop) break;
        S.M[r[0]] = r[1];
        S.M[r[1]] = r[0];
        cur.x = r[1];
      } else {
        if (depth == 256) {
          std::fprintf(stderr, "XBError: blossom nesting deeper than 256 while flipping a path\n");
          std::abort();
        }
        stack[depth++] = Walk{r[1], cur.stop};
        S.M[r[0]] = r[1];
        S.M[r[1]] = r[0];
        cur = Walk{r[0], cur.x};
      }
    }
    if (depth == 0) break;
    cur = stack[--depth];
  }
  gacge_walk(visited);
}

// v, then the nodes of its path table entries up to its tree's root, appended to out.
static void AppendPathToRoot(const State &S, int v, std::vector<int> &out) {
  out.push_back(v);
  int cur = v;
  while (!S.path_table[cur].empty()) {
    out.insert(out.end(), S.path_table[cur].begin(), S.path_table[cur].end());
    cur = out.back();
  }
}

static std::vector<int> PathToRoot(const State &S, int v) {
  std::vector<int> path{v};
  int cur = v;
  while (!S.path_table[cur].empty()) {
    path.insert(path.end(), S.path_table[cur].begin(), S.path_table[cur].end());
    cur = path.back();
  }
  return path;
}

// 1. Augmenting paths.
struct Augment_F {
  static constexpr const char *kPhase = "augment";
  State &S;
  explicit Augment_F(State &s) : S(s) {}
  inline bool cond(uintE w) { return S.is_even[w]; }
  inline bool update(uintE v, uintE w) { return updateAtomic(v, w); }
  inline bool updateAtomic(uintE v, uintE w) {
    int tv = S.belongs[v], tw = S.belongs[w];
    if (tv == tw || tv == -1 || tw == -1) return false;
    if (XBP_OPT) {
      // Claim the two trees in id order. A failed claim means another worker
      // holds that tree; following who holds what only climbs in tree id, so
      // whenever a pass has cross-tree edges some worker gets both trees and
      // records a path: no pass ends empty because of mutual back-off, and the
      // serial rescan below is not needed.
      const int lo = std::min(tv, tw), hi = std::max(tv, tw);
      if (!TryClaim(S.select_tree[lo])) return false;
      if (!TryClaim(S.select_tree[hi])) {
        S.select_tree[lo].store(0);
        return false;
      }
    } else {
      int expected = 0;
      if (!S.select_tree[tv].compare_exchange_strong(expected, 1)) {
        S.dropped[WorkerId()].push_back(v);
        return false;
      }
      if (!S.select_tree[tw].compare_exchange_strong(expected, 1)) {
        int taken = 1;
        S.select_tree[tv].compare_exchange_strong(taken, 0);
        S.dropped[WorkerId()].push_back(v);
        return false;
      }
    }
    if (XBP_TRIPLET) {
      std::vector<int> &ap = S.aug_pairs[WorkerId()];
      ap.push_back(v);
      ap.push_back(w);
      return false;
    }
    if (XBP_OPT) {
      // root(v) .. v, w .. root(w), written straight into this worker's buffer
      const int t = WorkerId();
      std::vector<int> &sc = S.scratch[t], &pn = S.path_nodes[t];
      sc.clear();
      AppendPathToRoot(S, v, sc);
      S.path_start[t].push_back(static_cast<int>(pn.size()));
      pn.insert(pn.end(), sc.rbegin(), sc.rend());
      AppendPathToRoot(S, w, pn);
      return false;
    }
    std::vector<int> pv = PathToRoot(S, v), pw = PathToRoot(S, w);
    std::vector<int> path(pv.rbegin(), pv.rend());
    path.insert(path.end(), pw.begin(), pw.end());
    std::lock_guard<std::mutex> g(S.paths_mutex);
    S.paths.push_back(std::move(path));
    return false;
  }
};

// 2. Expansion: w (in no tree) becomes odd, its mate x even.
struct Expand_F {
  static constexpr const char *kPhase = "expand";
  State &S;
  explicit Expand_F(State &s) : S(s) {}
  inline bool cond(uintE w) { return S.belongs[w] == -1; }
  inline bool update(uintE v, uintE w) { return updateAtomic(v, w); }
  inline bool updateAtomic(uintE v, uintE w) {
    int x = S.M[w];
    if (!TryClaim(S.select_match[std::min<int>(w, x)])) return false;
    if (XBP_TRIPLET) {
      SetRec(S, x, w, v, -1);
    } else {
      S.path_table[x].push_back(w);
      S.path_table[x].push_back(v);
    }
    S.is_even[w] = 0;
    S.is_even[x] = 1;
    S.belongs[w] = S.belongs[v];
    S.belongs[x] = S.belongs[v];
    S.emit(x);
    return false;
  }
};

// 3. Blossoms (XB-Pro's findBlossom / parBlossom).
static void FindBlossom(const State &S, int v, int w, std::vector<int> &blossom) {
  blossom.clear();
  int cv = v, cw = w, it = 0;
  while (cv != cw) {
    cv = S.path_table[cv].empty() ? w : S.path_table[cv].back();
    cw = S.path_table[cw].empty() ? v : S.path_table[cw].back();
    if (++it > static_cast<int>(S.n)) return;  // inconsistent read under concurrency
  }
  int base = cv;
  int cur = v;
  if (XBP_OPT) {
    // the v side front to back, reversed once: the same sequence as the front inserts below
    blossom.push_back(cur);
    while (cur != base && !S.path_table[cur].empty()) {
      blossom.insert(blossom.end(), S.path_table[cur].begin(), S.path_table[cur].end());
      cur = S.path_table[cur].back();
    }
    std::reverse(blossom.begin(), blossom.end());
    cur = w;
    blossom.push_back(cur);
    while (cur != base && !S.path_table[cur].empty()) {
      blossom.insert(blossom.end(), S.path_table[cur].begin(), S.path_table[cur].end());
      cur = blossom.back();
    }
    return;
  }
  blossom.insert(blossom.begin(), cur);
  while (cur != base && !S.path_table[cur].empty()) {
    blossom.insert(blossom.begin(), S.path_table[cur].rbegin(), S.path_table[cur].rend());
    cur = S.path_table[cur].back();
  }
  cur = w;
  blossom.push_back(cur);
  while (cur != base && !S.path_table[cur].empty()) {
    blossom.insert(blossom.end(), S.path_table[cur].begin(), S.path_table[cur].end());
    cur = blossom.back();
  }
}

static bool HasDuplicate(const std::vector<int> &p) {
  for (size_t c = 0; c < p.size(); c++)
    for (size_t d = c + 1; d < p.size(); d++)
      if (p[c] == p[d]) return true;
  return false;
}

struct Blossom_F {
  static constexpr const char *kPhase = "blossom";
  State &S;
  explicit Blossom_F(State &s) : S(s) {}
  inline bool cond(uintE w) { return S.is_even[w]; }
  inline bool update(uintE v, uintE w) { return updateAtomic(v, w); }

  inline void TryMakeEven(const std::vector<int> &blossom, int k, bool anticlockwise) {
    int current = blossom[k];
    if (S.blossom_to_base[current] == -1) S.blossom_to_base[current] = blossom[0];
    if (S.is_even[current] || !S.path_table[current].empty()) return;
    int expected = 0;
    if (!S.select_blossom[current].compare_exchange_strong(expected, 1)) return;
    // Build the path locally and publish it in one step: other workers may walk
    // this node's path as soon as it is even. A path with a repeated node comes
    // from an inconsistent read of the tree under concurrency; then release the
    // node so another blossom can still make it even in this phase (keeping the
    // claim would leave it odd for the rest of the phase and can lose the last
    // augmenting path).
    std::vector<int> pt;
    if (anticlockwise) {
      for (size_t m = k + 1; m < blossom.size(); m++) pt.push_back(blossom[m]);
    } else {
      for (int m = k - 1; m >= 0; m--) pt.push_back(blossom[m]);
    }
    if (HasDuplicate(pt)) {
      S.select_blossom[current].store(0);
      return;
    }
    S.path_table[current].swap(pt);
    std::atomic_thread_fence(std::memory_order_release);
    S.is_even[current] = 1;
    S.emit(current);
  }

  inline bool updateAtomic(uintE v, uintE w) {
    if (S.belongs[w] != S.belongs[v] || static_cast<int>(w) == S.M[v] || S.belongs[w] == -1) return false;
    if (XBP_TRIPLET) {
      long steps = 0;
      const int base = FindBaseByAnchor(S, v, w, &steps);
      if (base < 0) return false;
      steps += ClaimOddOnChain(S, v, w, base);
      steps += ClaimOddOnChain(S, w, v, base);
      gacge_walk(static_cast<unsigned long>(steps));  // one serial strand: base walk + both claims
      return false;
    }
    std::vector<int> local;
    std::vector<int> &blossom = XBP_OPT ? S.scratch[WorkerId()] : local;
    FindBlossom(S, v, w, blossom);
    if (blossom.size() < 3) return false;
    S.blossom_to_base[blossom[0]] = blossom[0];
    for (int k = static_cast<int>(blossom.size()) - 3; k >= 0; k -= 2) TryMakeEven(blossom, k, true);
    for (int k = 2; k < static_cast<int>(blossom.size()) - 1; k += 2) TryMakeEven(blossom, k, false);
    return false;
  }
};

// src[0, n) to dst, in parallel when it is large. An element loop: written as
// a cilk_for over 16 K blocks of std::copy it ran serially inside the search
// phase (Cilkscale on 48 workers: parallelism 1.0, a quarter of the matching's
// span on Wikipedia and LiveJournal; 2026-10-14), parallel only in isolation.
static void ParCopy(const uintE *src, size_t n, uintE *dst) {
  if (n <= (1 << 16)) {
    std::copy(src, src + n, dst);
    return;
  }
  parallel_for(size_t i = 0; i < n; i++) dst[i] = src[i];
}

// One edgeMap over `frontier` with functor f; returns nothing (results are in S).
template <class vertex, class F>
static void Step(graph<vertex> &GA, const Frontier &frontier, F f) {
  if (frontier.empty()) return;
  gacge_phase() = F::kPhase;
  // XBP_OPT: the vertex subset reads the frontier in place (edgeMap neither
  // changes nor frees a sparse subset's indices); XB-Pro copies it into a
  // fresh array every step.
  uintE *idx = XBP_OPT ? const_cast<uintE *>(frontier.data()) : newA(uintE, frontier.size());
  if (!XBP_OPT) std::copy(frontier.begin(), frontier.end(), idx);
  vertexSubset vs(GA.n, frontier.size(), idx);
  // Always sparse (push from the frontier): letting Ligra pull large steps
  // dense made no graph faster and Amazon, Patent, Hyperlink and LiveJournal
  // 1.2-2.7x slower (2026-10-13): augment and blossom keep their condition
  // after an update, so a pulled vertex scans all its edges, and the changed
  // claim order costs more phases on large graphs.
  edgeMap(GA, vs, f, -1, no_output | no_dense);
  if (!XBP_OPT) vs.del();
}

// dst = a followed by b (XBP_OPT), reusing dst's storage.
static void Concat(Frontier &dst, const Frontier &a, const Frontier &b) {
  dst.resize(a.size() + b.size());
  ParCopy(a.data(), a.size(), dst.data());
  ParCopy(b.data(), b.size(), dst.data() + a.size());
}

// The workers' emitted nodes into v (XBP_OPT): one parallel copy per worker.
static void TakeOutInto(State &S, Frontier &v) {
  const size_t W = S.out.size();
  std::vector<size_t> off(W + 1, 0);
  for (size_t t = 0; t < W; t++) off[t + 1] = off[t] + S.out[t].size();
  v.resize(off[W]);
  parallel_for(size_t t = 0; t < W; t++) {
    std::copy(S.out[t].begin(), S.out[t].end(), v.begin() + off[t]);
    S.out[t].clear();
  }
}

static Frontier TakeOut(State &S) {
  Frontier v;
  for (auto &o : S.out) {
    v.insert(v.end(), o.begin(), o.end());
    o.clear();
  }
  return v;
}

// One search phase; fills S.paths or leaves it empty when the forest is exhausted.
template <class vertex>
static void SearchPhase(graph<vertex> &GA, State &S) {
  const long n = S.n;
  // Reuse trees whose root is still exposed (XB-Pro's parReuseRemainingTrees).
  bool *roots = newA(bool, n), *kept_even = newA(bool, n);
  parallel_for(long x = 0; x < n; x++) {
    S.select_tree[x] = 0;
    S.select_match[x] = 0;
    S.select_blossom[x] = 0;
    S.blossom_to_base[x] = -1;
    roots[x] = kept_even[x] = false;
    int root = S.belongs[x];
    if (S.M[x] == -1) {
      // Case 1: exposed node, root of a new tree.
      S.is_even[x] = 1;
      S.belongs[x] = x;
      if (XBP_TRIPLET) SetRec(S, x, -1, -1, -1);
      else S.path_table[x].clear();
      roots[x] = true;
      if (XBP_COUNTERS) S.reset_tree_nodes.fetch_add(1, std::memory_order_relaxed);
    } else if (!XBP_REUSE || root == -1 || S.M[root] != -1) {
      // Case 2: in no tree, or in a tree whose root got matched (or no reuse).
      S.is_even[x] = 0;
      S.belongs[x] = -1;
      if (XBP_TRIPLET) SetRec(S, x, -1, -1, -1);
      else S.path_table[x].clear();
      if (XBP_COUNTERS) S.reset_tree_nodes.fetch_add(1, std::memory_order_relaxed);
    } else {
      // Case 3: in a tree whose root is still exposed: the tree is kept.
      if (S.is_even[x]) kept_even[x] = true;
      if (XBP_COUNTERS) S.reused_tree_nodes.fetch_add(1, std::memory_order_relaxed);
    }
  }
  _seq<uintE> r = sequence::packIndex<uintE>(roots, n);
  _seq<uintE> k = sequence::packIndex<uintE>(kept_even, n);
  Frontier v1, v2;
  if (XBP_OPT) {  // parallel copies (XB-Pro: serial ones, as std::vector's range constructor)
    v1.resize(r.n);
    v2.resize(k.n);
    ParCopy(r.A, r.n, v1.data());
    ParCopy(k.A, k.n, v2.data());
  } else {
    v1.assign(r.A, r.A + r.n);
    v2.assign(k.A, k.A + k.n);
  }
  r.del();
  k.del();
  free(roots);
  free(kept_even);

  Frontier frontier, blossom_frontier;
  if (XBP_OPT) {
    Concat(frontier, v1, v2);
  } else {
    frontier = v1;
    frontier.insert(frontier.end(), v2.begin(), v2.end());
  }
  gacge_cut("forest");  // GACGE_REGIONS: resetting the forest and packing its frontier
  while (true) {
    // 1. Augment. Two workers on the two ends of one edge can each claim one
    // tree, fail on the other and both back off, so neither records the path.
    // If the pass found nothing, every claim was released: rescanning the nodes
    // that gave up an edge on one worker (which cannot contend) settles it.
    // (XBP_OPT claims trees in id order instead, so a pass with cross-tree
    // edges always records a path and there is no rescan.)
    for (auto &d : S.dropped) d.clear();
    Step(GA, frontier, Augment_F(S));
    if (!XBP_OPT && !S.have_paths()) {
      std::vector<uintE> rescan;
      for (auto &d : S.dropped) {
        rescan.insert(rescan.end(), d.begin(), d.end());
        d.clear();
      }
      if (!rescan.empty()) {
        std::sort(rescan.begin(), rescan.end());
        rescan.erase(std::unique(rescan.begin(), rescan.end()), rescan.end());
        // One worker, written out rather than via setWorkers(1), which only the
        // OpenMP runtime can do at run time.
        Augment_F f(S);
        for (uintE v : rescan) {
          for (uintE j = 0; j < GA.V[v].getOutDegree(); j++) {
            uintE w = GA.V[v].getOutNeighbor(j);
            if (f.cond(w)) f.updateAtomic(v, w);
          }
        }
      }
    }
    if (S.have_paths()) return;

    // 2. Expand from every even node of the frontier: v1 = the new even nodes.
    Step(GA, frontier, Expand_F(S));
    if (XBP_OPT) TakeOutInto(S, v1);
    else v1 = TakeOut(S);
    // 3. Blossoms from the new even nodes and the previous blossom step's
    // (XB-Pro scans v1 + v2 here; v2 starts as the reused trees' even nodes).
    if (XBP_OPT) {
      Concat(blossom_frontier, v1, v2);
    } else {
      blossom_frontier = v1;
      blossom_frontier.insert(blossom_frontier.end(), v2.begin(), v2.end());
    }
    Step(GA, blossom_frontier, Blossom_F(S));
    if (XBP_TRIPLET) {
      // the claimed odd nodes become even around their blossoms
      parallel_for(size_t t = 0; t < S.tasks.size(); t++) {
        const std::vector<int> &tk = S.tasks[t];
        parallel_for(size_t i = 0; i < tk.size() / 4; i++) {
          const int u = tk[4 * i];
          SetRec(S, u, tk[4 * i + 1], tk[4 * i + 2], tk[4 * i + 3]);
          S.is_even[u] = 1;
          S.emit(u);
        }
      }
      for (auto &tk : S.tasks) tk.clear();
    }
    if (XBP_OPT) TakeOutInto(S, v2);
    else v2 = TakeOut(S);

    if (XBP_OPT) {
      Concat(frontier, v1, v2);
    } else {
      frontier = v1;
      frontier.insert(frontier.end(), v2.begin(), v2.end());
    }
    gacge_cut("apply");  // the blossom's claims applied, the next frontier packed
    if (frontier.empty()) return;
  }
}

template <class vertex>
static void MaximumMatching(graph<vertex> &GA, State &S) {
  while (true) {
    SearchPhase(GA, S);
    if (!S.have_paths()) return;
    if (XBP_TRIPLET) {
      // each augmenting edge (v, w) owns its two trees: the flips do not overlap
      parallel_for(size_t t = 0; t < S.aug_pairs.size(); t++) {
        const std::vector<int> &ap = S.aug_pairs[t];
        parallel_for(size_t i = 0; i < ap.size() / 2; i++) {
          const int v = ap[2 * i], w = ap[2 * i + 1];
          S.M[v] = w;
          S.M[w] = v;
          FlipRoute(S, v);
          FlipRoute(S, w);
        }
      }
      for (auto &ap : S.aug_pairs) ap.clear();
      gacge_cut("flip");  // GACGE_REGIONS: the phase's augmenting paths flipped
      continue;
    }
    if (XBP_OPT) {
      parallel_for(size_t t = 0; t < S.path_start.size(); t++) {
        const std::vector<int> &pn = S.path_nodes[t], &ps = S.path_start[t];
        parallel_for(size_t i = 0; i < ps.size(); i++) {
          const size_t end = i + 1 < ps.size() ? ps[i + 1] : pn.size();
          for (size_t j = ps[i]; j + 1 < end; j += 2) {
            S.M[pn[j]] = pn[j + 1];
            S.M[pn[j + 1]] = pn[j];
          }
        }
      }
      for (auto &p : S.path_nodes) p.clear();
      for (auto &p : S.path_start) p.clear();
      continue;
    }
    const long np = S.paths.size();
    parallel_for(long i = 0; i < np; i++) {
      const std::vector<int> &p = S.paths[i];
      for (size_t j = 0; j + 1 < p.size(); j += 2) {
        S.M[p[j]] = p[j + 1];
        S.M[p[j + 1]] = p[j];
      }
    }
    S.paths.clear();
  }
}

// Matching size if mate is a valid matching of GA (mutual, along edges), else -1.
template <class vertex>
static long CheckMatching(graph<vertex> &GA, const Ints &M) {
  const long n = GA.n;
  long *ok = newA(long, n);
  parallel_for(long v = 0; v < n; v++) {
    ok[v] = 0;
    int u = M[v];
    if (u == -1) continue;
    if (u < 0 || u >= n || M[u] != v) { ok[v] = -1; continue; }
    bool edge = false;
    for (uintE j = 0; j < GA.V[v].getOutDegree(); j++)
      if (static_cast<int>(GA.V[v].getOutNeighbor(j)) == u) { edge = true; break; }
    ok[v] = edge ? 1 : -1;
  }
  long matched = 0;
  bool valid = true;
  for (long v = 0; v < n; v++) {
    if (ok[v] < 0) valid = false;
    matched += ok[v] > 0;
  }
  free(ok);
  return valid ? matched / 2 : -1;
}

struct Run {
  long calls = 0, timed = 0, first_size = -2;
  double total = 0;
  bool all_valid = true, all_same = true;
};

}  // namespace xbp

template <class vertex>
void Compute(graph<vertex> &GA, commandLine P) {
  static xbp::Run R;
  const long warmup = P.getOptionLongValue("-warmup", 1);
  const long rounds = P.getOptionLongValue("-rounds", 3);
  if (R.calls == 0) {
    long max_degree = 0;
    for (long v = 0; v < GA.n; v++) max_degree = std::max<long>(max_degree, GA.V[v].getOutDegree());
    std::printf("XBConfig: arm=ligra variant=%s reuse=%d counters=%d opt=%d triplet=%d lb=ligra-edgemap dataset=%s nodes=%ld "
                "edges=%ld max_degree=%ld warmup=%ld rounds=%ld threads=%d load_first=%d\n",
                !XBP_OPT ? "xb-pro-ligra-xbpro" : !XBP_TRIPLET ? "xb-pro-ligra-vector" : XBP_REUSE ? "xb-pro-ligra" : "xb-ligra",
                XBP_REUSE, XBP_COUNTERS, XBP_OPT, XBP_TRIPLET, P.getOptionValue("-dataset", "unknown").c_str(), GA.n, GA.m / 2, max_degree, warmup,
                rounds - warmup, getWorkers(), kLoadFirst ? 1 : 0);
  }
  const long call = R.calls++;

  auto t0 = std::chrono::steady_clock::now();
  XBP_WS_BEGIN(ws_init);
  xbp::State S(GA.n);
  gacge_cut("init");  // GACGE_REGIONS: the state's allocation
  XBP_WS_END(ws_init, "init", call);
  XBP_WS_BEGIN(ws_match);
  xbp::MaximumMatching(GA, S);
  gacge_cut("tail");  // the last phase's end; the check follows as the run's last region
  XBP_WS_END(ws_match, "match", call);
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  if (call < warmup) return;
  XBP_WS_BEGIN(ws_check);
  long size = xbp::CheckMatching(GA, S.M);
  XBP_WS_END(ws_check, "check", call);
  bool valid = size >= 0;
  R.all_valid = R.all_valid && valid;
  if (R.first_size == -2) R.first_size = size;
  else if (size != R.first_size) R.all_same = false;
  R.total += secs;
  const long index = R.timed++;
  std::printf("XBRound: index=%ld runtime_s=%.9f matching_size=%ld valid=%d\n", index, secs, size,
              valid ? 1 : 0);
  if (XBP_COUNTERS)
    std::printf("XBCounters: index=%ld reused_tree_nodes=%llu reset_tree_nodes=%llu\n", index,
                S.reused_tree_nodes.load(), S.reset_tree_nodes.load());
  if (call == rounds - 1) {
    const char *status = !R.all_valid ? "invalid_matching" : (!R.all_same ? "matching_size_varies" : "ok");
    std::printf("XBResult: status=%s matching_size=%ld valid=%d rounds=%ld mean_runtime_s=%.9f\n", status,
                R.first_size, R.all_valid ? 1 : 0, R.timed, R.timed ? R.total / R.timed : 0.0);
  }
  std::fflush(stdout);
}
