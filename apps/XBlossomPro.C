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
// XBP_SPAN 1 = (on top of XBP_OPT; XBlossomPro_span_cilk, not the default) also
//              shorten what is serial inside a step: a blossom is checked for a
//              repeated node once (sorted copy) instead of every node's new path
//              being checked pair by pair (O(L^3) per blossom edge), the phase's
//              initial frontiers and each worker's emitted nodes are copied in
//              parallel, and frontier vectors skip the serial zero fill on resize
#ifndef XBP_SPAN
#define XBP_SPAN 0
#endif
#if XBP_SPAN && !XBP_OPT
#error "XBP_SPAN builds on XBP_OPT"
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
#include <mutex>
#include <string>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

// Cilkscale (GACGE_CILKSCALE): the work and span of a part of one Compute(),
// one CSV row tagged <name><call> (init, match, check) next to ligra.h's run<r>.
#ifdef GACGE_CILKSCALE
#define XBP_WS_BEGIN(v) wsp_t v = wsp_getworkspan()
#define XBP_WS_END(v, name, call)                                     \
  do {                                                                \
    char t_[32];                                                      \
    std::snprintf(t_, sizeof t_, "%s%ld", name, static_cast<long>(call)); \
    wsp_dump(wsp_sub(wsp_getworkspan(), v), t_);                      \
  } while (0)
#else
#define XBP_WS_BEGIN(v) (void)0
#define XBP_WS_END(v, name, call) (void)0
#endif

namespace xbp {

struct State {
  long n = 0;
  std::vector<int> M;         // mate, -1 = exposed
  std::vector<int> is_even;   // 1 = even node of a tree
  std::vector<int> belongs;   // tree (root) of a node, -1 = none
  std::vector<int> blossom_to_base;
  std::vector<std::vector<int>> path_table;  // node -> nodes towards its root
  std::vector<std::atomic<int>> select_tree, select_match, select_blossom;

  std::vector<std::vector<int>> paths;  // augmenting paths of this phase
  std::mutex paths_mutex;
  // XBP_OPT: the augmenting paths per worker, path i of worker t being
  // path_nodes[t][path_start[t][i] .. path_start[t][i + 1] or the end); and a
  // per-worker scratch buffer for one root path or one blossom.
  std::vector<std::vector<int>> path_nodes, path_start, scratch;
  std::vector<std::vector<int>> scratch2;  // XBP_SPAN: a blossom's sorted copy
  // Frontier nodes whose scan gave up a cross-tree edge under contention, per worker.
  std::vector<std::vector<uintE>> dropped;

  // Nodes produced by expand / blossom, one buffer per worker (as XB-Pro's
  // thread-local vectors: one shared counter would make every emit contend).
  std::vector<std::vector<uintE>> out;

  // Tree nodes reused / reset at the start of each search phase (XBP_COUNTERS).
  std::atomic<unsigned long long> reused_tree_nodes{0}, reset_tree_nodes{0};

  explicit State(long n_) : n(n_), M(n_, -1), is_even(n_, 0), belongs(n_, -1), blossom_to_base(n_, -1),
                            path_table(n_), select_tree(n_), select_match(n_), select_blossom(n_),
                            dropped(std::max(1, getWorkers())), out(std::max(1, getWorkers())),
                            path_nodes(XBP_OPT ? std::max(1, getWorkers()) : 0),
                            path_start(XBP_OPT ? std::max(1, getWorkers()) : 0),
                            scratch(XBP_OPT ? std::max(1, getWorkers()) : 0),
                            scratch2(XBP_SPAN ? std::max(1, getWorkers()) : 0) {
    // as XB-Pro: limits reallocation under concurrent readers
    if (XBP_OPT) {
      parallel_for(long i = 0; i < n_; i++) path_table[i].reserve(100);
    } else {
      for (auto &p : path_table) p.reserve(100);
    }
  }

  bool have_paths() const {
    if (!XBP_OPT) return !paths.empty();
    for (auto &p : path_start)
      if (!p.empty()) return true;
    return false;
  }

  void emit(int v) { out[WorkerId()].push_back(v); }
};

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
      int e1 = 0, e2 = 0;
      if (!S.select_tree[lo].compare_exchange_strong(e1, 1)) return false;
      if (!S.select_tree[hi].compare_exchange_strong(e2, 1)) {
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
  State &S;
  explicit Expand_F(State &s) : S(s) {}
  inline bool cond(uintE w) { return S.belongs[w] == -1; }
  inline bool update(uintE v, uintE w) { return updateAtomic(v, w); }
  inline bool updateAtomic(uintE v, uintE w) {
    int x = S.M[w];
    int expected = 0;
    if (!S.select_match[std::min<int>(w, x)].compare_exchange_strong(expected, 1)) return false;
    S.path_table[x].push_back(w);
    S.path_table[x].push_back(v);
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

// XBP_SPAN: whether p repeats a node, from a sorted copy in tmp (O(L log L)).
static bool HasDuplicateSorted(const std::vector<int> &p, std::vector<int> &tmp) {
  tmp.assign(p.begin(), p.end());
  std::sort(tmp.begin(), tmp.end());
  return std::adjacent_find(tmp.begin(), tmp.end()) != tmp.end();
}

static bool HasDuplicate(const std::vector<int> &p) {
  for (size_t c = 0; c < p.size(); c++)
    for (size_t d = c + 1; d < p.size(); d++)
      if (p[c] == p[d]) return true;
  return false;
}

struct Blossom_F {
  State &S;
  explicit Blossom_F(State &s) : S(s) {}
  inline bool cond(uintE w) { return S.is_even[w]; }
  inline bool update(uintE v, uintE w) { return updateAtomic(v, w); }

  // check: test the node's new path for a repeated node (always, except under
  // XBP_SPAN when the whole blossom has none: the path is a part of it).
  inline void TryMakeEven(const std::vector<int> &blossom, int k, bool anticlockwise, bool check) {
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
    if (check && HasDuplicate(pt)) {
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
    std::vector<int> local;
    std::vector<int> &blossom = XBP_OPT ? S.scratch[WorkerId()] : local;
    FindBlossom(S, v, w, blossom);
    if (blossom.size() < 3) return false;
    const bool check = !XBP_SPAN || HasDuplicateSorted(blossom, S.scratch2[WorkerId()]);
    S.blossom_to_base[blossom[0]] = blossom[0];
    for (int k = static_cast<int>(blossom.size()) - 3; k >= 0; k -= 2) TryMakeEven(blossom, k, true, check);
    for (int k = 2; k < static_cast<int>(blossom.size()) - 1; k += 2) TryMakeEven(blossom, k, false, check);
    return false;
  }
};

// XBP_SPAN: an allocator whose resize leaves new elements uninitialized (the
// parallel copy that follows fills them), so frontier vectors skip a serial
// zero fill per step.
template <class T>
struct NoInitAlloc : std::allocator<T> {
  template <class U>
  struct rebind { using other = NoInitAlloc<U>; };
  NoInitAlloc() = default;
  template <class U>
  NoInitAlloc(const NoInitAlloc<U> &) {}
  template <class U>
  void construct(U *p) noexcept { ::new (static_cast<void *>(p)) U; }
  template <class U, class... A>
  void construct(U *p, A &&...a) { ::new (static_cast<void *>(p)) U(std::forward<A>(a)...); }
};
using Frontier = std::conditional<XBP_SPAN, std::vector<uintE, NoInitAlloc<uintE>>, std::vector<uintE>>::type;

// src[0, n) to dst, in parallel blocks when it is large.
static void ParCopy(const uintE *src, size_t n, uintE *dst) {
  const size_t B = 1 << 14;
  if (n <= 4 * B) {
    std::copy(src, src + n, dst);
    return;
  }
  parallel_for(size_t b = 0; b < (n + B - 1) / B; b++) {
    const size_t lo = b * B, hi = std::min(n, lo + B);
    std::copy(src + lo, src + hi, dst + lo);
  }
}

// One edgeMap over `frontier` with functor f; returns nothing (results are in S).
template <class vertex, class F>
static void Step(graph<vertex> &GA, const Frontier &frontier, F f) {
  if (frontier.empty()) return;
  uintE *idx = newA(uintE, frontier.size());
  if (XBP_OPT) ParCopy(frontier.data(), frontier.size(), idx);
  else std::copy(frontier.begin(), frontier.end(), idx);
  vertexSubset vs(GA.n, frontier.size(), idx);
  edgeMap(GA, vs, f, -1, no_output | no_dense);
  vs.del();
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
    if (XBP_SPAN) ParCopy(S.out[t].data(), S.out[t].size(), v.data() + off[t]);
    else std::copy(S.out[t].begin(), S.out[t].end(), v.begin() + off[t]);
    S.out[t].clear();
  }
}

// v = TakeOut's result: a move when v is a plain vector, a copy otherwise.
static void Assign(std::vector<uintE> &d, std::vector<uintE> &&s) { d = std::move(s); }
template <class A>
static void Assign(std::vector<uintE, A> &d, std::vector<uintE> &&s) { d.assign(s.begin(), s.end()); }

static std::vector<uintE> TakeOut(State &S) {
  std::vector<uintE> v;
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
      S.path_table[x].clear();
      roots[x] = true;
      if (XBP_COUNTERS) S.reset_tree_nodes.fetch_add(1, std::memory_order_relaxed);
    } else if (!XBP_REUSE || root == -1 || S.M[root] != -1) {
      // Case 2: in no tree, or in a tree whose root got matched (or no reuse).
      S.is_even[x] = 0;
      S.belongs[x] = -1;
      S.path_table[x].clear();
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
  if (XBP_SPAN) {
    v1.resize(r.n);
    ParCopy(r.A, r.n, v1.data());
    v2.resize(k.n);
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
    else Assign(v1, TakeOut(S));
    // 3. Blossoms from the new even nodes and the previous blossom step's
    // (XB-Pro scans v1 + v2 here; v2 starts as the reused trees' even nodes).
    if (XBP_OPT) {
      Concat(blossom_frontier, v1, v2);
    } else {
      blossom_frontier = v1;
      blossom_frontier.insert(blossom_frontier.end(), v2.begin(), v2.end());
    }
    Step(GA, blossom_frontier, Blossom_F(S));
    if (XBP_OPT) TakeOutInto(S, v2);
    else Assign(v2, TakeOut(S));

    if (XBP_OPT) {
      Concat(frontier, v1, v2);
    } else {
      frontier = v1;
      frontier.insert(frontier.end(), v2.begin(), v2.end());
    }
    if (frontier.empty()) return;
  }
}

template <class vertex>
static void MaximumMatching(graph<vertex> &GA, State &S) {
  while (true) {
    SearchPhase(GA, S);
    if (!S.have_paths()) return;
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
static long CheckMatching(graph<vertex> &GA, const std::vector<int> &M) {
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
    std::printf("XBConfig: arm=ligra variant=%s reuse=%d counters=%d opt=%d lb=ligra-edgemap dataset=%s nodes=%ld "
                "edges=%ld max_degree=%ld warmup=%ld rounds=%ld threads=%d\n",
                !XBP_OPT ? "xb-pro-ligra-xbpro" : XBP_SPAN ? "xb-pro-ligra-span" : XBP_REUSE ? "xb-pro-ligra" : "xb-ligra",
                XBP_REUSE, XBP_COUNTERS, XBP_OPT, P.getOptionValue("-dataset", "unknown").c_str(), GA.n, GA.m / 2, max_degree, warmup,
                rounds - warmup, getWorkers());
  }
  const long call = R.calls++;

  auto t0 = std::chrono::steady_clock::now();
  XBP_WS_BEGIN(ws_init);
  xbp::State S(GA.n);
  XBP_WS_END(ws_init, "init", call);
  XBP_WS_BEGIN(ws_match);
  xbp::MaximumMatching(GA, S);
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
