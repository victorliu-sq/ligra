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
//   ./XBlossomPro -s -rounds 11 -warmup 1 -dataset Amazon <adj file>
#include "ligra.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <mutex>
#include <string>
#include <vector>

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
  std::atomic<bool> contended{false};

  std::vector<uintE> out;  // nodes produced by expand / blossom
  std::atomic<long> out_n{0};

  explicit State(long n_) : n(n_), M(n_, -1), is_even(n_, 0), belongs(n_, -1), blossom_to_base(n_, -1),
                            path_table(n_), select_tree(n_), select_match(n_), select_blossom(n_), out(n_) {
    for (auto &p : path_table) p.reserve(100);  // as XB-Pro: limits reallocation under concurrent readers
  }

  void emit(int v) { out[out_n.fetch_add(1, std::memory_order_relaxed)] = v; }
};

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
    int expected = 0;
    if (!S.select_tree[tv].compare_exchange_strong(expected, 1)) {
      S.contended.store(true, std::memory_order_relaxed);
      return false;
    }
    if (!S.select_tree[tw].compare_exchange_strong(expected, 1)) {
      int taken = 1;
      S.select_tree[tv].compare_exchange_strong(taken, 0);
      S.contended.store(true, std::memory_order_relaxed);
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
    std::vector<int> &pt = S.path_table[current];
    if (anticlockwise) {
      for (size_t m = k + 1; m < blossom.size(); m++) pt.push_back(blossom[m]);
    } else {
      for (int m = k - 1; m >= 0; m--) pt.push_back(blossom[m]);
    }
    if (HasDuplicate(pt)) {
      pt.clear();
      return;
    }
    S.emit(current);
    S.is_even[current] = 1;
  }

  inline bool updateAtomic(uintE v, uintE w) {
    if (S.belongs[w] != S.belongs[v] || static_cast<int>(w) == S.M[v] || S.belongs[w] == -1) return false;
    std::vector<int> blossom;
    FindBlossom(S, v, w, blossom);
    if (blossom.size() < 3) return false;
    S.blossom_to_base[blossom[0]] = blossom[0];
    for (int k = static_cast<int>(blossom.size()) - 3; k >= 0; k -= 2) TryMakeEven(blossom, k, true);
    for (int k = 2; k < static_cast<int>(blossom.size()) - 1; k += 2) TryMakeEven(blossom, k, false);
    return false;
  }
};

// One edgeMap over `frontier` with functor f; returns nothing (results are in S).
template <class vertex, class F>
static void Step(graph<vertex> &GA, const std::vector<uintE> &frontier, F f) {
  if (frontier.empty()) return;
  uintE *idx = newA(uintE, frontier.size());
  std::copy(frontier.begin(), frontier.end(), idx);
  vertexSubset vs(GA.n, frontier.size(), idx);
  edgeMap(GA, vs, f, -1, no_output | no_dense);
  vs.del();
}

static std::vector<uintE> TakeOut(State &S) {
  std::vector<uintE> v(S.out.begin(), S.out.begin() + S.out_n.load());
  S.out_n = 0;
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
      S.is_even[x] = 1;
      S.belongs[x] = x;
      S.path_table[x].clear();
      roots[x] = true;
    } else if (root == -1 || S.M[root] != -1) {
      S.is_even[x] = 0;
      S.belongs[x] = -1;
      S.path_table[x].clear();
    } else if (S.is_even[x]) {
      kept_even[x] = true;
    }
  }
  _seq<uintE> r = sequence::packIndex<uintE>(roots, n);
  _seq<uintE> k = sequence::packIndex<uintE>(kept_even, n);
  std::vector<uintE> v1(r.A, r.A + r.n), v2(k.A, k.A + k.n);
  r.del();
  k.del();
  free(roots);
  free(kept_even);

  std::vector<uintE> frontier(v1);
  frontier.insert(frontier.end(), v2.begin(), v2.end());
  while (true) {
    // 1. Augment. A pass that found nothing but gave up contended edges is not
    // conclusive (two workers on the two ends of one edge can both back off);
    // rerun it, the third time on one worker, which cannot contend.
    S.contended = false;
    Step(GA, frontier, Augment_F(S));
    for (int retry = 0; S.paths.empty() && S.contended.exchange(false); retry++) {
      int workers = getWorkers();
      if (retry >= 2) setWorkers(1);
      Step(GA, frontier, Augment_F(S));
      if (retry >= 2) setWorkers(workers);
    }
    if (!S.paths.empty()) return;

    // 2. Expand from every even node of the frontier: v1 = the new even nodes.
    Step(GA, frontier, Expand_F(S));
    v1 = TakeOut(S);
    // 3. Blossoms from the new even nodes and the previous blossom step's
    // (XB-Pro scans v1 + v2 here; v2 starts as the reused trees' even nodes).
    std::vector<uintE> blossom_frontier(v1);
    blossom_frontier.insert(blossom_frontier.end(), v2.begin(), v2.end());
    Step(GA, blossom_frontier, Blossom_F(S));
    v2 = TakeOut(S);

    frontier = v1;
    frontier.insert(frontier.end(), v2.begin(), v2.end());
    if (frontier.empty()) return;
  }
}

template <class vertex>
static void MaximumMatching(graph<vertex> &GA, State &S) {
  while (true) {
    SearchPhase(GA, S);
    if (S.paths.empty()) return;
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
    std::printf("XBConfig: arm=ligra variant=xb-pro-ligra reuse=1 lb=ligra-edgemap dataset=%s nodes=%ld "
                "edges=%ld max_degree=%ld warmup=%ld rounds=%ld threads=%d\n",
                P.getOptionValue("-dataset", "unknown").c_str(), GA.n, GA.m / 2, max_degree, warmup,
                rounds - warmup, getWorkers());
  }
  const long call = R.calls++;

  auto t0 = std::chrono::steady_clock::now();
  xbp::State S(GA.n);
  xbp::MaximumMatching(GA, S);
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  if (call < warmup) return;
  long size = xbp::CheckMatching(GA, S.M);
  bool valid = size >= 0;
  R.all_valid = R.all_valid && valid;
  if (R.first_size == -2) R.first_size = size;
  else if (size != R.first_size) R.all_same = false;
  R.total += secs;
  std::printf("XBRound: index=%ld runtime_s=%.9f matching_size=%ld valid=%d\n", R.timed++, secs, size,
              valid ? 1 : 0);
  if (call == rounds - 1) {
    const char *status = !R.all_valid ? "invalid_matching" : (!R.all_same ? "matching_size_varies" : "ok");
    std::printf("XBResult: status=%s matching_size=%ld valid=%d rounds=%ld mean_runtime_s=%.9f\n", status,
                R.first_size, R.all_valid ? 1 : 0, R.timed, R.timed ? R.total / R.timed : 0.0);
  }
  std::fflush(stdout);
}
