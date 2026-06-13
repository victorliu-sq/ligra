// Multi-source SSSP throughput variant for standard weighted adjacency graphs.
// Sources are provided either by --src 0 1 3 or randomly by --src-count 10.
#define WEIGHTED 1
#include "ligra.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>

using std::string;

typedef float weightT;
static double source_setup_seconds = 0.0;

struct SSSP_F {
  weightT* Dist;
  int* Visited;
  SSSP_F(weightT* _Dist, int* _Visited) : Dist(_Dist), Visited(_Visited) {}

  inline bool update(uintE s, uintE d, intE edgeLen) {
    weightT newDist = Dist[s] + static_cast<weightT>(edgeLen);
    if (Dist[d] > newDist) {
      Dist[d] = newDist;
      if (Visited[d] == 0) {
        Visited[d] = 1;
        return 1;
      }
    }
    return 0;
  }

  inline bool updateAtomic(uintE s, uintE d, intE edgeLen) {
    weightT newDist = Dist[s] + static_cast<weightT>(edgeLen);
    bool improved = writeMin(&Dist[d], newDist);
    return improved && CAS(&Visited[d], 0, 1);
  }

  inline bool cond(uintE d) { return cond_true(d); }
};

struct Reset_F {
  int* Visited;
  Reset_F(int* _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 0;
    return 1;
  }
};

static void add_source(long source, const string& context, bool* sources,
                       weightT* Dist, long n, long* source_count) {
  if (source < 0 || source >= n) {
    std::cerr << "Source vertex out of range in " << context << ": "
              << source << std::endl;
    abort();
  }
  if (!sources[source]) {
    sources[source] = true;
    Dist[source] = 0;
    (*source_count)++;
  }
}

static void parse_source_text(string text, const string& context, bool* sources,
                              weightT* Dist, long n, long* source_count) {
  std::replace(text.begin(), text.end(), ',', ' ');
  std::istringstream in(text);
  string token;
  while (in >> token) {
    char* end = NULL;
    long source = std::strtol(token.c_str(), &end, 10);
    if (end == token.c_str() || *end != '\0') {
      std::cerr << "Invalid source vertex in " << context << ": "
                << token << std::endl;
      abort();
    }
    add_source(source, context, sources, Dist, n, source_count);
  }
}

static void add_random_sources(long requested_count, bool* sources,
                               weightT* Dist, long n, long* source_count) {
  if (requested_count <= 0 || requested_count > n) {
    std::cerr << "--src-count must be between 1 and the number of graph vertices ("
              << n << "): " << requested_count << std::endl;
    abort();
  }

  std::unordered_set<long> selected;
  selected.reserve(requested_count);
  std::random_device seed;
  std::mt19937_64 engine(seed());
  std::uniform_int_distribution<long> dist(0, n - 1);

  while (static_cast<long>(selected.size()) < requested_count) {
    selected.insert(dist(engine));
  }

  for (long source : selected) {
    add_source(source, "--src-count", sources, Dist, n, source_count);
  }
}

static void add_random_sources_untimed(long requested_count, bool* sources,
                                       weightT* Dist, long n,
                                       long* source_count) {
  bool timer_was_on = _tm.on;
  if (timer_was_on) _tm.stop();
  auto setup_start = std::chrono::high_resolution_clock::now();
  add_random_sources(requested_count, sources, Dist, n, source_count);
  auto setup_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> setup_time = setup_end - setup_start;
  source_setup_seconds += setup_time.count();
  if (timer_was_on) _tm.start();
}

static vertexSubset make_initial_frontier(commandLine& P, weightT* Dist, long n) {
  bool* sources = newA(bool, n);
  parallel_for(long i = 0; i < n; i++) sources[i] = false;
  long source_count = 0;

  char* source_count_spec = P.getOptionValue("--src-count");
  bool explicit_sources = false;
  for (int i = 1; i < P.argc - 1; i++) {
    string arg = P.argv[i];
    if (arg == "--src" || arg == "-srcs") {
      explicit_sources = true;
      i++;
      while (i < P.argc - 1 && P.argv[i][0] != '-') {
        parse_source_text(P.argv[i], "--src", sources, Dist, n, &source_count);
        i++;
      }
      i--;
    }
  }

  if (explicit_sources && source_count_spec != NULL) {
    std::cerr << "Use either --src or --src-count, not both" << std::endl;
    abort();
  }

  if (source_count_spec != NULL) {
    char* end = NULL;
    long requested_count = std::strtol(source_count_spec, &end, 10);
    if (end == source_count_spec || *end != '\0') {
      std::cerr << "Invalid --src-count value: " << source_count_spec << std::endl;
      abort();
    }
    add_random_sources_untimed(requested_count, sources, Dist, n, &source_count);
  } else if (!explicit_sources) {
    long source = P.getOptionLongValue("-r", 0);
    add_source(source, "-r", sources, Dist, n, &source_count);
  }

  if (source_count == 0) {
    std::cerr << "MultiSSSP requires at least one source vertex" << std::endl;
    abort();
  }

  return vertexSubset(n, source_count, sources);
}

template <class vertex>
long CountFrontierEdges(graph<vertex>& GA, vertexSubset& Frontier) {
  long frontier_size = Frontier.numNonzeros();
  if (frontier_size == 0) return 0;

  if (Frontier.dense()) {
    long* Degrees = newA(long, GA.n);
    parallel_for(long i = 0; i < GA.n; i++) {
      Degrees[i] = Frontier.isIn(i) ? GA.V[i].getOutDegree() : 0;
    }
    long total_degree = sequence::plusReduce(Degrees, GA.n);
    free(Degrees);
    return total_degree;
  }

  long* Degrees = newA(long, frontier_size);
  parallel_for(long i = 0; i < frontier_size; i++) {
    Degrees[i] = GA.V[Frontier.vtx(i)].getOutDegree();
  }
  long total_degree = sequence::plusReduce(Degrees, frontier_size);
  free(Degrees);
  return total_degree;
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  source_setup_seconds = 0.0;
  auto throughput_start = std::chrono::high_resolution_clock::now();
  long edges_processed = 0;
  long n = GA.n;

  weightT* Dist = newA(weightT, n);
  weightT infinity = std::numeric_limits<weightT>::max() / 4;
  parallel_for(long i = 0; i < n; i++) Dist[i] = infinity;

  vertexSubset Frontier = make_initial_frontier(P, Dist, n);

  int* Visited = newA(int, n);
  parallel_for(long i = 0; i < n; i++) Visited[i] = 0;

  long round = 0;
  while (!Frontier.isEmpty()) {
    if (round == n) {
      parallel_for(long i = 0; i < n; i++) Dist[i] = -infinity;
      break;
    }
    edges_processed += CountFrontierEdges(GA, Frontier);
    vertexSubset output = edgeMap(GA, Frontier, SSSP_F(Dist, Visited));
    vertexMap(output, Reset_F(Visited));
    Frontier.del();
    Frontier = output;
    round++;
  }

  Frontier.del();
  free(Visited);
  free(Dist);

  auto throughput_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime = throughput_end - throughput_start;
  double algorithm_runtime = runtime.count() - source_setup_seconds;
  if (algorithm_runtime <= 0.0) algorithm_runtime = runtime.count();
  std::cout << "Average runtime: " << algorithm_runtime << std::endl;
  std::cout << "Number of Processed Edges: " << edges_processed << std::endl;
  std::cout << "Processed Edges per Second: "
            << static_cast<double>(edges_processed) / algorithm_runtime
            << std::endl;
}
