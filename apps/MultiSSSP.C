// Multi-source SSSP for standard weighted adjacency graphs.
// Sources are read from -srcs <file>, one vertex id per line; lines beginning
// with '#' are ignored. If -srcs is omitted, -r <source> is used.
#define WEIGHTED 1
#include "ligra.h"
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

using std::string;

typedef float weightT;

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
    return writeMin(&Dist[d], newDist) && CAS(&Visited[d], 0, 1);
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

static void write_distances(const char* output_path, weightT* dist, long n) {
  std::ofstream out(output_path);
  if (!out) {
    std::cerr << "Could not open MultiSSSP output file: " << output_path << std::endl;
    abort();
  }
  for (long i = 0; i < n; i++) {
    out << i << " " << dist[i] << "\n";
  }
}

static bool parse_source_line(const string& line, long* source) {
  string trimmed = line;
  size_t first = trimmed.find_first_not_of(" \t\r\n");
  if (first == string::npos || trimmed[first] == '#') return false;

  std::istringstream in(trimmed.substr(first));
  return static_cast<bool>(in >> *source);
}

static vertexSubset make_initial_frontier(commandLine& P, weightT* Dist, long n) {
  bool* sources = newA(bool, n);
  parallel_for(long i = 0; i < n; i++) sources[i] = false;
  long source_count = 0;

  char* sources_path = P.getOptionValue("-srcs");
  if (sources_path != NULL) {
    std::ifstream in(sources_path);
    if (!in) {
      std::cerr << "Could not open MultiSSSP source file: " << sources_path << std::endl;
      abort();
    }

    string line;
    while (std::getline(in, line)) {
      long source = 0;
      if (!parse_source_line(line, &source)) continue;
      if (source < 0 || source >= n) {
        std::cerr << "Source vertex out of range in " << sources_path << ": "
                  << source << std::endl;
        abort();
      }
      if (!sources[source]) {
        sources[source] = true;
        Dist[source] = 0;
        source_count++;
      }
    }
  } else {
    long source = P.getOptionLongValue("-r", 0);
    if (source < 0 || source >= n) {
      std::cerr << "Source vertex out of range: " << source << std::endl;
      abort();
    }
    sources[source] = true;
    Dist[source] = 0;
    source_count = 1;
  }

  if (source_count == 0) {
    std::cerr << "MultiSSSP requires at least one source vertex" << std::endl;
    abort();
  }

  return vertexSubset(n, source_count, sources);
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
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
    vertexSubset output =
        edgeMap(GA, Frontier, SSSP_F(Dist, Visited), GA.m / 20, dense_forward);
    vertexMap(output, Reset_F(Visited));
    Frontier.del();
    Frontier = output;
    round++;
  }

  char* output_path = P.getOptionValue("-o");
  if (output_path != NULL) {
    write_distances(output_path, Dist, n);
  }

  Frontier.del();
  free(Visited);
  free(Dist);
}
