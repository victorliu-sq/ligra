// This code is adapted from Ligra's BellmanFord and HyperSSSP examples.
// It runs single-source shortest path on a standard weighted adjacency graph.
#define WEIGHTED 1
#include "ligra.h"
#include <fstream>
#include <limits>

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
    std::cerr << "Could not open SSSP output file: " << output_path << std::endl;
    abort();
  }
  for (long i = 0; i < n; i++) {
    out << i << " " << dist[i] << "\n";
  }
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long start = P.getOptionLongValue("-r", 0);
  long n = GA.n;
  if (start < 0 || start >= n) {
    std::cerr << "Source vertex out of range: " << start << std::endl;
    abort();
  }

  weightT* Dist = newA(weightT, n);
  parallel_for(long i = 0; i < n; i++) Dist[i] = std::numeric_limits<weightT>::max() / 4;
  Dist[start] = 0;

  int* Visited = newA(int, n);
  parallel_for(long i = 0; i < n; i++) Visited[i] = 0;

  vertexSubset Frontier(n, start);
  long round = 0;
  while (!Frontier.isEmpty()) {
    if (round == n) {
      parallel_for(long i = 0; i < n; i++) Dist[i] = -(std::numeric_limits<weightT>::max() / 4);
      break;
    }
    vertexSubset output = edgeMap(GA, Frontier, SSSP_F(Dist, Visited));
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
