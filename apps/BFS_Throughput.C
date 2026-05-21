// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra.h"
#include <chrono>

struct BFS_F {
  uintE* Parents;
  BFS_F(uintE* _Parents) : Parents(_Parents) {}
  inline bool update (uintE s, uintE d) { //Update
    if(Parents[d] == UINT_E_MAX) { Parents[d] = s; return 1; }
    else return 0;
  }
  inline bool updateAtomic (uintE s, uintE d){ //atomic version of Update
    return (CAS(&Parents[d],UINT_E_MAX,s));
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return (Parents[d] == UINT_E_MAX); } 
};

template <class vertex>
long CountFrontierEdges(graph<vertex>& GA, vertexSubset& Frontier) {
  long frontier_size = Frontier.numNonzeros();
  if (frontier_size == 0) return 0;

  if (Frontier.dense()) {
    long* Degrees = newA(long, GA.n);
    parallel_for(long i=0;i<GA.n;i++) {
      Degrees[i] = Frontier.isIn(i) ? GA.V[i].getOutDegree() : 0;
    }
    long total_degree = sequence::plusReduce(Degrees, GA.n);
    free(Degrees);
    return total_degree;
  }

  long* Degrees = newA(long, frontier_size);
  parallel_for(long i=0;i<frontier_size;i++) {
    Degrees[i] = GA.V[Frontier.vtx(i)].getOutDegree();
  }
  long total_degree = sequence::plusReduce(Degrees, frontier_size);
  free(Degrees);
  return total_degree;
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  auto throughput_start = std::chrono::high_resolution_clock::now();
  long edges_processed = 0;
  long start = P.getOptionLongValue("-r",0);
  long n = GA.n;
  //creates Parents array, initialized to all -1, except for start
  uintE* Parents = newA(uintE,n);
  parallel_for(long i=0;i<n;i++) Parents[i] = UINT_E_MAX;
  Parents[start] = start;
  vertexSubset Frontier(n,start); //creates initial frontier
  while(!Frontier.isEmpty()){ //loop until frontier is empty
    edges_processed += CountFrontierEdges(GA, Frontier);
    vertexSubset output = edgeMap(GA, Frontier, BFS_F(Parents));
    Frontier.del();
    Frontier = output; //set new frontier
  } 
  Frontier.del();
  free(Parents); 

  auto throughput_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime = throughput_end - throughput_start;
  std::cout << "Average runtime: " << runtime.count() << std::endl;
  std::cout << "Number of Processed Edges: " << edges_processed << std::endl;
  std::cout << "Processed Edges per Second: "
            << static_cast<double>(edges_processed) / runtime.count()
            << std::endl;
}
