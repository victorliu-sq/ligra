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
#define HYPER 1
#define WEIGHTED 1
#include "hygra.h"
#include <chrono>

struct BF_Relax_F {
  intE* ShortestPathLenSrc, *ShortestPathLenDest;
  int* Visited;
  BF_Relax_F(intE* _ShortestPathLenSrc, intE* _ShortestPathLenDest, int* _Visited) :
    ShortestPathLenSrc(_ShortestPathLenSrc), ShortestPathLenDest(_ShortestPathLenDest), Visited(_Visited) {}
  inline bool update (uintE s, uintE d, intE edgeLen) { //Update ShortestPathLen if found a shorter path
    intE newDist = ShortestPathLenSrc[s] + edgeLen;
    if(ShortestPathLenDest[d] > newDist) {
      ShortestPathLenDest[d] = newDist;
      if(Visited[d] == 0) { Visited[d] = 1 ; return 1;}
    }
    return 0;
  }
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){ //atomic Update
    intE newDist = ShortestPathLenSrc[s] + edgeLen;
    bool improved = writeMin(&ShortestPathLenDest[d],newDist);
    return (improved && CAS(&Visited[d],0,1));
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

//reset visited elements
struct BF_Reset_F {
  int* Visited;
  BF_Reset_F(int* _Visited) : Visited(_Visited) {}
  inline bool operator() (uintE i){
    Visited[i] = 0;
    return 1;
  }
};

template <class vertex>
long CountFrontierEdges(vertex* G, long n, vertexSubset& Frontier) {
  long frontier_size = Frontier.numNonzeros();
  if (frontier_size == 0) return 0;

  if (Frontier.dense()) {
    long* Degrees = newA(long, n);
    parallel_for(long i=0;i<n;i++) {
      Degrees[i] = Frontier.isIn(i) ? G[i].getOutDegree() : 0;
    }
    long total_degree = sequence::plusReduce(Degrees, n);
    free(Degrees);
    return total_degree;
  }

  long* Degrees = newA(long, frontier_size);
  parallel_for(long i=0;i<frontier_size;i++) {
    Degrees[i] = G[Frontier.vtx(i)].getOutDegree();
  }
  long total_degree = sequence::plusReduce(Degrees, frontier_size);
  free(Degrees);
  return total_degree;
}

template <class vertex>
void Compute(hypergraph<vertex>& GA, commandLine P) {
  auto throughput_start = std::chrono::high_resolution_clock::now();
  long edges_processed = 0;
  long start = P.getOptionLongValue("-r",0);
  long nv = GA.nv, nh = GA.nh;
  //initialize ShortestPathLen to "infinity"
  intE* ShortestPathLenV = newA(intE,nv);
  intE* ShortestPathLenH = newA(intE,nh);
  {parallel_for(long i=0;i<nv;i++) ShortestPathLenV[i] = INT_MAX/2;}
  {parallel_for(long i=0;i<nh;i++) ShortestPathLenH[i] = INT_MAX/2;}
  ShortestPathLenV[start] = 0;

  int* Visited = newA(int,max(nv,nh));
  {parallel_for(long i=0;i<max(nv,nh);i++) Visited[i] = 0;}

  vertexSubset Frontier(nv,start); //initial frontier

  long round = 0;
  while(1){
    if(round == nv-1) {
      //negative weight cycle
      {parallel_for(long i=0;i<nv;i++) ShortestPathLenV[i] = -(INT_E_MAX/2);}
      break;
    }
    //cout << Frontier.numNonzeros() << endl;
    edges_processed += CountFrontierEdges(GA.V, GA.nv, Frontier);
    hyperedgeSubset output = vertexProp(GA, Frontier, BF_Relax_F(ShortestPathLenV,ShortestPathLenH,Visited),-1,dense_forward);
    hyperedgeMap(output,BF_Reset_F(Visited));
    Frontier.del();
    Frontier = output;
    if(Frontier.isEmpty()) break;
    //cout << Frontier.numNonzeros() << endl;

    // only count once for each csr edge instead of hyper edge
    // edges_processed += CountFrontierEdges(GA.H, GA.nh, Frontier);
    output = hyperedgeProp(GA, Frontier, BF_Relax_F(ShortestPathLenH,ShortestPathLenV,Visited),-1,dense_forward);
    vertexMap(output,BF_Reset_F(Visited));
    Frontier.del();
    Frontier = output;
    if(Frontier.isEmpty()) break;
    round++;
  }
  Frontier.del(); free(Visited);
  free(ShortestPathLenV); free(ShortestPathLenH);

  auto throughput_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime = throughput_end - throughput_start;
  std::cout << "Average runtime: " << runtime.count() << std::endl;
  std::cout << "Number of Processed Edges: " << edges_processed << std::endl;
  std::cout << "Processed Edges per Second: "
            << static_cast<double>(edges_processed) / runtime.count()
            << std::endl;
}
