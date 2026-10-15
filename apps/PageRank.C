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
#include "math.h"

// GACGE (2026-10-14): the same PageRank as Gunrock's pr (GACGE_PR=pull, its
// default): float ranks, pull (edgeMap threshold 0: the dense edgeMap sums each
// vertex's in-neighbors), contrib[u] = p[u] / deg(u) computed once per
// iteration, p'[v] = (1 - d) / n + d * sum of contrib over v's neighbors, no
// redistribution of dangling (degree-0) mass. With -epsilon 0 the L1 norm is
// computed for the last iteration only (no convergence test).
typedef float pr_t;

struct PR_F {
  pr_t* contrib, *p_next;
  PR_F(pr_t* _contrib, pr_t* _p_next) : contrib(_contrib), p_next(_p_next) {}
  inline bool update(uintE s, uintE d){ //update function applies PageRank equation
    p_next[d] += contrib[s];
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d) { //atomic Update
    writeAdd(&p_next[d],contrib[s]);
    return 1;
  }
  inline bool cond (intT d) { return cond_true(d); }};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
  pr_t damping;
  pr_t addedConstant;
  pr_t* p_next;
  PR_Vertex_F(pr_t* _p_next, pr_t _damping, intE n) :
    damping(_damping), addedConstant((1-_damping)*(1/(pr_t)n)), p_next(_p_next) {}
  inline bool operator () (uintE i) {
    p_next[i] = damping*p_next[i] + addedConstant;
    return 1;
  }
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long maxIters = P.getOptionLongValue("-maxiters",100);
  const intE n = GA.n;
  // -epsilon 0: no convergence test, exactly -maxiters iterations (GACGE runs
  // 20, as Gunrock's pr driver with GACGE_PR_ITERS=20)
  const pr_t damping = 0.85;
  const double epsilon = P.getOptionDoubleValue("-epsilon",0.0000001);
  vertex* V = GA.V;

  pr_t one_over_n = 1/(pr_t)n;
  pr_t* p_curr = newA(pr_t,n);
  {parallel_for(long i=0;i<n;i++) p_curr[i] = one_over_n;}
  pr_t* p_next = newA(pr_t,n);
  {parallel_for(long i=0;i<n;i++) p_next[i] = 0;}
  pr_t* contrib = newA(pr_t,n);
  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;}

  vertexSubset Frontier(n,n,frontier);

  long iter = 0;
  double L1 = 0;
  while(iter++ < maxIters) {
    {parallel_for(long i=0;i<n;i++) {
      const uintE d = V[i].getOutDegree();
      contrib[i] = d ? p_curr[i]/d : 0;
    }}
    edgeMap(GA,Frontier,PR_F(contrib,p_next),0, no_output);
    vertexMap(Frontier,PR_Vertex_F(p_next,damping,n));
    if (epsilon > 0 || iter == maxIters) {
      //L1-norm between p_curr and p_next (contrib is free until the next iteration)
      {parallel_for(long i=0;i<n;i++) contrib[i] = fabs(p_curr[i]-p_next[i]);}
      L1 = sequence::plusReduce(contrib,n);
    }
    swap(p_curr,p_next);
    if (epsilon > 0 && L1 < epsilon) break;
    {parallel_for(long i=0;i<n;i++) p_next[i] = 0;}
  }
  // p_curr holds the ranks of the last iteration
  const long done = iter > maxIters ? maxIters : iter;
  double sum = sequence::plusReduce(p_curr,n);
  std::printf("PRResult: iterations=%ld L1=%g rank_sum=%.9f rank0=%.9e\n", done, L1, sum, (double)p_curr[0]);
  Frontier.del(); free(p_curr); free(p_next); free(contrib);
}
