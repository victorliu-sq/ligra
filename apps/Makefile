ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

ifdef PD
PD = -DPD
endif

ifdef BYTE
CODE = -DBYTE
else ifdef NIBBLE
CODE = -DNIBBLE
else
CODE = -DBYTERLE
endif

ifdef LOWMEM
MEM = -DLOWMEM
endif

#compilers
ifdef CILK
PCC = g++
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -O3 -DCILK $(INTT) $(INTE) $(CODE) $(PD) $(MEM)
PLFLAGS = -fcilkplus -lcilkrts

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -std=c++14 -O3 -DCILKP $(INTT) $(INTE) $(CODE) $(PD) $(MEM)

else
PCC = g++
PCFLAGS = -std=c++14 -O3 $(INTT) $(INTE) $(CODE) $(PD) $(MEM)
endif

COMMON= ligra.h graph.h compressedVertex.h vertex.h utils.h IO.h parallel.h gettime.h index_map.h maybe.h sequence.h edgeMap_utils.h binary_search.h quickSort.h blockRadixSort.h transpose.h parseCommandLine.h byte.h byteRLE.h nibble.h byte-pd.h byteRLE-pd.h nibble-pd.h vertexSubset.h encoder.C decoder.C

ALL= encoder decoder BFS BC SSSP SSSP_Throughput MultiSSSP MultiSSSP_Throughput BellmanFord Components Components-Shortcut Radii PageRank PageRankDelta BFSCC BFS-Bitvector KCore MIS Triangle CF

all: $(ALL)

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

# <app>_cilk: the app built with OpenCilk (clang -fopencilk, work-stealing
# runtime), the build GACGE uses; OPENCILK_HOME is the install (deps/opencilk).
OPENCILK_HOME ?= ../../../deps/opencilk
CILKCXX = $(OPENCILK_HOME)/bin/clang++ -std=c++14 -fopencilk -march=native -O3 -DOPENCILK $(INTT) $(INTE) $(CODE) $(PD) $(MEM)
%_cilk : %.C $(COMMON)
	$(CILKCXX) -o $@ $<
# <app>_cilkscale: Cilkscale build, the work and span of every Compute() (ligra.h)
%_cilkscale : %.C $(COMMON)
	$(CILKCXX) -fcilktool=cilkscale -DGACGE_CILKSCALE -o $@ $<

# XBlossomPro variants for the reuse measurement: XB-Ligra (no tree reuse) and
# XB-Pro-Ligra with reused / reset tree-node counters.
XBlossomPro_noreuse_cilk : XBlossomPro.C $(COMMON)
	$(CILKCXX) -DXBP_REUSE=0 -o $@ $<
XBlossomPro_counters_cilk : XBlossomPro.C $(COMMON)
	$(CILKCXX) -DXBP_COUNTERS=1 -o $@ $<
# XB-Pro's own design (path mutex, per-edge allocations, serial copies, its
# path table) for measure/xb-cpu-opt and xb-llc; XBlossomPro_cilk is the
# optimized default (XBP_OPT, XBP_TRIPLET).
XBlossomPro_xbpro_cilk : XBlossomPro.C $(COMMON)
	$(CILKCXX) -DXBP_OPT=0 -o $@ $<
# XB-Ligra with XB-Pro's path table (a vector per node, blossoms written out)
# instead of the default triplet records (XBP_TRIPLET=0): the design
# comparisons (measure/xb-cpu-opt, xb-llc, xb-ligra-triplet).
XBlossomPro_vector_cilk : XBlossomPro.C $(COMMON)
	$(CILKCXX) -DXBP_TRIPLET=0 -o $@ $<
# XB-Ligra letting Ligra run large steps dense (XBP_DENSE), not the default:
# measure/xb-ligra-dense.
XBlossomPro_dense_cilk : XBlossomPro.C $(COMMON)
	$(CILKCXX) -DXBP_DENSE=1 -o $@ $<
XBlossomPro_dense_cilkscale : XBlossomPro.C $(COMMON)
	$(CILKCXX) -fcilktool=cilkscale -DGACGE_CILKSCALE -DXBP_DENSE=1 -o $@ $<

$(COMMON):
	ln -s ../ligra/$@ .

.PHONY : clean

clean :
	rm -f *.o $(ALL)

cleansrc :
	rm -f *.o $(ALL)
	rm $(COMMON)
