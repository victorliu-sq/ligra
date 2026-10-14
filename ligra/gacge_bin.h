// GACGE binary CSR (after GraphRTX's GRTXCSR, RT-Benchmarks). One canonical
// symmetric graph per dataset, read by Ligra and Gunrock alike:
//   <name>.bin    header (48 B) | uint32 row_ptr[n + 1] | uint32 col[m]
//   <name>.w.bin  header (48 B) | uint32 weight[m]      (SSSP: in [1, 64], the
//                 same splitmix64 of the edge's end points as the text files)
// m counts directed entries (both directions of every edge); every neighbor
// list is sorted, without self-loops or duplicates. Little-endian, as written.
// Written by data/targets/scripts/bin/xbcsr_to_bin.cpp; a copy of this header
// is src/gunrock/include/gunrock/io/gacge_bin.hxx.
#pragma once

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace gacge_bin {

struct Header {
  char magic[8];      // "GACGECSR"
  uint32_t version;   // 1
  uint32_t kind;      // 0 graph, 1 weights
  uint64_t n;         // vertices
  uint64_t m;         // directed entries
  uint64_t row_count; // n + 1 (graph), 0 (weights)
  uint64_t checksum;  // FNV-1a 64 of the arrays after the header (0: none)
};
static_assert(sizeof(Header) == 48, "GACGE binary header is 48 bytes");

enum : uint32_t { kGraph = 0, kWeights = 1 };
static const char kMagic[8] = {'G', 'A', 'C', 'G', 'E', 'C', 'S', 'R'};

inline bool is_bin(const char* path) {
  const size_t l = std::strlen(path);
  return l > 4 && std::strcmp(path + l - 4, ".bin") == 0;
}

// <name>.bin -> <name>.w.bin
inline std::string weights_path(const std::string& graph) {
  return graph.substr(0, graph.size() - 4) + ".w.bin";
}

[[noreturn]] inline void fail(const std::string& path, const char* why) {
  std::fprintf(stderr, "gacge_bin: %s: %s\n", path.c_str(), why);
  std::abort();
}

// A read-only mapping of a GACGE file; a (row_ptr or weights) and b (col) point
// into it. The caller copies what it needs and unmaps.
struct Mapped {
  Header h{};
  const uint32_t* a = nullptr;
  const uint32_t* b = nullptr;
  void* base = nullptr;
  size_t len = 0;
};

inline Mapped map_file(const std::string& path, uint32_t kind) {
  Mapped f;
  const int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0) fail(path, "cannot open");
  struct stat st;
  if (fstat(fd, &st) != 0) fail(path, "cannot stat");
  f.len = static_cast<size_t>(st.st_size);
  if (f.len < sizeof(Header)) fail(path, "shorter than the header");
  f.base = mmap(nullptr, f.len, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (f.base == MAP_FAILED) fail(path, "mmap failed");
  madvise(f.base, f.len, MADV_WILLNEED);
  std::memcpy(&f.h, f.base, sizeof(Header));
  if (std::memcmp(f.h.magic, kMagic, 8) != 0) fail(path, "bad magic (not a GACGE binary graph)");
  if (f.h.version != 1) fail(path, "unsupported version");
  if (f.h.kind != kind) fail(path, kind == kGraph ? "not a graph file" : "not a weights file");
  if (f.h.n > UINT32_MAX || f.h.m > UINT32_MAX) fail(path, "n or m beyond 32 bits");
  const uint64_t words = (kind == kGraph) ? f.h.row_count + f.h.m : f.h.m;
  if (kind == kGraph && f.h.row_count != f.h.n + 1) fail(path, "row_count != n + 1");
  if (f.len != sizeof(Header) + 4 * words) fail(path, "file size does not match the header");
  const uint32_t* p = reinterpret_cast<const uint32_t*>(static_cast<const char*>(f.base) + sizeof(Header));
  if (kind == kGraph) {
    f.a = p;
    f.b = p + f.h.row_count;
    if (f.a[0] != 0 || f.a[f.h.n] != f.h.m) fail(path, "row_ptr does not run from 0 to m");
  } else {
    f.a = p;
  }
  return f;
}

inline void unmap(Mapped& f) {
  if (f.base) munmap(f.base, f.len);
  f = Mapped();
}

inline uint64_t fnv1a(const void* data, size_t bytes, uint64_t h = 1469598103934665603ULL) {
  const unsigned char* p = static_cast<const unsigned char*>(data);
  for (size_t i = 0; i < bytes; i++) {
    h ^= p[i];
    h *= 1099511628211ULL;
  }
  return h;
}

}  // namespace gacge_bin
