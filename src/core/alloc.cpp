#include "slamcore/core/alloc.hpp"

// PoolAllocator and RingAllocator are header-only; this translation unit
// exists so the core library has a stable symbol for the allocator module and
// to host the out-of-line global budget check used by the app layer.

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace sc {

namespace {

Size parseVmHWM() {
  std::ifstream f("/proc/self/status");
  std::string key;
  while (f >> key) {
    if (key == "VmHWM:") {
      Size kb = 0;
      f >> kb;
      return kb * 1024;
    }
    std::string rest;
    std::getline(f, rest);
  }
  return 0;
}

}  // namespace

/// Resident-set high-water mark in bytes, or 0 if /proc is unavailable.
Size residentHighWaterBytes() { return parseVmHWM(); }

}  // namespace sc
