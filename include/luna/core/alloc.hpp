// LUNA -- deterministic allocation discipline.
//
// The 8 GB budget is enforced structurally, not by hope:
//   * PoolAllocator  -- bump allocator over a pre-reserved arena; every
//                       allocation in a real-time loop is O(1) and never
//                       touches the global heap.
//   * RingAllocator  -- fixed-capacity FIFO arena for per-frame scratch that
//                       must survive exactly one lag window.
//
// Both are noexcept-on-exhaustion (they throw once, at init, if the arena
// cannot be reserved) so hot loops contain no allocation failure branches.
#pragma once

#include <cstddef>
#include <new>
#include <vector>

#include "luna/core/types.hpp"

namespace luna {

/// Bump allocator over a caller-reserved arena. Not thread-safe by design:
/// each real-time thread owns one.
class PoolAllocator {
 public:
  explicit PoolAllocator(Size bytes) : buf_(bytes), head_(0) {}
  PoolAllocator(const PoolAllocator&) = delete;
  PoolAllocator& operator=(const PoolAllocator&) = delete;

  /// Aligns to `align` and returns `n` bytes, or nullptr if exhausted.
  void* allocate(Size n, Size align = alignof(std::max_align_t)) noexcept {
    const Size mask = align - 1;
    const Size start = (head_ + mask) & ~mask;
    if (start + n > buf_.size()) {
      ++misses_;
      return nullptr;
    }
    head_ = start + n;
    ++allocs_;
    if (head_ > high_water_) high_water_ = head_;
    return buf_.data() + start;
  }

  void reset() noexcept { head_ = 0; }

  Size used() const noexcept { return head_; }
  Size capacity() const noexcept { return buf_.size(); }
  Size highWater() const noexcept { return high_water_; }
  Size allocs() const noexcept { return allocs_; }
  Size misses() const noexcept { return misses_; }
  Scalar utilisation() const noexcept {
    return buf_.empty() ? 0.0 : static_cast<Scalar>(high_water_) / static_cast<Scalar>(buf_.size());
  }

 private:
  std::vector<std::byte> buf_;
  Size head_;
  Size high_water_ = 0;
  Size allocs_ = 0;
  Size misses_ = 0;
};

/// Fixed-capacity ring of trivially-copyable, fixed-stride slots.
class RingAllocator {
 public:
  RingAllocator(Size slot_bytes, Size slots)
      : buf_(slot_bytes * slots), slot_bytes_(slot_bytes), nslots_(slots) {}

  /// Oldest slot, which is recycled on wrap. Never fails after construction.
  void* next() noexcept {
    void* p = buf_.data() + slot_ * slot_bytes_;
    slot_ = (slot_ + 1) % nslots_;
    if (++filled_ > nslots_) filled_ = nslots_;
    return p;
  }

  Size slots() const noexcept { return nslots_; }
  Size live() const noexcept { return filled_; }
  Size bytes() const noexcept { return buf_.size(); }

 private:
  std::vector<std::byte> buf_;
  Size slot_bytes_;
  Size nslots_;
  Size slot_ = 0;
  Size filled_ = 0;
};

}  // namespace luna
