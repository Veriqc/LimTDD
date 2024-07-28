#pragma once

#include "Complex.hpp"
#include "ComplexTable.hpp"

#include <unordered_map>
#include <cassert>
#include <cstddef>
#include <vector>
#include <utility>

namespace dd {

template <std::size_t INITIAL_ALLOCATION_SIZE = 2048,
          std::size_t GROWTH_FACTOR = 2>
class ComplexCache {
  using Entry = ComplexTable<>::Entry;
  using ComplexKey = std::pair<fp*, fp*>;
  // Custom hash function for ComplexKey
  struct ComplexKeyHash {
    std::size_t operator()(const ComplexKey& key) const {
      auto hash1 = std::hash<fp*>{}(key.first);
      auto hash2 = std::hash<fp*>{}(key.second);
      return hash1 ^ (hash2 << 1);  // Shift and XOR for combining hash values
    }
  };

public:
  ComplexCache() : allocationSize(INITIAL_ALLOCATION_SIZE) {
    // allocate first chunk of cache entries
    chunks.emplace_back(allocationSize);
    allocations += allocationSize;
    allocationSize *= GROWTH_FACTOR;
    chunkIt = chunks[0].begin();
    chunkEndIt = chunks[0].end();
  }

  ~ComplexCache() = default;

  // access functions
  [[nodiscard]] std::size_t getCount() const { return count; }
  [[nodiscard]] std::size_t getPeakCount() const { return peakCount; }
  [[nodiscard]] std::size_t getAllocations() const { return allocations; }
  [[nodiscard]] std::size_t getGrowthFactor() const { return GROWTH_FACTOR; }

  [[nodiscard]] Complex getCachedComplex() {
    // an entry is available on the stack
    if (available != nullptr) {
      assert(available->next != nullptr);
      auto entry = Complex{available, available->next};
      available = entry.i->next;
      count += 2;
      std::cout << "53: "<< entry.r->value <<","<< entry.i->value << " get Cached Complex function in: " << &entry << std::endl;
      complexMap.insert({{&entry.r->value, &entry.i->value}, true});
      return entry;
    }

    // new chunk has to be allocated
    if (chunkIt == chunkEndIt) {
      chunks.emplace_back(allocationSize);
      allocations += allocationSize;
      allocationSize *= GROWTH_FACTOR;
      chunkID++;
      chunkIt = chunks[chunkID].begin();
      chunkEndIt = chunks[chunkID].end();
    }

    Complex c{};
    c.r = &(*chunkIt);
    ++chunkIt;
    c.i = &(*chunkIt);
    ++chunkIt;
    count += 2;
    std::cout << "74:"<< c.r->value <<","<< c.i->value<<" get Cached Complex function in: " << &c<< std::endl;
    complexMap.insert({{&c.r->value, &c.i->value}, true});
    return c;
  }

  [[nodiscard]] Complex getTemporaryComplex() {
    // an entry is available on the stack
    if (available != nullptr) {
      assert(available->next != nullptr);
      return {available, available->next};
    }

    // new chunk has to be allocated
    if (chunkIt == chunkEndIt) {
      chunks.emplace_back(allocationSize);
      allocations += allocationSize;
      allocationSize *= GROWTH_FACTOR;
      chunkID++;
      chunkIt = chunks[chunkID].begin();
      chunkEndIt = chunks[chunkID].end();
    }
    return {&(*chunkIt), &(*(chunkIt + 1))};
  }

  void returnToCache(Complex& c) {
    std::cout << c << " return to cache in " << &c << std::endl;
    assert(count >= 2);
    assert(c != Complex::zero);
    assert(c != Complex::one);
    assert(c.r->refCount == 0);
    assert(c.i->refCount == 0);
    c.i->next = available;
    c.r->next = c.i;
    available = c.r;
    count -= 2;
    complexMap.erase({&c.r->value, &c.i->value});
  }

  bool isInCache(fp real, fp imag) {
    for (const auto& [key, _] : complexMap) {
      if (*key.first == real && *key.second == imag) {
        return true;
      }
    }
    return false;
  }

  void clear() {
    // clear available stack
    available = nullptr;

    // release memory of all but the first chunk TODO: it could be desirable to
    // keep the memory
    while (chunkID > 0) {
      chunks.pop_back();
      chunkID--;
    }
    // restore initial chunk setting
    chunkIt = chunks[0].begin();
    chunkEndIt = chunks[0].end();
    allocationSize = INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR;
    allocations = INITIAL_ALLOCATION_SIZE;

    count = 0;
    peakCount = 0;
    complexMap.clear();
  };

private:
  Entry* available{};
  std::vector<std::vector<Entry>> chunks{};
  std::size_t chunkID{0};
  typename std::vector<Entry>::iterator chunkIt;
  typename std::vector<Entry>::iterator chunkEndIt;
  std::size_t allocationSize;

  std::size_t allocations = 0;
  std::size_t count = 0;
  std::size_t peakCount = 0;
  std::unordered_map<ComplexKey, bool, ComplexKeyHash> complexMap;
};
} // namespace dd
