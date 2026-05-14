#pragma once

#include "Complex.hpp"
#include "ComplexTable.hpp"

#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <cstddef>
#include <vector>
#include <utility>

namespace dd {

template <std::size_t INITIAL_ALLOCATION_SIZE = 2048,
          std::size_t GROWTH_FACTOR = 2>
class ComplexCache {
  using Entry = ComplexTable<>::Entry;
  // using ComplexKey = std::pair<fp*, fp*>;
  // // Custom hash function for ComplexKey
  using ComplexKey = std::pair<fp, fp>;
  using ActiveKey = std::pair<Entry*, Entry*>;

  // Hash for value-based complex key
  struct ComplexKeyHash {
    std::size_t operator()(const ComplexKey& key) const {
      // auto hash1 = std::hash<fp*>{}(key.first);
      // auto hash2 = std::hash<fp*>{}(key.second);
      auto hash1 = std::hash<fp>{}(key.first);
      auto hash2 = std::hash<fp>{}(key.second);
      return hash1 ^ (hash2 << 1);  // Shift and XOR for combining hash values
    }
  };

  // Hash for active pointer-pair bookkeeping
  struct ActiveKeyHash {
    std::size_t operator()(const ActiveKey& key) const {
      auto hash1 = std::hash<Entry*>{}(key.first);
      auto hash2 = std::hash<Entry*>{}(key.second);
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
      // // std::cout << "53:  get Cached Complex function in: " << entry.r << " " << entry.i << " " << available << std::endl;
      // complexMap.insert({{&entry.r->value, &entry.i->value}, true});
      peakCount = std::max(peakCount, count);
      activeComplexes.insert({entry.r, entry.i});
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
    // // std::cout << "74:"<< c.r->value <<","<< c.i->value<<" get Cached Complex function in: " << &c<< std::endl;
    // complexMap.insert({{&c.r->value, &c.i->value}, true});
    peakCount = std::max(peakCount, count);
    activeComplexes.insert({c.r, c.i});
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
    assert(count >= 2);
    assert(c != Complex::zero);
    assert(c != Complex::one);
    assert(c.r->refCount == 0);
    assert(c.i->refCount == 0);
    activeComplexes.erase({c.r, c.i});
    c.i->next = available;
    c.r->next = c.i;
    available = c.r;
    count -= 2;
    // complexMap.erase({&c.r->value, &c.i->value});
    // // Remove these debugging lines later!!
    // if(available->next->next == available) {
    //   std::cout << available << " " << available->next << " " << std::endl;
    //   assert(1 == 0);
    // }
  }

  bool isInCache(fp real, fp imag) {
    // for (const auto& [key, _] : complexMap) {
    //   if (*key.first == real && *key.second == imag) {
    //     return true;
    //   }
    // }
    // return false;
    rebuildValueIndex();
    const ComplexKey key{real, imag};
    const auto it = valueIndex.find(key);
    return it != valueIndex.end() && it->second > 0;
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
    // complexMap.clear();
    activeComplexes.clear();
    valueIndex.clear();
  };

private:
  void rebuildValueIndex() {
    valueIndex.clear();
    for (const auto& key : activeComplexes) {
      const ComplexKey valueKey{key.first->value, key.second->value};
      auto it = valueIndex.find(valueKey);
      if (it == valueIndex.end()) {
        valueIndex.emplace(valueKey, 1U);
      } else {
        ++(it->second);
      }
    }
  }

  Entry* available{};
  std::vector<std::vector<Entry>> chunks{};
  std::size_t chunkID{0};
  typename std::vector<Entry>::iterator chunkIt;
  typename std::vector<Entry>::iterator chunkEndIt;
  std::size_t allocationSize;

  std::size_t allocations = 0;
  std::size_t count = 0;
  std::size_t peakCount = 0;
  // std::unordered_map<ComplexKey, bool, ComplexKeyHash> complexMap;
  std::unordered_set<ActiveKey, ActiveKeyHash> activeComplexes;
  std::unordered_map<ComplexKey, std::size_t, ComplexKeyHash> valueIndex;
};
} // namespace dd
