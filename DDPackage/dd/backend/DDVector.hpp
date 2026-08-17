#pragma once

// QReach backend — DDVector (quantum state vectors). Contract §2.
//
// Index convention (see the "observable contract" in the API doc §6.2):
//   - a vector over n qubits has n index entries, key "q<qubit>" with idx 0;
//   - qubit 0 is the most significant (big-endian): basis integer
//       index = sum_{q} bit_q * 2^(n-1-q);
//   - the dense amplitude array is row-major with dimension d = qubit d, so the
//     flat array position IS the big-endian basis integer.
//
// Every vector is built through dd::Tensor::to_tdd, which fixes varOrder /
// reOrder / key_2_index consistently, so gate·vector contraction (DDMatrix)
// will match keys exactly.

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "dd/Tensor.hpp"
#include "dd/backend/DDTypes.hpp"

namespace DDVector {

using limtdd::DD;
using limtdd::DDComplex;

// Idempotent. Funnels into limtdd::Initialize().
inline void Initialize() { limtdd::Initialize(); }

namespace detail {

// Dense-array / enumeration are inherently O(2^n); guard against overflow.
constexpr unsigned int kMaxQubits = 31;
// Compact constructions (basis/constant states) are O(n) in memory, so they
// may go well past kMaxQubits. Bound by Initialize()'s key pre-registration
// (kMaxKey = 256) so the interleaved variable-order invariant still holds.
constexpr unsigned int kMaxCompactQubits = 256;

inline std::vector<dd::Index> qubitIndices(unsigned int n) {
    std::vector<dd::Index> idx;
    idx.reserve(n);
    for (unsigned int q = 0; q < n; ++q) {
        idx.push_back({"q" + std::to_string(q), 0});
    }
    return idx;
}

// Build an n-qubit basis state |index> (big-endian) as a tensor product of n
// single-qubit basis tensors via Package::cont. Compact O(n) — no dense 2^n
// array. TensorNetwork::cont returns an edge with refcount 1; release it back
// to 0 so the DD constructor adopts it with a single reference.
inline DD basisState(unsigned int n, unsigned int index) {
    if (n == 0) {
        // 0-qubit "basis state" is the scalar 1 (index must be 0, enforced by caller).
        dd::TDD res;
        res.e = dd::Edge<dd::mNode>::one;
        return DD(std::move(res));
    }
    auto& pkg = limtdd::backendPackage();
    dd::TensorNetwork tn;
    for (unsigned int q = 0; q < n; ++q) {
        const bool bit = ((index >> (n - 1 - q)) & 1u) != 0u;
        xt::xarray<dd::ComplexValue> arr(std::vector<std::size_t>{2});
        arr[0] = dd::ComplexValue{bit ? 0.0 : 1.0, 0.0};
        arr[1] = dd::ComplexValue{bit ? 1.0 : 0.0, 0.0};
        tn.add_ts(dd::Tensor(arr, {{"q" + std::to_string(q), 0}}, "basis"));
    }
    dd::TDD res = tn.cont(&pkg);
    pkg.decRef(res.e);
    return DD(std::move(res));
}

// Build an n-qubit all-ones vector (every amplitude 1) as a tensor product of
// n single-qubit [1,1] tensors. Compact O(n).
inline DD allOnesState(unsigned int n) {
    if (n == 0) {
        // 0-qubit all-ones vector is the scalar 1.
        dd::TDD res;
        res.e = dd::Edge<dd::mNode>::one;
        return DD(std::move(res));
    }
    auto& pkg = limtdd::backendPackage();
    dd::TensorNetwork tn;
    for (unsigned int q = 0; q < n; ++q) {
        xt::xarray<dd::ComplexValue> arr(std::vector<std::size_t>{2});
        arr[0] = dd::ComplexValue{1.0, 0.0};
        arr[1] = dd::ComplexValue{1.0, 0.0};
        tn.add_ts(dd::Tensor(arr, {{"q" + std::to_string(q), 0}}, "ones"));
    }
    dd::TDD res = tn.cont(&pkg);
    pkg.decRef(res.e);
    return DD(std::move(res));
}

// Recursively walk the DD and collect every non-zero (index, amplitude).
// Big-endian: slicing variable v with value c sets bit (n-1-v) of the index.
inline void enumerate(const dd::Edge<dd::mNode>& e, unsigned int n,
                      unsigned int index,
                      std::vector<std::pair<unsigned int, DDComplex>>& out) {
    if (e.p == nullptr) {
        return;
    }
    const auto w = limtdd::fromComplex(e.w);
    if (w.real() == 0.0 && w.imag() == 0.0) {
        return;  // exact zero subtree — contributes nothing
    }
    if (e.p->v == -1) {  // terminal
        out.emplace_back(index, w);
        return;
    }
    const int v = static_cast<int>(e.p->v);
    for (int c = 0; c < 2; ++c) {
        const auto child = limtdd::backendPackage().backendSlice(e, v, c);
        enumerate(child, n, index | (static_cast<unsigned int>(c) << (n - 1 - v)),
                  out);
    }
}

inline std::vector<std::pair<unsigned int, DDComplex>> enumerateAll(const DD& c) {
    std::vector<std::pair<unsigned int, DDComplex>> out;
    const unsigned int n = static_cast<unsigned int>(c.tdd.index_set.size());
    enumerate(c.tdd.e, n, 0u, out);
    return out;
}

// Sum over all remaining variables of a scalar DD. Package::cont leaves a
// residual node (e.g. children [+1, -1]) for an all-variables-contracted zero
// result instead of collapsing to a terminal, so "exactly one non-zero leaf"
// is wrong for the zero inner-product case. Recursively summing the branches
// yields the correct scalar.
inline DDComplex sumRemaining(const dd::Edge<dd::mNode>& e) {
    if (e.p == nullptr) {
        return DDComplex(0.0, 0.0);
    }
    const auto w = limtdd::fromComplex(e.w);
    if (w.real() == 0.0 && w.imag() == 0.0) {
        return DDComplex(0.0, 0.0);
    }
    if (e.p->v == -1) {
        return w;
    }
    const int v = static_cast<int>(e.p->v);
    const auto b0 = limtdd::backendPackage().backendSlice(e, v, 0);
    const auto b1 = limtdd::backendPackage().backendSlice(e, v, 1);
    return sumRemaining(b0) + sumRemaining(b1);
}

}  // namespace detail

inline DD MkBasisVector(unsigned int n, unsigned int index) {
    if (n > detail::kMaxCompactQubits) {
        throw std::invalid_argument("MkBasisVector: too many qubits");
    }
    // `index` is 32-bit, so it can only address n <= 32; skip the (otherwise
    // overflowing) shift check beyond that.
    if (n <= 32 && index >= (std::size_t{1} << n)) {
        throw std::invalid_argument("MkBasisVector: index out of range");
    }
    return detail::basisState(n, index);
}

inline DD MkBasisVector(unsigned int n, std::string bitstring) {
    if (bitstring.size() != n) {
        throw std::invalid_argument("MkBasisVector: bitstring length != n");
    }
    unsigned int index = 0;
    for (const char ch : bitstring) {
        index = (index << 1) | (ch == '1' ? 1u : 0u);
    }
    return MkBasisVector(n, index);
}

inline DD NoDistinctionNode(unsigned int n, DDComplex val) {
    if (n > detail::kMaxCompactQubits) {
        throw std::invalid_argument("NoDistinctionNode: too many qubits");
    }
    return val * detail::allOnesState(n);
}

inline DD InitializeWithAmplitudes(unsigned int qnum, std::vector<double> amps) {
    if (qnum > detail::kMaxQubits) {
        throw std::invalid_argument("InitializeWithAmplitudes: too many qubits");
    }
    const std::size_t dim = std::size_t{1} << qnum;
    if (amps.size() != 2 * dim) {
        throw std::invalid_argument("InitializeWithAmplitudes: amps size != 2 * 2^qnum");
    }
    std::vector<std::size_t> shape(qnum, 2);
    xt::xarray<dd::ComplexValue> arr(shape);
    for (std::size_t i = 0; i < dim; ++i) {
        arr.data()[i] = dd::ComplexValue{amps[i], amps[dim + i]};
    }
    dd::Tensor tensor(arr, detail::qubitIndices(qnum), "amplitudes");
    return DD(tensor.to_tdd(&limtdd::backendPackage()));
}

// No-op: this backend multiplies gate matrices directly against plain vectors.
inline DD VectorToMatrixInterleaved(DD vec) { return vec; }

inline int GetLevel(const DD& c) {
    // The first constructor argument is now the real qubit count n (QReach no
    // longer pads to a power of 2). Return n directly so SingleVecTerm can
    // derive qNum without the lossy ceil(log2(n)).
    return static_cast<int>(c.tdd.index_set.size());
}

inline bool IsApproximatelyZero(const DD& c, double threshold = 1e-8) {
    // A canonical DD is zero iff its root weight is zero — check the weight
    // directly and never slice (Slicing asserts on zero edges).
    const auto w = limtdd::fromComplex(c.tdd.e.w);
    return std::abs(w.real()) <= threshold && std::abs(w.imag()) <= threshold;
}

inline DDComplex ExtractSingleAmplitude(const DD& c) {
    // Sum over any residual variables left by Package::cont. This returns the
    // scalar value for dot()/normalize() results, including the zero case
    // (orthogonal states). sumRemaining short-circuits zero edges before
    // slicing.
    return detail::sumRemaining(c.tdd.e);
}

inline std::vector<std::pair<unsigned int, DDComplex>> GetNonZeroAmplitudes(
    const DD& c, double threshold = 1e-8) {
    // Short-circuit a zero DD (root weight zero) so we never slice it.
    if (IsApproximatelyZero(c, 0.0)) {
        return {};
    }
    auto amps = detail::enumerateAll(c);
    std::vector<std::pair<unsigned int, DDComplex>> out;
    out.reserve(amps.size());
    for (const auto& [index, amp] : amps) {
        if (std::abs(amp.real()) > threshold || std::abs(amp.imag()) > threshold) {
            out.emplace_back(index, amp);
        }
    }
    return out;
}

inline DD Normalize(const DD& c) {
    double norm2 = 0.0;
    for (const auto& [index, amp] : detail::enumerateAll(c)) {
        (void)index;
        norm2 += amp.norm();
    }
    if (norm2 <= 0.0) {
        return c;  // zero vector — nothing to normalize
    }
    return DDComplex(1.0 / std::sqrt(norm2), 0.0) * c;
}

// Inner product ⟨a|b⟩ = conj(a)·b = Σ aᵢ* bᵢ.
//
// Implemented as a sparse dot product over the two vectors' non-zero
// amplitudes (joined on the big-endian basis index), NOT via Package::cont.
// Package::cont's fully-contracted (scalar) path leaves a spurious residual
// node and scales the root weight by a data-dependent factor — see
// docs/limtdd-innerproduct-bug-report.md. Concretely <00|00> came out as 2 and
// <v1|v3> as -3·2³·⟨v1|v3⟩, so the Gram-Schmidt coefficients in QReach's
// disjunction/span_qops were wrong and linearly-dependent vectors were kept.
//
// Cost is O(#non-zero amplitudes) — the same order as the dense conjugation
// this function previously performed, so no regression; a DD-level contraction
// (fixing cont's scalar path) is the future path to sub-exponential cost.
inline DDComplex InnerProduct(DD a, DD b) {
    const auto va = GetNonZeroAmplitudes(a, 0.0);
    const auto vb = GetNonZeroAmplitudes(b, 0.0);
    std::unordered_map<unsigned int, DDComplex> bmap;
    bmap.reserve(vb.size());
    for (const auto& [idx, val] : vb) {
        bmap[idx] = val;
    }
    DDComplex acc(0.0, 0.0);
    for (const auto& [idx, aval] : va) {
        const auto it = bmap.find(idx);
        if (it == bmap.end()) {
            continue;
        }
        const auto& bval = it->second;
        // conj(a)·b = (a.re·b.re + a.im·b.im) + i(a.re·b.im - a.im·b.re)
        acc += DDComplex(aval.real() * bval.real() + aval.imag() * bval.imag(),
                         aval.real() * bval.imag() - aval.imag() * bval.real());
    }
    return acc;
}

inline void VectorPrintColumnHead(const DD& c, std::ostream& out) {
    for (const auto& [index, amp] : GetNonZeroAmplitudes(c)) {
        out << '|' << index << ">\t" << amp.real() << (amp.imag() >= 0 ? " + " : " - ")
            << std::abs(amp.imag()) << "i\n";
    }
}

}  // namespace DDVector
