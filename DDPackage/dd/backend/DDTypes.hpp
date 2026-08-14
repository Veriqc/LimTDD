#pragma once

// QReach backend — core types.
//
// This is the LimTDD side of the "DD backend replacement" API contract
// (see LimTDD/docs/backend-replacement-api-contract.md). It defines the two
// public value types (`DDComplex`, `DD`) plus the package singleton and the
// idempotent initialization hook. The `DDVector` / `DDMatrix` namespaces are
// in the sibling headers DDVector.hpp / DDMatrix.hpp.
//
// Internal variable order is our own choice (the contract only pins down the
// observable big-endian basis/bitstring convention). Every gate matrix and
// every vector is built through Package::xarray_2_edge / cont with a fixed
// index naming scheme, so gate construction and MatrixMultiplyWithVector stay
// mutually consistent.
//
// Memory model: LimTDD uses explicit reference counting (Package::incRef /
// decRef / garbageCollect over a unique table). `DD` wraps that convention in
// RAII — a `DD` owns exactly one reference to its `tdd.e`; copying increments
// it, destroying decrements it. Do NOT shallow-copy a dd::TDD and let it dangle.

#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

#include "dd/Package.hpp"
#include "dd/Tdd.hpp"

namespace limtdd {

// ---------------------------------------------------------------------------
// DDComplex — the scalar amplitude type (contract §1 + §6.3).
//
// Wraps std::complex<double> and satisfies every member the QReach semantic
// layer requires: real()/imag()/abs()/norm(), + - * /, ==/!=, comparison
// against integer literals 0 and 1, and construction from (re, im) / double.
// ---------------------------------------------------------------------------
struct DDComplex {
    std::complex<double> v{0.0, 0.0};

    DDComplex() = default;
    DDComplex(double re, double im) : v(re, im) {}
    DDComplex(double re) : v(re, 0.0) {}  // non-explicit: allows `amp = 0.0`
    DDComplex(const std::complex<double>& c) : v(c) {}

    [[nodiscard]] double real() const { return v.real(); }
    [[nodiscard]] double imag() const { return v.imag(); }
    [[nodiscard]] double abs() const { return std::abs(v); }
    [[nodiscard]] double norm() const { return std::norm(v); }

    DDComplex& operator+=(const DDComplex& o) { v += o.v; return *this; }
    DDComplex& operator-=(const DDComplex& o) { v -= o.v; return *this; }
    DDComplex& operator*=(const DDComplex& o) { v *= o.v; return *this; }
    DDComplex& operator/=(const DDComplex& o) { v /= o.v; return *this; }

    friend DDComplex operator+(DDComplex a, const DDComplex& b) { a += b; return a; }
    friend DDComplex operator-(DDComplex a, const DDComplex& b) { a -= b; return a; }
    friend DDComplex operator*(DDComplex a, const DDComplex& b) { a *= b; return a; }
    friend DDComplex operator/(DDComplex a, const DDComplex& b) { a /= b; return a; }

    friend DDComplex operator*(DDComplex a, double s) { a.v *= s; return a; }
    friend DDComplex operator*(double s, DDComplex a) { a.v *= s; return a; }

    bool operator==(const DDComplex& o) const { return v == o.v; }
    bool operator!=(const DDComplex& o) const { return v != o.v; }
    bool operator==(int rhs) const { return v == std::complex<double>(rhs, 0.0); }
    bool operator!=(int rhs) const { return !(*this == rhs); }
};

// Forward declaration — defined below (after DD, which uses it in its
// copy/move/destroy operations).
dd::Package<>& backendPackage();

// ---------------------------------------------------------------------------
// DD — the vector/matrix value type (contract §1).
//
// RAII wrapper over dd::TDD. Copy is cheap (shared-node edge + two small index
// vectors) and refcount-correct.
// ---------------------------------------------------------------------------
struct DD {
    dd::TDD tdd;

    // A default DD is the zero vector/matrix (terminal zero edge). Its
    // destructor safely decrements a terminal, which is a no-op for the node.
    DD() { tdd.e = dd::Edge<dd::mNode>::zero; }

    // Adopt a freshly built TDD. LimTDD's xarray_2_edge / cont / T_add2 return
    // edges whose node refcount is 0 — the caller must incRef to own a
    // reference (this is why test_fidelity incRefs every cont result). The DD
    // then owns exactly one reference, released by the destructor.
    explicit DD(dd::TDD&& t) : tdd(std::move(t)) { backendPackage().incRef(tdd.e); }

    DD(const DD& o) : tdd(o.tdd) { backendPackage().incRef(tdd.e); }
    DD(DD&& o) noexcept : tdd(std::move(o.tdd)) {
        o.tdd.e = dd::Edge<dd::mNode>::zero;  // leave source safe to destroy
    }

    DD& operator=(const DD& o) {
        if (this != &o) {
            backendPackage().decRef(tdd.e);
            tdd = o.tdd;
            backendPackage().incRef(tdd.e);
        }
        return *this;
    }
    DD& operator=(DD&& o) noexcept {
        if (this != &o) {
            backendPackage().decRef(tdd.e);
            tdd = std::move(o.tdd);
            o.tdd.e = dd::Edge<dd::mNode>::zero;
        }
        return *this;
    }

    ~DD() { backendPackage().decRef(tdd.e); }

    bool operator==(const DD& o) const {
        return tdd.e == o.tdd.e && sameIndexSet(tdd.index_set, o.tdd.index_set);
    }
    bool operator!=(const DD& o) const { return !(*this == o); }

  private:
    static bool sameIndexSet(const std::vector<dd::Index>& a,
                             const std::vector<dd::Index>& b) {
        if (a.size() != b.size()) {
            return false;
        }
        for (std::size_t i = 0; i < a.size(); ++i) {
            if (a[i].key != b[i].key || a[i].idx != b[i].idx) {
                return false;
            }
        }
        return true;
    }
};

// The shared package instance. DDVector::Initialize() and
// DDMatrix::Initialize() funnel through this; a function-local static makes
// repeated Initialize() calls naturally idempotent.
inline dd::Package<>& backendPackage() {
    static dd::Package<> pkg;  // default 300 qubits
    return pkg;
}

// Idempotent initialization. Called once at startup and again in multiple
// QReach constructors; must be safe to call repeatedly. The package is lazily
// constructed on first use and never needs a reset between runs.
//
// We also pre-register the qubit key names in a fixed INTERLEAVED order
// ("o0", "q0", "o1", "q1", …) so a matrix's row/col variable convention is
// deterministic regardless of build order. The interleaved row/col order is
// required for Package::cont to contract matrices correctly (the grouped
// "o0..oN, q0..qN" order produces a factor-of-2 error in cont2's key mapping).
inline void Initialize() {
    auto& pkg = backendPackage();
    constexpr unsigned int kMaxKey = 256;
    for (unsigned int q = 0; q < kMaxKey; ++q) {
        (void)pkg.varOrder.try_emplace("o" + std::to_string(q), static_cast<int>(2 * q));
        (void)pkg.varOrder.try_emplace("q" + std::to_string(q), static_cast<int>(2 * q + 1));
    }
}

// Conversion helpers between the public scalar and the package's ref-counted
// Complex (edge weights) / value ComplexValue.
inline dd::Complex toComplex(const DDComplex& c) {
    return backendPackage().cn.lookup(c.real(), c.imag());
}
inline DDComplex fromComplex(const dd::Complex& c) {
    return DDComplex(dd::ComplexTable<>::Entry::val(c.r),
                     dd::ComplexTable<>::Entry::val(c.i));
}
inline DDComplex fromComplexValue(const dd::ComplexValue& c) {
    return DDComplex(c.r, c.i);
}

// Scalar * DD — scale every amplitude by a scalar (contract §1).
//
// Refcount accounting (weights are canonical ComplexTable entries that start
// at refcount 0 and must be incRef'd to be owned):
//   - `DD r = d` copies and incRefs both the node and the old root weight;
//   - we release r's reference to the old weight, look up a fresh canonical
//     weight for the scaled value, and incRef it so r owns it.
inline DD operator*(const DDComplex& s, const DD& d) {
    DD r = d;  // shares node + old weight (both refcounts incremented)
    const auto cur = fromComplex(d.tdd.e.w);
    const double re = cur.real() * s.real() - cur.imag() * s.imag();
    const double im = cur.real() * s.imag() + cur.imag() * s.real();
    auto& cn = backendPackage().cn;
    cn.decRef(r.tdd.e.w);           // release the copy's ref to the old weight
    r.tdd.e.w = cn.lookup(re, im);  // fresh canonical weight (refcount 0)
    cn.incRef(r.tdd.e.w);           // r owns it (refcount 1)
    return r;
}

// DD + DD — pointwise (element-wise) addition (contract §1, used by
// InitializeWithAmplitudes / span). Operands must share the same index set;
// the result keeps a's index metadata.
inline DD operator+(const DD& a, const DD& b) {
    dd::TDD res;
    res.index_set = a.tdd.index_set;
    res.key_2_index = a.tdd.key_2_index;
    res.e = backendPackage().backendAdd(a.tdd.e, b.tdd.e);  // fresh edge, refcount 1
    return DD(std::move(res));
}

}  // namespace limtdd
