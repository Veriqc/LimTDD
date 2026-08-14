#pragma once

// QReach backend — DDMatrix (gate matrices + matrix algebra). Contract §3.
//
// Index convention (internal; must agree with DDVector):
//   - state vector qubit q : key "q<q>", idx 0;
//   - gate matrix: row (output) key "o<q>", column (input) key "q<q>", idx 0.
// `cont(gate, vec)` therefore contracts the gate's column with the vector and
// leaves the row. MatrixMultiplyWithVector renames the surviving row keys
// "o<q>" back to "q<q>" so the result is a well-formed state vector for the
// next application. (The rename is safe: cont renumbers output variables
// densely by varOrder rank, and "o<q>" / "q<q>" sort in the same qubit order.)
//
// KroneckerProduct(a, b) with disjoint index sets is exactly cont(a, b) —
// this is how TensorNetwork::cont builds tensor products of gate tensors.

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "dd/Tensor.hpp"
#include "dd/backend/DDTypes.hpp"
#include "dd/backend/DDVector.hpp"

namespace DDMatrix {

using limtdd::DD;
using limtdd::DDComplex;

inline void Initialize() { limtdd::Initialize(); }

// Forward declaration — used by detail::liftSingleQubitGate before its definition.
inline DD KroneckerProduct(DD a, DD b);

namespace detail {

// Angles are in units of π (contract §6.4).
inline double cospi(double x) { return std::cos(dd::PI * x); }
inline double sinpi(double x) { return std::sin(dd::PI * x); }

// Dense matrix construction is O(4^n) memory (O(8^n) for MatrixMultiply).
// Guard against overflow / OOM — small n only (contract §6.7).
constexpr unsigned int kMaxDenseQubits = 12;

// Build a 2x2 matrix TDD from row-major complex entries, with the given row /
// column index keys.
inline DD matrixFrom2x2(const std::array<std::complex<double>, 4>& m,
                        const std::string& rowKey, const std::string& colKey) {
    xt::xarray<dd::ComplexValue> arr(std::vector<std::size_t>{2, 2});
    for (std::size_t r = 0; r < 2; ++r) {
        for (std::size_t c = 0; c < 2; ++c) {
            arr(r, c) = dd::ComplexValue{m[r * 2 + c].real(), m[r * 2 + c].imag()};
        }
    }
    dd::Tensor tensor(arr, {{rowKey, 0}, {colKey, 0}}, "gate2x2");
    return DD(tensor.to_tdd(&limtdd::backendPackage()));
}

inline DD identity2x2(const std::string& rowKey, const std::string& colKey) {
    return matrixFrom2x2({{{1, 0}, {0, 0}, {0, 0}, {1, 0}}}, rowKey, colKey);
}

// Rename "o<q>" -> "q<q>" in a TDD's index metadata. The node variable numbers
// are left untouched (they are already dense and match the qubit order).
inline void renameRowsToState(dd::TDD& tdd) {
    for (auto& ix : tdd.index_set) {
        if (!ix.key.empty() && ix.key[0] == 'o') {
            ix.key[0] = 'q';
        }
    }
    for (auto& k : tdd.key_2_index) {
        if (!k.empty() && k[0] == 'o') {
            k[0] = 'q';
        }
    }
}

// Enumerate the four entries of a 2x2 matrix DD (row-major: [row*2+col]).
// Because the backend pre-registers "o<q>" before "q<q>", the row key has a
// lower varOrder rank than the column key, so the matrix's root variable is the
// COLUMN and its child variable is the ROW.
inline std::array<std::complex<double>, 4> extract2x2(const DD& m) {
    std::array<std::complex<double>, 4> out{};
    const auto& e = m.tdd.e;
    if (e.p == nullptr || e.p->v == -1) {
        return out;
    }
    const int colVar = static_cast<int>(e.p->v);
    for (int c = 0; c < 2; ++c) {
        const auto colEdge = limtdd::backendPackage().backendSlice(e, colVar, c);
        if (colEdge.p == nullptr || colEdge.p->v == -1) {
            continue;
        }
        const int rowVar = static_cast<int>(colEdge.p->v);
        for (int r = 0; r < 2; ++r) {
            const auto leaf = limtdd::backendPackage().backendSlice(colEdge, rowVar, r);
            const auto w = limtdd::fromComplex(leaf.w);
            out[r * 2 + c] = std::complex<double>{w.real(), w.imag()};
        }
    }
    return out;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Single-qubit gate matrices (2x2, built on qubit 0).
// ---------------------------------------------------------------------------
inline DD MkIdRelation(unsigned int /*level*/) {
    return detail::identity2x2("o0", "q0");
}
inline DD MkWalsh(unsigned int /*level*/) {
    const double s = 1.0 / std::sqrt(2.0);
    return detail::matrixFrom2x2({{{s, 0}, {s, 0}, {s, 0}, {-s, 0}}}, "o0", "q0");
}
inline DD MkNegation(unsigned int /*level*/) {
    return detail::matrixFrom2x2({{{0, 0}, {1, 0}, {1, 0}, {0, 0}}}, "o0", "q0");
}
inline DD MkPauliY(unsigned int /*level*/) {
    return detail::matrixFrom2x2({{{0, 0}, {0, -1}, {0, 1}, {0, 0}}}, "o0", "q0");
}
inline DD MkPauliZ(unsigned int /*level*/) {
    return detail::matrixFrom2x2({{{1, 0}, {0, 0}, {0, 0}, {-1, 0}}}, "o0", "q0");
}
inline DD MkSGate(unsigned int /*level*/) {
    return detail::matrixFrom2x2({{{1, 0}, {0, 0}, {0, 0}, {0, 1}}}, "o0", "q0");
}
inline DD MkPhaseShift(unsigned int /*level*/, double theta) {
    const double re = detail::cospi(theta);
    const double im = detail::sinpi(theta);
    return detail::matrixFrom2x2({{{1, 0}, {0, 0}, {0, 0}, {re, im}}}, "o0", "q0");
}
inline DD MkU3(unsigned int /*level*/, std::vector<double> params) {
    const double th = params[0], ph = params[1], la = params[2];
    const double c = detail::cospi(th / 2.0);
    const double s = detail::sinpi(th / 2.0);
    return detail::matrixFrom2x2(
        {{{c, 0},
          {-s * detail::cospi(la), -s * detail::sinpi(la)},
          {s * detail::cospi(ph), s * detail::sinpi(ph)},
          {c * detail::cospi(ph + la), c * detail::sinpi(ph + la)}}},
        "o0", "q0");
}
inline DD MkArbitrary(unsigned int /*level*/, std::vector<double> params) {
    return detail::matrixFrom2x2(
        {{{params[0], params[1]},
          {params[2], params[3]},
          {params[4], params[5]},
          {params[6], params[7]}}},
        "o0", "q0");
}

// ---------------------------------------------------------------------------
// Single-qubit gate lifted onto an n-qubit system.
// ---------------------------------------------------------------------------
namespace detail {
// Build I ⊗ … ⊗ G(target) ⊗ … ⊗ I as a dense 2n-dim matrix, then convert to a
// DD via Tensor::to_tdd. This avoids composing gates through cont (which
// produces a malformed DD on re-contraction — see the row/col variable bug).
// NOTE: O(4^n) dense — only suitable for small n. A compact DD-level lifting
// is a future optimization.
inline DD denseSingleQubitGateOnN(unsigned int n, unsigned int target,
                                  const std::array<std::complex<double>, 4>& G) {
    if (target >= n) {
        throw std::invalid_argument("denseSingleQubitGateOnN: target out of range");
    }
    if (n > kMaxDenseQubits) {
        throw std::invalid_argument("denseSingleQubitGateOnN: too many qubits (dense O(4^n))");
    }
    std::vector<std::size_t> shape(2 * n, 2);
    xt::xarray<dd::ComplexValue> arr(shape);
    arr.fill(dd::ComplexValue{0.0, 0.0});
    const std::size_t dim = std::size_t{1} << n;
    for (std::size_t r = 0; r < dim; ++r) {
        for (std::size_t c = 0; c < dim; ++c) {
            bool match = true;
            for (unsigned int q = 0; q < n; ++q) {
                if (q == target) {
                    continue;
                }
                if (((r >> (n - 1 - q)) & 1u) != ((c >> (n - 1 - q)) & 1u)) {
                    match = false;
                    break;
                }
            }
            if (!match) {
                continue;
            }
            const unsigned int rt = (r >> (n - 1 - target)) & 1u;
            const unsigned int ct = (c >> (n - 1 - target)) & 1u;
            arr.data()[r * dim + c] =
                dd::ComplexValue{G[rt * 2 + ct].real(), G[rt * 2 + ct].imag()};
        }
    }
    std::vector<dd::Index> iset;
    iset.reserve(2 * n);
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"o" + std::to_string(q), 0});
    }
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"q" + std::to_string(q), 0});
    }
    dd::Tensor tensor(arr, iset, "lifted");
    return DD(tensor.to_tdd(&limtdd::backendPackage()));
}
}  // namespace detail

inline DD MkSingleQubitGateOnN(unsigned int n, unsigned int target,
                               DD (*gate1q)(unsigned int)) {
    return detail::denseSingleQubitGateOnN(n, target, detail::extract2x2(gate1q(1)));
}

inline DD MkSingleQubitGateOnNWithParam(unsigned int n, unsigned int target,
                                        DD (*gate1q)(unsigned int, double),
                                        double theta) {
    return detail::denseSingleQubitGateOnN(n, target, detail::extract2x2(gate1q(1, theta)));
}

inline DD MkSingleQubitGateOnNWithParamVec(unsigned int n, unsigned int target,
                                           DD (*gate1q)(unsigned int, std::vector<double>),
                                           std::vector<double> v) {
    return detail::denseSingleQubitGateOnN(n, target, detail::extract2x2(gate1q(1, v)));
}

// ---------------------------------------------------------------------------
// Multi-qubit gate matrices (n = number of qubits).
// ---------------------------------------------------------------------------
namespace detail {
// Build a dense n-qubit gate matrix from a per-(row, col) amplitude function.
// O(4^n) — small n only.
template <typename F>
inline DD denseGateFromFunction(unsigned int n, F&& amplitude) {
    if (n > kMaxDenseQubits) {
        throw std::invalid_argument("denseGateFromFunction: too many qubits (dense O(4^n))");
    }
    std::vector<std::size_t> shape(2 * n, 2);
    xt::xarray<dd::ComplexValue> arr(shape);
    arr.fill(dd::ComplexValue{0.0, 0.0});
    const std::size_t dim = std::size_t{1} << n;
    for (std::size_t r = 0; r < dim; ++r) {
        for (std::size_t c = 0; c < dim; ++c) {
            const auto a = amplitude(r, c);
            if (a.real() == 0.0 && a.imag() == 0.0) {
                continue;
            }
            arr.data()[r * dim + c] = dd::ComplexValue{a.real(), a.imag()};
        }
    }
    std::vector<dd::Index> iset;
    iset.reserve(2 * n);
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"o" + std::to_string(q), 0});
    }
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"q" + std::to_string(q), 0});
    }
    dd::Tensor tensor(arr, iset, "gate");
    return DD(tensor.to_tdd(&limtdd::backendPackage()));
}
}  // namespace detail

inline DD MkCNOT(unsigned int /*level*/, unsigned int n, long ctrl, long tgt) {
    return detail::denseGateFromFunction(n, [&](std::size_t r, std::size_t c) {
        std::size_t flipped = c;
        if (((c >> (n - 1 - ctrl)) & 1u) != 0u) {
            flipped ^= (std::size_t{1} << (n - 1 - tgt));
        }
        return (r == flipped) ? std::complex<double>{1.0, 0.0}
                              : std::complex<double>{0.0, 0.0};
    });
}

inline DD MkCCNOT(unsigned int /*level*/, unsigned int n, long c1, long c2, long tgt) {
    return detail::denseGateFromFunction(n, [&](std::size_t r, std::size_t c) {
        std::size_t flipped = c;
        if (((c >> (n - 1 - c1)) & 1u) != 0u && ((c >> (n - 1 - c2)) & 1u) != 0u) {
            flipped ^= (std::size_t{1} << (n - 1 - tgt));
        }
        return (r == flipped) ? std::complex<double>{1.0, 0.0}
                              : std::complex<double>{0.0, 0.0};
    });
}

inline DD MkSwap(unsigned int level, long i, long j) {
    const unsigned int n = 1u << (level - 1);
    return detail::denseGateFromFunction(n, [&](std::size_t r, std::size_t c) {
        const unsigned int bi = (c >> (n - 1 - i)) & 1u;
        const unsigned int bj = (c >> (n - 1 - j)) & 1u;
        std::size_t swapped = c;
        swapped &= ~((std::size_t{1} << (n - 1 - i)) | (std::size_t{1} << (n - 1 - j)));
        swapped |= (static_cast<std::size_t>(bi) << (n - 1 - j)) |
                   (static_cast<std::size_t>(bj) << (n - 1 - i));
        return (r == swapped) ? std::complex<double>{1.0, 0.0}
                              : std::complex<double>{0.0, 0.0};
    });
}

inline DD MkiSwap(unsigned int level, long i, long j) {
    const unsigned int n = 1u << (level - 1);
    return detail::denseGateFromFunction(n, [&](std::size_t r, std::size_t c) {
        const unsigned int bi = (c >> (n - 1 - i)) & 1u;
        const unsigned int bj = (c >> (n - 1 - j)) & 1u;
        if (bi == bj) {
            return (r == c) ? std::complex<double>{1.0, 0.0}
                            : std::complex<double>{0.0, 0.0};
        }
        const std::size_t swapped =
            c ^ ((std::size_t{1} << (n - 1 - i)) | (std::size_t{1} << (n - 1 - j)));
        return (r == swapped) ? std::complex<double>{0.0, 1.0}
                              : std::complex<double>{0.0, 0.0};
    });
}

inline DD MkCP(unsigned int level, long ctrl, long tgt, double theta) {
    const unsigned int n = 1u << (level - 1);
    const std::complex<double> phase{detail::cospi(theta), detail::sinpi(theta)};
    return detail::denseGateFromFunction(n, [&](std::size_t r, std::size_t c) {
        if (r != c) {
            return std::complex<double>{0.0, 0.0};
        }
        if (((c >> (n - 1 - ctrl)) & 1u) != 0u && ((c >> (n - 1 - tgt)) & 1u) != 0u) {
            return phase;
        }
        return std::complex<double>{1.0, 0.0};
    });
}

// ---------------------------------------------------------------------------
// Matrix algebra.
// ---------------------------------------------------------------------------
namespace detail {
// Enumerate an n-qubit matrix DD (row keys "o<q>", col keys "q<q>") into a
// dense row-major vector of 4^n complex entries. With the pre-registered
// variable order, column variables are 0..n-1 and row variables are n..2n-1.
// NOTE: O(4^n) — small n only.
inline std::vector<std::complex<double>> matrixToDense(const DD& m, unsigned int n) {
    if (n > kMaxDenseQubits) {
        throw std::invalid_argument("matrixToDense: too many qubits (dense O(4^n))");
    }
    std::vector<std::complex<double>> out(std::size_t{1} << (2 * n),
                                          std::complex<double>{0.0, 0.0});
    struct Walker {
        unsigned int n;
        std::vector<std::complex<double>>& out;
        void go(const dd::Edge<dd::mNode>& e, unsigned long long row, unsigned long long col) {
            const auto w = limtdd::fromComplex(e.w);
            if (w.real() == 0.0 && w.imag() == 0.0) {
                return;
            }
            if (e.p == nullptr || e.p->v == -1) {
                out[row * (1ull << n) + col] = {w.real(), w.imag()};
                return;
            }
            const int v = static_cast<int>(e.p->v);
            for (int c = 0; c < 2; ++c) {
                const auto ch = limtdd::backendPackage().backendSlice(e, v, c);
                // Interleaved order: row qubit q -> var 2q (even),
                // column qubit q -> var 2q+1 (odd).
                if (v % 2 == 0) {
                    const int q = v / 2;
                    go(ch, row | (static_cast<unsigned long long>(c) << (n - 1 - q)), col);
                } else {
                    const int q = (v - 1) / 2;
                    go(ch, row, col | (static_cast<unsigned long long>(c) << (n - 1 - q)));
                }
            }
        }
    } walker{n, out};
    walker.go(m.tdd.e, 0, 0);
    return out;
}

inline unsigned int matrixQubitCount(const DD& m) {
    unsigned int n = 0;
    for (const auto& ix : m.tdd.index_set) {
        if (!ix.key.empty() && ix.key[0] == 'o') {
            ++n;
        }
    }
    return n;
}

inline bool isMatrix(const DD& m) {
    for (const auto& ix : m.tdd.index_set) {
        if (!ix.key.empty() && ix.key[0] == 'o') {
            return true;
        }
    }
    return false;
}

// Build an n-qubit matrix DD from a dense row-major vector (2^n x 2^n entries).
inline DD matrixFromDense(const std::vector<std::complex<double>>& dense, unsigned int n) {
    const std::size_t dim = std::size_t{1} << n;
    std::vector<std::size_t> shape(2 * n, 2);
    xt::xarray<dd::ComplexValue> arr(shape);
    for (std::size_t r = 0; r < dim; ++r) {
        for (std::size_t c = 0; c < dim; ++c) {
            const auto& v = dense[r * dim + c];
            arr.data()[r * dim + c] = dd::ComplexValue{v.real(), v.imag()};
        }
    }
    std::vector<dd::Index> iset;
    iset.reserve(2 * n);
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"o" + std::to_string(q), 0});
    }
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"q" + std::to_string(q), 0});
    }
    dd::Tensor tensor(arr, iset, "matrix");
    return DD(tensor.to_tdd(&limtdd::backendPackage()));
}
}  // namespace detail

inline DD KroneckerProduct(DD a, DD b) {
    const unsigned int na = detail::matrixQubitCount(a);
    const unsigned int nb = detail::matrixQubitCount(b);
    const unsigned int n = na + nb;
    const auto da = detail::matrixToDense(a, na);
    const auto db = detail::matrixToDense(b, nb);
    const std::size_t dim = std::size_t{1} << n;
    const std::size_t dima = std::size_t{1} << na;
    const std::size_t dimb = std::size_t{1} << nb;

    std::vector<std::size_t> shape(2 * n, 2);
    xt::xarray<dd::ComplexValue> arr(shape);
    arr.fill(dd::ComplexValue{0.0, 0.0});
    for (std::size_t ra = 0; ra < dima; ++ra) {
        for (std::size_t ca = 0; ca < dima; ++ca) {
            const auto va = da[ra * dima + ca];
            if (va.real() == 0.0 && va.imag() == 0.0) {
                continue;
            }
            for (std::size_t rb = 0; rb < dimb; ++rb) {
                for (std::size_t cb = 0; cb < dimb; ++cb) {
                    const auto vb = db[rb * dimb + cb];
                    if (vb.real() == 0.0 && vb.imag() == 0.0) {
                        continue;
                    }
                    const std::size_t r = (ra << nb) | rb;
                    const std::size_t c = (ca << nb) | cb;
                    const auto p = va * vb;
                    arr.data()[r * dim + c] = dd::ComplexValue{p.real(), p.imag()};
                }
            }
        }
    }
    std::vector<dd::Index> iset;
    iset.reserve(2 * n);
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"o" + std::to_string(q), 0});
    }
    for (unsigned int q = 0; q < n; ++q) {
        iset.push_back({"q" + std::to_string(q), 0});
    }
    dd::Tensor tensor(arr, iset, "kron");
    return DD(tensor.to_tdd(&limtdd::backendPackage()));
}

// gate · vec — the hot path. Contracts the gate's column with the vector and
// renames the surviving row keys back to the state convention.
inline DD MatrixMultiplyWithVector(DD gate, DD vec) {
    dd::TDD res = limtdd::backendPackage().cont(gate.tdd, vec.tdd);
    detail::renameRowsToState(res);
    return DD(std::move(res));
}

// Complex conjugate of a DD (matrix or vector). Implemented via dense
// enumeration (weights and map phases are unfolded by the walk, conjugated,
// and re-derived on rebuild). O(2^n) — a compact DD-level conjugate (negating
// each edge weight and each map's rotate/extra_phase) is a future optimization.
inline DD Conjugate(DD c) {
    if (detail::isMatrix(c)) {
        const unsigned int n = detail::matrixQubitCount(c);
        auto dense = detail::matrixToDense(c, n);
        for (auto& v : dense) {
            v = std::conj(v);
        }
        return detail::matrixFromDense(dense, n);
    }
    const unsigned int n = static_cast<unsigned int>(c.tdd.index_set.size());
    if (n == 0) {
        return c;
    }
    const std::size_t dim = std::size_t{1} << n;
    std::vector<double> amps(2 * dim, 0.0);
    for (const auto& [idx, a] : DDVector::GetNonZeroAmplitudes(c, 0.0)) {
        amps[idx] = a.real();
        amps[dim + idx] = -a.imag();  // conjugate
    }
    return DDVector::InitializeWithAmplitudes(n, amps);
}

// Transpose. For a matrix this swaps row and column; for a vector it is the
// identity (a vector has no second index to swap). NOTE: the exact resetall
// semantics (Conjugate ∘ Transpose on a vector) need a dense cross-check.
inline DD Transpose(DD c) {
    if (!detail::isMatrix(c)) {
        return c;
    }
    const unsigned int n = detail::matrixQubitCount(c);
    auto dense = detail::matrixToDense(c, n);
    const std::size_t dim = std::size_t{1} << n;
    std::vector<std::complex<double>> t(dense.size(), std::complex<double>{0.0, 0.0});
    for (std::size_t r = 0; r < dim; ++r) {
        for (std::size_t c = 0; c < dim; ++c) {
            t[c * dim + r] = dense[r * dim + c];
        }
    }
    return detail::matrixFromDense(t, n);
}

// Matrix · matrix (both square n-qubit matrices). O(8^n) dense — small n only.
inline DD MatrixMultiply(DD a, DD b) {
    const unsigned int na = detail::matrixQubitCount(a);
    const unsigned int nb = detail::matrixQubitCount(b);
    if (na != nb) {
        throw std::invalid_argument("MatrixMultiply: size mismatch");
    }
    const unsigned int n = na;
    auto da = detail::matrixToDense(a, n);
    auto db = detail::matrixToDense(b, n);
    const std::size_t dim = std::size_t{1} << n;
    std::vector<std::complex<double>> prod(dim * dim, std::complex<double>{0.0, 0.0});
    for (std::size_t r = 0; r < dim; ++r) {
        for (std::size_t c = 0; c < dim; ++c) {
            std::complex<double> acc{0.0, 0.0};
            for (std::size_t k = 0; k < dim; ++k) {
                acc += da[r * dim + k] * db[k * dim + c];
            }
            prod[r * dim + c] = acc;
        }
    }
    return detail::matrixFromDense(prod, n);
}

}  // namespace DDMatrix
