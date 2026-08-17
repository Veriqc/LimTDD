#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "dd/backend/DDMatrix.hpp"
#include "dd/backend/DDVector.hpp"

namespace {

int failures = 0;

void check(bool cond, const std::string& msg) {
    if (cond) {
        std::printf("  PASS  %s\n", msg.c_str());
    } else {
        std::printf("  FAIL  %s\n", msg.c_str());
        ++failures;
    }
}

// Check a 2-qubit vector has exactly the expected (index, re, im) entries.
void check2Qubit(const DDVector::DD& dd,
                 const std::vector<std::pair<unsigned int, std::pair<double, double>>>& expected,
                 const std::string& msg) {
    const auto amps = DDVector::GetNonZeroAmplitudes(dd);
    if (amps.size() != expected.size()) {
        std::printf("    got %zu entries (want %zu): ", amps.size(), expected.size());
        for (const auto& [i, a] : amps) {
            std::printf("[%u -> %g%+gi] ", i, a.real(), a.imag());
        }
        std::printf("\n");
        check(false, msg);
        return;
    }
    bool ok = true;
    for (std::size_t k = 0; k < expected.size(); ++k) {
        const auto& [wantIdx, wantVal] = expected[k];
        const auto& [gotIdx, gotVal] = amps[k];
        if (gotIdx != wantIdx ||
            std::abs(gotVal.real() - wantVal.first) > 1e-9 ||
            std::abs(gotVal.imag() - wantVal.second) > 1e-9) {
            ok = false;
            std::printf("    entry %zu: got [%u -> %g%+gi], want [%u -> %g%+gi]\n",
                        k, gotIdx, gotVal.real(), gotVal.imag(),
                        wantIdx, wantVal.first, wantVal.second);
        }
    }
    check(ok, msg);
}

// Check a vector has exactly the expected (index, re, im) non-zero entries.
void checkAmps(const DDVector::DD& dd,
               const std::vector<std::pair<unsigned int, std::pair<double, double>>>& expected,
               const std::string& msg) {
    const auto amps = DDVector::GetNonZeroAmplitudes(dd);
    if (amps.size() != expected.size()) {
        std::printf("    got %zu entries (want %zu): ", amps.size(), expected.size());
        for (const auto& [i, a] : amps) {
            std::printf("[%u -> %g%+gi] ", i, a.real(), a.imag());
        }
        std::printf("\n");
        check(false, msg);
        return;
    }
    bool ok = true;
    for (std::size_t k = 0; k < expected.size(); ++k) {
        const auto& [wantIdx, wantVal] = expected[k];
        const auto& [gotIdx, gotVal] = amps[k];
        if (gotIdx != wantIdx ||
            std::abs(gotVal.real() - wantVal.first) > 1e-9 ||
            std::abs(gotVal.imag() - wantVal.second) > 1e-9) {
            ok = false;
            std::printf("    entry %zu: got [%u -> %g%+gi], want [%u -> %g%+gi]\n",
                        k, gotIdx, gotVal.real(), gotVal.imag(),
                        wantIdx, wantVal.first, wantVal.second);
        }
    }
    check(ok, msg);
}

}  // namespace

int main() {
    DDVector::Initialize();
    DDMatrix::Initialize();

    const double s = 1.0 / std::sqrt(2.0);
    const auto zero2 = DDVector::MkBasisVector(2, 0);   // |00>
    const auto ten2 = DDVector::MkBasisVector(2, 2);    // |10>

    std::printf("== single-qubit gate on n (MkSingleQubitGateOnN) ==\n");
    {
        auto gate = DDMatrix::MkSingleQubitGateOnN(2, 0, DDMatrix::MkWalsh);
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(gate, zero2),
                    {{0, {s, 0}}, {2, {s, 0}}}, "H on qubit 0 of |00>");
    }
    {
        auto gate = DDMatrix::MkSingleQubitGateOnN(2, 1, DDMatrix::MkWalsh);
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(gate, zero2),
                    {{0, {s, 0}}, {1, {s, 0}}}, "H on qubit 1 of |00>");
    }
    {
        auto gate = DDMatrix::MkSingleQubitGateOnN(2, 0, DDMatrix::MkNegation);
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(gate, zero2),
                    {{2, {1.0, 0}}}, "X on qubit 0 of |00> = |10>");
    }
    {
        auto gate = DDMatrix::MkSingleQubitGateOnN(2, 1, DDMatrix::MkNegation);
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(gate, zero2),
                    {{1, {1.0, 0}}}, "X on qubit 1 of |00> = |01>");
    }
    {
        auto gate = DDMatrix::MkSingleQubitGateOnN(2, 0, DDMatrix::MkPauliZ);
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(gate, ten2),
                    {{2, {-1.0, 0}}}, "Z on qubit 0 of |10> = -|10>");
    }

    std::printf("== sequential application (renames rows -> state keys) ==\n");
    {
        auto h0 = DDMatrix::MkSingleQubitGateOnN(2, 0, DDMatrix::MkWalsh);
        auto z0 = DDMatrix::MkSingleQubitGateOnN(2, 0, DDMatrix::MkPauliZ);
        auto afterH = DDMatrix::MatrixMultiplyWithVector(h0, zero2);
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(z0, afterH),
                    {{0, {s, 0}}, {2, {-s, 0}}}, "Z·H on qubit 0 of |00>");
    }

    std::printf("== parameterized single-qubit gate on n ==\n");
    {
        // S = diag(1, i) on qubit 0; S|10> = i|10>
        auto s0 = DDMatrix::MkSingleQubitGateOnNWithParam(2, 0, DDMatrix::MkPhaseShift, 0.5);
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(s0, ten2),
                    {{2, {0.0, 1.0}}}, "S on qubit 0 of |10> = i|10>");
    }
    {
        // U3(pi/2, 0, 0) on qubit 0 ~ H up to global phase applied to |00>
        std::vector<double> u3p = {1.0, 0.0, 0.0};  // theta=pi (units of pi), phi=0, lambda=0
        auto g = DDMatrix::MkSingleQubitGateOnNWithParamVec(2, 0, DDMatrix::MkU3, u3p);
        // U3(theta=1 (pi), 0, 0) = [[cos(pi/2), -sin(pi/2)],[sin(pi/2), cos(pi/2)]] = X
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(g, zero2),
                    {{2, {1.0, 0}}}, "U3(pi,0,0) on qubit 0 of |00> = |10>");
    }

    std::printf("== KroneckerProduct builds the lifted gate ==\n");
    {
        // H on qubit 0, I on qubit 1 == MkSingleQubitGateOnN(2,0,MkWalsh)
        auto h0 = DDMatrix::KroneckerProduct(DDMatrix::MkWalsh(1), DDMatrix::MkIdRelation(1));
        check2Qubit(DDMatrix::MatrixMultiplyWithVector(h0, zero2),
                    {{0, {s, 0}}, {2, {s, 0}}}, "H0 ⊗ I1 on |00>");
    }

    std::printf("== multi-qubit gates ==\n");
    {
        auto gate = DDMatrix::MkCNOT(2, 2, 0, 1);
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, ten2),
                  {{3, {1.0, 0}}}, "CNOT(0,1) on |10> = |11>");
    }
    {
        auto gate = DDMatrix::MkCNOT(2, 2, 0, 1);
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, zero2),
                  {{0, {1.0, 0}}}, "CNOT(0,1) on |00> = |00>");
    }
    {
        // Toffoli on 4 qubits: control 0,1 target 2
        auto gate = DDMatrix::MkCCNOT(3, 4, 0, 1, 2);
        auto vec = DDVector::MkBasisVector(4, 12);  // |1100>
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, vec),
                  {{14, {1.0, 0}}}, "CCNOT(0,1,2) on |1100> = |1110>");
    }
    {
        auto gate = DDMatrix::MkSwap(2, 0, 1);
        auto vec = DDVector::MkBasisVector(2, 1);  // |01>
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, vec),
                  {{2, {1.0, 0}}}, "SWAP(0,1) on |01> = |10>");
    }
    {
        auto gate = DDMatrix::MkCP(2, 0, 1, 1.0);  // theta = pi -> phase -1
        auto vec = DDVector::MkBasisVector(2, 3);  // |11>
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, vec),
                  {{3, {-1.0, 0}}}, "CP(0,1,pi) on |11> = -|11>");
    }
    {
        auto gate = DDMatrix::MkiSwap(2, 0, 1);
        auto vec = DDVector::MkBasisVector(2, 1);  // |01>
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, vec),
                  {{2, {0.0, 1.0}}}, "iSWAP(0,1) on |01> = i|10>");
    }

    std::printf("== large-n compact gates (no dense O(4^n)) ==\n");
    {
        // 16-qubit system: X on target 15 (LSB) of |0...0> -> |0...01> (index 1).
        auto gate = DDMatrix::MkSingleQubitGateOnN(16, 15, DDMatrix::MkNegation);
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, DDVector::MkBasisVector(16, 0)),
                  {{1, {1.0, 0}}}, "X on qubit 15 of 16-qubit |0...0>");
    }
    {
        // 16-qubit CNOT(0, 15): control = qubit 0 (MSB), target = qubit 15 (LSB).
        auto gate = DDMatrix::MkCNOT(5, 16, 0, 15);
        // |10...0> = index 2^15; CNOT flips qubit 15 -> |10...01> = 2^15 + 1.
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, DDVector::MkBasisVector(16, 32768)),
                  {{32769, {1.0, 0}}}, "CNOT(0,15) on |10...0> -> |10...01>");
    }
    {
        // 16-qubit CCNOT(0, 1, 15): controls 0,1 (MSBs), target 15 (LSB).
        auto gate = DDMatrix::MkCCNOT(5, 16, 0, 1, 15);
        // |11...0> = 2^15 + 2^14; Toffoli flips qubit 15 -> |11...01> = 2^15 + 2^14 + 1.
        checkAmps(DDMatrix::MatrixMultiplyWithVector(gate, DDVector::MkBasisVector(16, 49152)),
                  {{49153, {1.0, 0}}}, "CCNOT(0,1,15) on |11...0> -> |11...01>");
    }

    std::printf("== Conjugate / Transpose / MatrixMultiply ==\n");
    {
        // (|0> + i|1>)/sqrt(2) -> conjugate -> (|0> - i|1>)/sqrt(2)
        const double s = 1.0 / std::sqrt(2.0);
        std::vector<double> amps = {s, 0.0, 0.0, s};
        auto v = DDVector::InitializeWithAmplitudes(1, amps);
        checkAmps(DDMatrix::Conjugate(v),
                  {{0, {s, 0.0}}, {1, {0.0, -s}}}, "Conjugate of (|0>+i|1>)/sqrt(2)");
    }
    {
        // Conjugate(S) = diag(1, -i), applied to |1> (1 qubit) -> -i|1>
        auto sc = DDMatrix::Conjugate(DDMatrix::MkSGate(1));
        checkAmps(DDMatrix::MatrixMultiplyWithVector(sc, DDVector::MkBasisVector(1, 1)),
                  {{1, {0.0, -1.0}}}, "Conjugate(S) on |1> = -i|1>");
    }
    {
        // X · X = I on |0> (1 qubit)
        auto xx = DDMatrix::MatrixMultiply(DDMatrix::MkNegation(1), DDMatrix::MkNegation(1));
        checkAmps(DDMatrix::MatrixMultiplyWithVector(xx, DDVector::MkBasisVector(1, 0)),
                  {{0, {1.0, 0}}}, "X·X = I on |0>");
    }
    {
        // H · H = I on |0> (1 qubit)
        auto hh = DDMatrix::MatrixMultiply(DDMatrix::MkWalsh(1), DDMatrix::MkWalsh(1));
        checkAmps(DDMatrix::MatrixMultiplyWithVector(hh, DDVector::MkBasisVector(1, 0)),
                  {{0, {1.0, 0}}}, "H·H = I on |0>");
    }
    {
        // Transpose([[1,2],[3,4]]) = [[1,3],[2,4]], on |0> (1 qubit) -> |0> + 2|1>
        std::vector<double> params = {1, 0, 2, 0, 3, 0, 4, 0};
        auto mt = DDMatrix::Transpose(DDMatrix::MkArbitrary(1, params));
        checkAmps(DDMatrix::MatrixMultiplyWithVector(mt, DDVector::MkBasisVector(1, 0)),
                  {{0, {1.0, 0}}, {1, {2.0, 0}}}, "Transpose([[1,2],[3,4]]) on |0>");
    }

    std::printf("\n%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
