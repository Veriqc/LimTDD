#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

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

void expectSingle(const DDVector::DD& dd, unsigned int index, double re, double im,
                  const std::string& msg) {
    const auto amps = DDVector::GetNonZeroAmplitudes(dd);
    bool ok = amps.size() == 1 && amps[0].first == index &&
              std::abs(amps[0].second.real() - re) < 1e-9 &&
              std::abs(amps[0].second.imag() - im) < 1e-9;
    if (!ok) {
        std::printf("    got %zu entries: ", amps.size());
        for (const auto& [i, a] : amps) {
            std::printf("[%u -> %g%+gi] ", i, a.real(), a.imag());
        }
        std::printf("\n");
    }
    check(ok, msg);
}

}  // namespace

int main() {
    DDVector::Initialize();

    std::printf("== MkBasisVector (integer) ==\n");
    expectSingle(DDVector::MkBasisVector(1, 0), 0, 1.0, 0.0, "level1 idx0 = |00>");
    expectSingle(DDVector::MkBasisVector(1, 1), 1, 1.0, 0.0, "level1 idx1 = |01>");
    expectSingle(DDVector::MkBasisVector(1, 2), 2, 1.0, 0.0, "level1 idx2 = |10>");
    expectSingle(DDVector::MkBasisVector(1, 3), 3, 1.0, 0.0, "level1 idx3 = |11>");
    expectSingle(DDVector::MkBasisVector(0, 0), 0, 1.0, 0.0, "level0 idx0 = |0>");
    expectSingle(DDVector::MkBasisVector(0, 1), 1, 1.0, 0.0, "level0 idx1 = |1>");

    std::printf("== MkBasisVector (bitstring, big-endian) ==\n");
    expectSingle(DDVector::MkBasisVector(1, std::string("10")), 2, 1.0, 0.0,
                 "level1 \"10\" = |10> == idx2");
    expectSingle(DDVector::MkBasisVector(1, std::string("01")), 1, 1.0, 0.0,
                 "level1 \"01\" = |01> == idx1");
    expectSingle(DDVector::MkBasisVector(2, std::string("1010")), 10, 1.0, 0.0,
                 "level2 \"1010\" = idx10");

    std::printf("== GetLevel ==\n");
    check(DDVector::GetLevel(DDVector::MkBasisVector(0, 0)) == 0, "level0 vector -> level 0");
    check(DDVector::GetLevel(DDVector::MkBasisVector(1, 0)) == 1, "level1 vector -> level 1");
    check(DDVector::GetLevel(DDVector::MkBasisVector(2, 0)) == 2, "level2 vector -> level 2");

    std::printf("== NoDistinctionNode ==\n");
    {
        auto dd = DDVector::NoDistinctionNode(1, DDVector::DDComplex(1.0, 0.0));
        auto amps = DDVector::GetNonZeroAmplitudes(dd);
        check(amps.size() == 4, "level1 constant vector has 4 entries");
        bool allOne = true;
        for (const auto& [i, a] : amps) {
            (void)i;
            if (std::abs(a.real() - 1.0) > 1e-9 || std::abs(a.imag()) > 1e-9) {
                allOne = false;
            }
        }
        check(allOne, "all amplitudes == 1");
    }

    std::printf("== IsApproximatelyZero ==\n");
    check(!DDVector::IsApproximatelyZero(DDVector::MkBasisVector(1, 0)),
          "basis vector is not zero");
    check(!DDVector::IsApproximatelyZero(DDVector::NoDistinctionNode(1, DDVector::DDComplex(1.0))),
          "constant vector is not zero");
    {
        auto zero = DDVector::DDComplex(0.0, 0.0) * DDVector::MkBasisVector(1, 0);
        check(DDVector::IsApproximatelyZero(zero), "scaled-to-zero vector is zero");
        check(DDVector::GetNonZeroAmplitudes(zero).empty(),
              "GetNonZeroAmplitudes(zero) short-circuits (empty)");
        auto sa = DDVector::ExtractSingleAmplitude(zero);
        check(std::abs(sa.real()) < 1e-9 && std::abs(sa.imag()) < 1e-9,
              "ExtractSingleAmplitude(zero) = 0");
    }

    std::printf("== Normalize ==\n");
    {
        auto dd = DDVector::Normalize(DDVector::NoDistinctionNode(1, DDVector::DDComplex(1.0, 0.0)));
        auto amps = DDVector::GetNonZeroAmplitudes(dd);
        check(amps.size() == 4, "normalized constant has 4 entries");
        bool allHalf = true;
        for (const auto& [i, a] : amps) {
            (void)i;
            if (std::abs(a.real() - 0.5) > 1e-9 || std::abs(a.imag()) > 1e-9) {
                allHalf = false;
            }
        }
        check(allHalf, "each amplitude == 1/sqrt(4) == 0.5");
    }

    std::printf("== InitializeWithAmplitudes ==\n");
    {
        // |00> on 2 qubits: re = [1,0,0,0], im = [0,0,0,0]
        std::vector<double> amps = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        auto dd = DDVector::InitializeWithAmplitudes(2, amps);
        expectSingle(dd, 0, 1.0, 0.0, "qnum=2 amplitudes -> |00>");
    }
    {
        // (|00> + |11>)/sqrt(2) on 2 qubits
        const double s = 1.0 / std::sqrt(2.0);
        std::vector<double> amps = {s, 0.0, 0.0, s, 0.0, 0.0, 0.0, 0.0};
        auto dd = DDVector::InitializeWithAmplitudes(2, amps);
        auto got = DDVector::GetNonZeroAmplitudes(dd);
        check(got.size() == 2 && got[0].first == 0 && got[1].first == 3,
              "Bell-like amplitudes at idx 0 and 3");
    }

    std::printf("== VectorToMatrixInterleaved (no-op) ==\n");
    {
        auto v = DDVector::MkBasisVector(1, 2);
        auto m = DDVector::VectorToMatrixInterleaved(v);
        check(m == v, "no-op returns the same vector");
    }

    std::printf("== InnerProduct ==\n");
    {
        const double s = 1.0 / std::sqrt(2.0);
        auto z0 = DDVector::MkBasisVector(0, 0);   // |0>
        auto z1 = DDVector::MkBasisVector(0, 1);   // |1>
        auto plus = DDVector::InitializeWithAmplitudes(1, {s, s, 0.0, 0.0});    // (|0>+|1>)/sqrt2
        auto minus = DDVector::InitializeWithAmplitudes(1, {s, -s, 0.0, 0.0});  // (|0>-|1>)/sqrt2
        auto iphase = DDVector::InitializeWithAmplitudes(1, {s, 0.0, 0.0, s});  // (|0>+i|1>)/sqrt2

        check(std::abs(DDVector::InnerProduct(z0, z0).real() - 1.0) < 1e-9,
              "<0|0> = 1");
        check(std::abs(DDVector::InnerProduct(z0, z1).real()) < 1e-9 &&
              std::abs(DDVector::InnerProduct(z0, z1).imag()) < 1e-9,
              "<0|1> = 0");
        check(std::abs(DDVector::InnerProduct(plus, plus).real() - 1.0) < 1e-9,
              "<+|+> = 1");
        check(std::abs(DDVector::InnerProduct(plus, minus).real()) < 1e-9,
              "<+|-> = 0");
        check(std::abs(DDVector::InnerProduct(iphase, iphase).real() - 1.0) < 1e-9 &&
              std::abs(DDVector::InnerProduct(iphase, iphase).imag()) < 1e-9,
              "<(0+i1)/sqrt2 | itself> = 1 (phase cancels)");
    }
    {
        // 2-qubit orthogonal states: <(00+11)/sqrt2 | (00-11)/sqrt2> = 0
        const double s = 1.0 / std::sqrt(2.0);
        std::vector<double> a = {s, 0.0, 0.0, s, 0.0, 0.0, 0.0, 0.0};
        std::vector<double> b = {s, 0.0, 0.0, -s, 0.0, 0.0, 0.0, 0.0};
        auto va = DDVector::InitializeWithAmplitudes(2, a);
        auto vb = DDVector::InitializeWithAmplitudes(2, b);
        auto ip = DDVector::InnerProduct(va, vb);
        check(std::abs(ip.real()) < 1e-9 && std::abs(ip.imag()) < 1e-9,
              "orthogonal 2-qubit inner product = 0 (no throw)");
    }

    std::printf("== InnerProduct scale (regression: limtdd-innerproduct-bug-report) ==\n");
    {
        // <v|v> must ALWAYS be exactly 1, for any qubit count. The old
        // cont-based implementation returned 2^(n-1) for n >= 2 qubits.
        check(std::abs(DDVector::InnerProduct(DDVector::MkBasisVector(1, 0),
                                              DDVector::MkBasisVector(1, 0)).real() - 1.0) < 1e-9,
              "<00|00> (2q) == 1");
        check(std::abs(DDVector::InnerProduct(DDVector::MkBasisVector(2, 10),
                                              DDVector::MkBasisVector(2, 10)).real() - 1.0) < 1e-9,
              "<1010|1010> (4q) == 1");

        // Gram-Schmidt projection coefficients: v3 is a linear combination of
        // v1 and v2, so <v1|v3> etc. must match the exact dense values (the
        // cont-based version returned -3*2^3 * <v1|v3> and kept v3 as an
        // extra support vector).
        const double s2 = 1.0 / std::sqrt(2.0);
        const double s10 = 1.0 / std::sqrt(10.0);
        auto mk = [](const std::vector<double>& amps) {
            return DDVector::InitializeWithAmplitudes(4, amps);
        };
        std::vector<double> a1(32, 0.0), a2(32, 0.0), a3(32, 0.0);
        a1[10] = s2; a1[16 + 10] = s2;                 // v1 = e^{i pi/4}|1010>
        a2[16 + 10] = s2; a2[16 + 14] = -s2;           // v2 = (i/sqrt2)(|1010>-|1110>)
        a3[16 + 10] = s10; a3[16 + 14] = -3.0 * s10;   // v3 = (i/sqrt10)|1010>-(3i/sqrt10)|1110>
        auto v1 = mk(a1), v2 = mk(a2), v3 = mk(a3);

        auto near = [](DDVector::DDComplex c, double re, double im) {
            return std::abs(c.real() - re) < 1e-6 && std::abs(c.imag() - im) < 1e-6;
        };
        check(near(DDVector::InnerProduct(v1, v2), 0.5, 0.5), "<v1|v2> = 0.5+0.5i");
        check(near(DDVector::InnerProduct(v1, v3), 0.22360679774997896, 0.22360679774997896),
              "<v1|v3> = (1+i)/sqrt(20)");
        check(near(DDVector::InnerProduct(v2, v3), 2.0 / std::sqrt(5.0), 0.0),
              "<v2|v3> = 2/sqrt(5)");
    }

    std::printf("\n%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
