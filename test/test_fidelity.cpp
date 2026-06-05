#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
#include "dd/Tensor.hpp"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace dd;

namespace {

struct InjectedPauliError {
    std::size_t afterGate = 0;
    std::size_t qubit = 0;
    char pauli = 'X';
};

struct ContractionStats {
    TDD tdd{};
    unsigned int maxNode = 0;
    double timeSeconds = 0.0;
};

xt::xarray<dd::ComplexValue> stateToArray(const BasisStates& state) {
    switch (state) {
        case BasisStates::zero:
            return {complex_one, complex_zero};
        case BasisStates::one:
            return {complex_zero, complex_one};
        case BasisStates::plus:
            return {complex_SQRT2_2, complex_SQRT2_2};
        case BasisStates::minus:
            return {complex_SQRT2_2, complex_mSQRT2_2};
        case BasisStates::right:
            return {complex_SQRT2_2, complex_iSQRT2_2};
        case BasisStates::left:
            return {complex_SQRT2_2, complex_miSQRT2_2};
    }
    throw std::invalid_argument("Unsupported basis state");
}

std::string readFile(const std::string& filename) {
    std::ifstream fileStream(filename);
    if (!fileStream.is_open()) {
        throw std::runtime_error("Unable to open QASM file: " + filename);
    }

    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    return buffer.str();
}

std::vector<InjectedPauliError> sampleInjectedErrors(
    const qc::QuantumComputation& circuit,
    const std::size_t errorCount,
    const std::uint64_t seed) {
    if (errorCount == 0) {
        return {};
    }
    if (circuit.getNops() == 0) {
        throw std::invalid_argument("Cannot inject errors into an empty circuit");
    }

    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<std::size_t> gateDistribution(0, circuit.getNops() - 1);
    std::uniform_int_distribution<std::size_t> qubitDistribution(0, circuit.getNqubits() - 1);
    std::uniform_int_distribution<int> pauliDistribution(0, 2);

    std::vector<InjectedPauliError> errors;
    errors.reserve(errorCount);

    for (std::size_t index = 0; index < errorCount; ++index) {
        const int pauliIndex = pauliDistribution(generator);
        const char pauli = pauliIndex == 0 ? 'X' : (pauliIndex == 1 ? 'Y' : 'Z');
        errors.push_back({gateDistribution(generator), qubitDistribution(generator), pauli});
    }

    std::stable_sort(errors.begin(), errors.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.afterGate < rhs.afterGate;
    });
    return errors;
}

qc::QuantumComputation buildFaultyCircuit(
    const qc::QuantumComputation& original,
    const std::vector<InjectedPauliError>& errors) {
    qc::QuantumComputation faulty(original.getNqubits(), original.getNcbits());
    faulty = qc::QuantumComputation(original);
    faulty.clear();

    std::size_t errorIndex = 0;
    for (std::size_t gateIndex = 0; gateIndex < original.getNops(); ++gateIndex) {
        faulty.emplace_back(original.at(gateIndex)->clone());
        while (errorIndex < errors.size() && errors[errorIndex].afterGate == gateIndex) {
            const auto qubit = static_cast<qc::Qubit>(errors[errorIndex].qubit);
            switch (errors[errorIndex].pauli) {
                case 'X':
                    faulty.x(qubit);
                    break;
                case 'Y':
                    faulty.y(qubit);
                    break;
                case 'Z':
                    faulty.z(qubit);
                    break;
                default:
                    throw std::logic_error("Unsupported sampled Pauli");
            }
            ++errorIndex;
        }
    }

    return faulty;
}

void appendShiftedOperation(
    qc::QuantumComputation& destination,
    const qc::Operation& operation,
    const qc::Qubit shift) {
    auto shiftedOperation = operation.clone();

    auto shiftedTargets = shiftedOperation->getTargets();
    for (auto& target : shiftedTargets) {
        target = static_cast<qc::Qubit>(target + shift);
    }
    shiftedOperation->setTargets(shiftedTargets);

    qc::Controls shiftedControls;
    for (const auto& control : shiftedOperation->getControls()) {
        shiftedControls.emplace(static_cast<qc::Qubit>(control.qubit + shift), control.type);
    }
    shiftedOperation->setControls(shiftedControls);
    shiftedOperation->setNqubits(destination.getNqubits());
    destination.emplace_back(std::move(shiftedOperation));
}

qc::QuantumComputation buildEvaluationCircuit(
    const qc::QuantumComputation& original,
    const qc::QuantumComputation& faulty) {
    qc::QuantumComputation inverseOriginal(original);
    inverseOriginal.invert();

    const auto systemQubits = static_cast<qc::Qubit>(original.getNqubits());
    qc::QuantumComputation evaluation(2 * original.getNqubits());

    for (qc::Qubit qubit = 0; qubit < systemQubits; ++qubit) {
        evaluation.h(qubit);
        evaluation.cx(qubit, static_cast<qc::Qubit>(systemQubits + qubit));
    }

    for (const auto& op : inverseOriginal) {
        appendShiftedOperation(evaluation, *op, systemQubits);
    }
    for (const auto& op : faulty) {
        appendShiftedOperation(evaluation, *op, systemQubits);
    }

    for (qc::Qubit qubit = 0; qubit < systemQubits; ++qubit) {
        evaluation.cx(qubit, static_cast<qc::Qubit>(systemQubits + qubit));
        evaluation.h(qubit);
    }

    return evaluation;
}

TDD makeZeroState(const int qubitCount, dd::Package<>* ddpackage) {
    TensorNetwork tensorNetwork;
    for (int qubit = 0; qubit < qubitCount; ++qubit) {
        const auto array = stateToArray(BasisStates::zero);
        tensorNetwork.add_ts(Tensor(array, {{"x" + std::to_string(qubit) + "_0", 0}}));
    }
    return tensorNetwork.cont(ddpackage);
}

ContractionStats contractWithStats(dd::TensorNetwork* tensorNetwork, dd::Package<>* ddpackage, const int qubitCount) {
    if (!ddpackage) {
        throw std::runtime_error("ddpackage is null");
    }
    if (tensorNetwork->tensors.empty()) {
        throw std::runtime_error("null tensor network");
    }

    clock_t start = clock();
    TDD result = makeZeroState(qubitCount, ddpackage);
    ddpackage->incRef(result.e);
    unsigned int maxNode = ddpackage->size(result.e);

    for (std::size_t index = 0; index < tensorNetwork->tensors.size(); ++index) {
        TDD current = tensorNetwork->tensors[index].to_tdd(ddpackage);
        TDD next = ddpackage->cont(result, current);
        ddpackage->incRef(next.e);
        ddpackage->decRef(result.e);
        ddpackage->garbageCollect();
        result = next;
        maxNode = std::max(maxNode, ddpackage->size(result.e));
    }

    clock_t end = clock();
    return {
        result,
        maxNode,
        static_cast<double>(end - start) / CLOCKS_PER_SEC,
    };
}

Edge<mNode> sliceStateEdge(const Edge<mNode>& edge, const int variable, const int value, dd::Package<>* ddpackage) {
    if (edge.w == Complex::zero) {
        return edge;
    }
    if (edge.p->v == -1 || edge.p->v < variable) {
        return edge;
    }
    if (edge.p->v != variable) {
        throw std::runtime_error("sliceStateEdge only supports direct slicing on the current variable");
    }

    if (edge.p->v != edge.map->level) {
        auto next = edge.p->e[value];
        if (next.w != Complex::zero) {
            next.w = ddpackage->cn.mulCached(next.w, edge.w);
            next.map = ddpackage->mapmul(edge.map, next.map);
            ddpackage->cn.mul(next.w, next.w, ddpackage->cn.getTemporary(cos(next.map->extra_phase * rotate_angle), sin(next.map->extra_phase * rotate_angle)));
        }
        return next;
    }

    if (edge.map->x == 0) {
        auto next = edge.p->e[value];
        if (next.w != Complex::zero) {
            next.w = ddpackage->cn.mulCached(next.w, edge.w);
            next.map = ddpackage->mapmul(edge.map->father, next.map);
            ddpackage->cn.mul(next.w, next.w, ddpackage->cn.getTemporary(cos(next.map->extra_phase * rotate_angle), sin(next.map->extra_phase * rotate_angle)));
            if (value == 1) {
                ddpackage->cn.mul(next.w, next.w, ddpackage->cn.getTemporary(cos(edge.map->rotate * rotate_angle), sin(edge.map->rotate * rotate_angle)));
            }
        }
        return next;
    }

    auto next = edge.p->e[1 - value];
    if (next.w != Complex::zero) {
        next.w = ddpackage->cn.mulCached(next.w, edge.w);
        next.map = ddpackage->mapmul(edge.map->father, next.map);
        ddpackage->cn.mul(next.w, next.w, ddpackage->cn.getTemporary(cos(next.map->extra_phase * rotate_angle), sin(next.map->extra_phase * rotate_angle)));
        if (value == 0) {
            ddpackage->cn.mul(next.w, next.w, ddpackage->cn.getTemporary(cos(edge.map->rotate * rotate_angle), sin(edge.map->rotate * rotate_angle)));
        }
    }
    return next;
}

Complex amplitudeForBitstring(const TDD& tdd, const std::string& basisState, dd::Package<>* ddpackage) {
    auto edge = tdd.e;
    for (const auto bitChar : basisState) {
        if (edge.p->v == -1) {
            break;
        }
        const auto bit = bitChar == '1' ? 1 : 0;
        edge = sliceStateEdge(edge, edge.p->v, bit, ddpackage);
    }
    return edge.w;
}

double squaredMagnitude(const dd::Complex& value) {
    const double real = dd::CTEntry::val(value.r);
    const double imag = dd::CTEntry::val(value.i);
    return real * real + imag * imag;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <qasm-file> <error-count> <seed>\n";
        return 1;
    }

    const std::string qasmPath = argv[1];
    const std::size_t errorCount = static_cast<std::size_t>(std::stoull(argv[2]));
    const std::uint64_t seed = static_cast<std::uint64_t>(std::stoull(argv[3]));

    try {
        const auto qasm = readFile(qasmPath);
        const auto originalCircuit = qc::QuantumComputation::fromQASM(qasm);
        const auto sampledErrors = sampleInjectedErrors(originalCircuit, errorCount, seed);
        const auto faultyCircuit = buildFaultyCircuit(originalCircuit, sampledErrors);
        const auto evaluationCircuit = buildEvaluationCircuit(originalCircuit, faultyCircuit);

        auto evaluationCircuitPtr = std::make_shared<qc::QuantumComputation>(evaluationCircuit);
        auto ddPack = std::make_shared<dd::Package<>>(3 * evaluationCircuitPtr->getNqubits());
        auto tensorNetwork = cir_2_tn(evaluationCircuitPtr, ddPack);
        const auto stats = contractWithStats(&tensorNetwork, ddPack.get(), static_cast<int>(evaluationCircuitPtr->getNqubits()));

        const auto overlap = amplitudeForBitstring(
            stats.tdd,
            std::string(evaluationCircuitPtr->getNqubits(), '0'),
            ddPack.get());

        const double dimension = static_cast<double>(std::uint64_t{1} << originalCircuit.getNqubits());
        const double traceReal = dd::CTEntry::val(overlap.r) * dimension;
        const double traceImag = dd::CTEntry::val(overlap.i) * dimension;
        const double fidelity = squaredMagnitude(overlap);

        std::cout << std::setprecision(17);
        std::cout << "FIDELITY_BEGIN\n";
        std::cout << "qasm\t" << qasmPath << "\n";
        std::cout << "qubits\t" << originalCircuit.getNqubits() << "\n";
        std::cout << "original_gates\t" << originalCircuit.getNops() << "\n";
        std::cout << "faulty_gates\t" << faultyCircuit.getNops() << "\n";
        std::cout << "evaluation_qubits\t" << evaluationCircuitPtr->getNqubits() << "\n";
        std::cout << "evaluation_gates\t" << evaluationCircuitPtr->getNops() << "\n";
        std::cout << "error_count\t" << errorCount << "\n";
        std::cout << "seed\t" << seed << "\n";
        for (const auto& error : sampledErrors) {
            std::cout << "ERROR\t" << error.afterGate << "\t" << error.qubit << "\t" << error.pauli << "\n";
        }
        std::cout << "trace_re\t" << traceReal << "\n";
        std::cout << "trace_im\t" << traceImag << "\n";
        std::cout << "fidelity\t" << fidelity << "\n";
        std::cout << "max_nodes\t" << stats.maxNode << "\n";
        std::cout << "time_s\t" << stats.timeSeconds << "\n";
        std::cout << "FIDELITY_END\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "test_fidelity failed: " << exception.what() << std::endl;
        return 2;
    }
}