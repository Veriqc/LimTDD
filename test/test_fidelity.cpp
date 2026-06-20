#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
#include "dd/Tensor.hpp"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstdlib>
#include <cmath>
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

bool envFlagEnabled(const char* name) {
    const auto* value = std::getenv(name);
    if (value == nullptr) {
        return false;
    }
    const std::string flag(value);
    return flag == "1" || flag == "true" || flag == "TRUE" || flag == "on" || flag == "ON";
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

qc::QuantumComputation takeCircuitPrefix(
    const qc::QuantumComputation& original,
    const std::size_t gateCount) {
    qc::QuantumComputation prefix(original.getNqubits(), original.getNcbits());
    prefix = qc::QuantumComputation(original);
    prefix.clear();

    const auto keptGates = std::min(gateCount, original.getNops());
    for (std::size_t gateIndex = 0; gateIndex < keptGates; ++gateIndex) {
        prefix.emplace_back(original.at(gateIndex)->clone());
    }

    return prefix;
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

    for (const auto& op : faulty) {
        appendShiftedOperation(evaluation, *op, systemQubits);
    }
    for (const auto& op : inverseOriginal) {
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

    std::size_t contractionLimit = tensorNetwork->tensors.size();
    if (const auto* limitEnv = std::getenv("LIMTDD_FIDELITY_PREFIX")) {
        contractionLimit = std::min(contractionLimit, static_cast<std::size_t>(std::stoull(limitEnv)));
    }
    const bool traceSteps = std::getenv("LIMTDD_FIDELITY_TRACE") != nullptr;
    const auto* traceStepEnv = std::getenv("LIMTDD_FIDELITY_TRACE_STEP");
    const auto tracedStep = traceStepEnv ? static_cast<std::size_t>(std::stoull(traceStepEnv)) : std::numeric_limits<std::size_t>::max();

    for (std::size_t index = 0; index < contractionLimit; ++index) {
        if (traceSteps) {
            std::cerr << "fidelity_step\t" << index << "\tcurrent_nodes\t" << ddpackage->size(result.e) << std::endl;
        }
        bool traceThisContraction = index == tracedStep;
        if (std::getenv("LIMTDD_CONT_INSERT_TRACE") != nullptr) {
            std::size_t insertTraceStart = 0;
            std::size_t insertTraceEnd = std::numeric_limits<std::size_t>::max();
            if (const auto* startEnv = std::getenv("LIMTDD_CONT_INSERT_TRACE_STEP_START")) {
                insertTraceStart = static_cast<std::size_t>(std::stoull(startEnv));
            }
            if (const auto* endEnv = std::getenv("LIMTDD_CONT_INSERT_TRACE_STEP_END")) {
                insertTraceEnd = static_cast<std::size_t>(std::stoull(endEnv));
            }
            traceThisContraction = traceThisContraction || (index >= insertTraceStart && index <= insertTraceEnd);
        }
        ddpackage->setContStageTrace(traceThisContraction, index);
        TDD current = tensorNetwork->tensors[index].to_tdd(ddpackage);
        TDD next = ddpackage->cont(result, current);
        ddpackage->incRef(next.e);
        ddpackage->decRef(result.e);
        ddpackage->garbageCollect();
        result = next;
        maxNode = std::max(maxNode, ddpackage->size(result.e));
    }

    ddpackage->setContStageTrace(false);

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

dd::Complex amplitudeForZeroState(const TDD& tdd, const std::size_t qubitCount, dd::Package<>* ddpackage) {
    auto edge = tdd.e;
    const bool traceAmplitude = std::getenv("LIMTDD_FIDELITY_SLICE_TRACE") != nullptr;
    std::size_t sliceCount = 0;
    for (std::size_t variable = 0; variable < qubitCount; ++variable) {
        if (edge.p->v == -1) {
            break;
        }
        if (traceAmplitude) {
            std::cerr << "slice_step\t" << sliceCount << "\tedge_var\t" << edge.p->v << "\tmap_level\t" << edge.map->level << '\n';
        }
        edge = sliceStateEdge(edge, edge.p->v, 0, ddpackage);
        ++sliceCount;
    }
    if (traceAmplitude) {
        std::cerr << "slice_done\tcount\t" << sliceCount << "\tterminal\t" << (edge.p->v == -1 ? 1 : 0) << "\tfinal_var\t" << edge.p->v << '\n';
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
        auto originalCircuit = qc::QuantumComputation::fromQASM(qasm);
        if (const auto* originalPrefixEnv = std::getenv("LIMTDD_ORIGINAL_GATE_PREFIX")) {
            originalCircuit = takeCircuitPrefix(originalCircuit, static_cast<std::size_t>(std::stoull(originalPrefixEnv)));
        }
        const auto sampledErrors = sampleInjectedErrors(originalCircuit, errorCount, seed);
        const auto faultyCircuit = buildFaultyCircuit(originalCircuit, sampledErrors);
        const auto evaluationCircuit = buildEvaluationCircuit(originalCircuit, faultyCircuit);

        if (std::getenv("LIMTDD_FIDELITY_TRACE") != nullptr) {
            if (const auto* traceStepEnv = std::getenv("LIMTDD_FIDELITY_TRACE_STEP")) {
                const auto tracedStep = static_cast<std::size_t>(std::stoull(traceStepEnv));
                if (tracedStep < evaluationCircuit.getNops()) {
                    std::cerr << "trace_gate\t" << tracedStep << "\t" << *evaluationCircuit.at(tracedStep) << '\n';
                }
                if (tracedStep > 0 && tracedStep - 1 < evaluationCircuit.getNops()) {
                    std::cerr << "trace_gate_prev\t" << (tracedStep - 1) << "\t" << *evaluationCircuit.at(tracedStep - 1) << '\n';
                }
            }
        }

        auto evaluationCircuitPtr = std::make_shared<qc::QuantumComputation>(evaluationCircuit);
        auto ddPack = std::make_shared<dd::Package<>>(3 * evaluationCircuitPtr->getNqubits());
        ddPack->enableRegressionDiagnostics = envFlagEnabled("LIMTDD_REGRESSION_DIAG");
        ddPack->disableMapdivLookupWriteback = envFlagEnabled("LIMTDD_DISABLE_MAPDIV_LOOKUP_WRITEBACK");
        ddPack->enableTailCxRenormExperiment = envFlagEnabled("LIMTDD_EXPERIMENTAL_TAIL_CX_RENORM");
        auto tensorNetwork = cir_2_tn(evaluationCircuitPtr, ddPack);
        const auto stats = contractWithStats(&tensorNetwork, ddPack.get(), static_cast<int>(evaluationCircuitPtr->getNqubits()));

        const auto overlap = amplitudeForZeroState(stats.tdd, evaluationCircuitPtr->getNqubits(), ddPack.get());

        const double dimension = std::ldexp(1.0, static_cast<int>(originalCircuit.getNqubits()));
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