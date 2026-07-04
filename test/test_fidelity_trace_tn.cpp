#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Tensor.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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

struct TraceContractionStats {
    TDD tdd{};
    unsigned int maxNode = 0;
    double timeSeconds = 0.0;
};

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

qc::QuantumComputation buildUDaggerVCircuit(
    const qc::QuantumComputation& original,
    const qc::QuantumComputation& faulty) {
    qc::QuantumComputation inverseOriginal(original);
    inverseOriginal.invert();

    qc::QuantumComputation traceCircuit(original.getNqubits(), original.getNcbits());
    traceCircuit.clear();

    for (const auto& op : inverseOriginal) {
        traceCircuit.emplace_back(op->clone());
    }
    for (const auto& op : faulty) {
        traceCircuit.emplace_back(op->clone());
    }

    return traceCircuit;
}

TraceContractionStats contractClosedNetworkStrict(dd::TensorNetwork& tensorNetwork, dd::Package<>* ddpackage) {
    if (!ddpackage) {
        throw std::runtime_error("ddpackage is null");
    }
    if (tensorNetwork.tensors.empty()) {
        throw std::runtime_error("null tensor network");
    }

    clock_t start = clock();
    TDD result = tensorNetwork.tensors[0].to_tdd(ddpackage);
    ddpackage->incRef(result.e);
    unsigned int maxNode = ddpackage->size(result.e);

    for (std::size_t index = 1; index < tensorNetwork.tensors.size(); ++index) {
        bool traceThisContraction = false;
        if (std::getenv("LIMTDD_CONT_INSERT_TRACE") != nullptr) {
            std::size_t insertTraceStart = 0;
            std::size_t insertTraceEnd = std::numeric_limits<std::size_t>::max();
            if (const auto* startEnv = std::getenv("LIMTDD_CONT_INSERT_TRACE_STEP_START")) {
                insertTraceStart = static_cast<std::size_t>(std::stoull(startEnv));
            }
            if (const auto* endEnv = std::getenv("LIMTDD_CONT_INSERT_TRACE_STEP_END")) {
                insertTraceEnd = static_cast<std::size_t>(std::stoull(endEnv));
            }
            traceThisContraction = (index >= insertTraceStart && index <= insertTraceEnd);
        }
        ddpackage->setContStageTrace(traceThisContraction, index);
        TDD current = tensorNetwork.tensors[index].to_tdd(ddpackage);
        TDD next = ddpackage->cont(result, current);
        ddpackage->incRef(next.e);
        ddpackage->decRef(result.e);
        ddpackage->garbageCollect();
        result = next;
        maxNode = std::max(maxNode, ddpackage->size(result.e));
        if (std::getenv("LIMTDD_TRACE_CONTRACT_STEPS")) {
            std::cerr << "TRACE_STEP\t" << index
                << "\tresult_indexset_sz\t" << result.index_set.size()
                << "\tresult_key2index_sz\t" << result.key_2_index.size()
                << "\tresult_root_var\t" << static_cast<int>(result.e.p->v)
                << "\tresult_w_re\t" << CTEntry::val(result.e.w.r)
                << "\tresult_w_im\t" << CTEntry::val(result.e.w.i)
                << "\tresult_map_level\t" << (result.e.map ? static_cast<int>(result.e.map->level) : -999)
                << "\tnodes\t" << ddpackage->size(result.e)
                << std::endl;
        }
    }

    ddpackage->setContStageTrace(false);

    clock_t end = clock();
    return {
        result,
        maxNode,
        static_cast<double>(end - start) / CLOCKS_PER_SEC,
    };
}

Edge<mNode> sliceEdge(const Edge<mNode>& edge, const int value, dd::Package<>* ddpackage) {
    if (edge.w == Complex::zero) {
        return edge;
    }
    if (edge.p->v == -1) {
        return edge;
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

Complex collapseResidualVariables(const Edge<mNode>& edge, dd::Package<>* ddpackage) {
    if (edge.w == Complex::zero || edge.p->v == -1) {
        return edge.w;
    }

    return collapseResidualVariables(sliceEdge(edge, 0, ddpackage), ddpackage);
}

dd::Complex scalarFromClosedTdd(const TDD& tdd, dd::Package<>* ddpackage) {
    if (!ddpackage) {
        throw std::runtime_error("ddpackage is null");
    }
    if (tdd.e.w == Complex::zero) {
        return tdd.e.w;
    }
    if (!tdd.index_set.empty() || !tdd.key_2_index.empty()) {
        std::ostringstream message;
        message << "trace tensor network did not close: index_set=" << tdd.index_set.size()
                << " key_2_index=" << tdd.key_2_index.size();
        throw std::runtime_error(message.str());
    }

    return tdd.e.w;
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
        if (originalCircuit.getNqubits() == 0) {
            throw std::invalid_argument("zero-qubit circuits are not supported by this prototype");
        }

        const auto sampledErrors = sampleInjectedErrors(originalCircuit, errorCount, seed);
        const auto faultyCircuit = buildFaultyCircuit(originalCircuit, sampledErrors);
        const auto traceCircuit = buildUDaggerVCircuit(originalCircuit, faultyCircuit);

        auto traceCircuitPtr = std::make_shared<qc::QuantumComputation>(traceCircuit);
        auto ddPack = std::make_shared<dd::Package<>>(3 * traceCircuitPtr->getNqubits());
        ddPack->enableRegressionDiagnostics = envFlagEnabled("LIMTDD_REGRESSION_DIAG");
        ddPack->disableMapdivLookupWriteback = envFlagEnabled("LIMTDD_DISABLE_MAPDIV_LOOKUP_WRITEBACK");
        ddPack->enableTailCxRenormExperiment = envFlagEnabled("LIMTDD_EXPERIMENTAL_TAIL_CX_RENORM");

        auto tnWithBoundary = cir_2_tn_with_boundary(traceCircuitPtr, ddPack, true);
        add_trace_delta_tensors(tnWithBoundary);
        const auto stats = contractClosedNetworkStrict(tnWithBoundary.tensorNetwork, ddPack.get());
        const auto trace = scalarFromClosedTdd(stats.tdd, ddPack.get());
        if (stats.tdd.e.p->v != -1) {
            std::cerr << "TRACE_TN_WARNING\tclosed network retained residual root_var\t"
                      << static_cast<int>(stats.tdd.e.p->v)
                      << "\tusing root edge weight as prototype scalar\n";
        }

        const double dimensionSquared = std::ldexp(1.0, 2 * static_cast<int>(originalCircuit.getNqubits()));
        const double fidelity = squaredMagnitude(trace) / dimensionSquared;

        std::cout << std::setprecision(17);
        std::cout << "TRACE_TN_FIDELITY_BEGIN\n";
        std::cout << "qasm\t" << qasmPath << "\n";
        std::cout << "qubits\t" << originalCircuit.getNqubits() << "\n";
        std::cout << "original_gates\t" << originalCircuit.getNops() << "\n";
        std::cout << "faulty_gates\t" << faultyCircuit.getNops() << "\n";
        std::cout << "trace_circuit_gates\t" << traceCircuit.getNops() << "\n";
        std::cout << "trace_tensors\t" << tnWithBoundary.tensorNetwork.tensors.size() << "\n";
        std::cout << "error_count\t" << errorCount << "\n";
        std::cout << "seed\t" << seed << "\n";
        for (std::size_t index = 0; index < sampledErrors.size(); ++index) {
            const auto& error = sampledErrors[index];
            std::cout << "ERROR\t" << index
                      << "\tafter_gate\t" << error.afterGate
                      << "\tqubit\t" << error.qubit
                      << "\tpauli\t" << error.pauli << "\n";
        }
        std::cout << "trace_re\t" << dd::CTEntry::val(trace.r) << "\n";
        std::cout << "trace_im\t" << dd::CTEntry::val(trace.i) << "\n";
        std::cout << "fidelity\t" << fidelity << "\n";
        std::cout << "max_nodes\t" << stats.maxNode << "\n";
        std::cout << "time_s\t" << stats.timeSeconds << "\n";
        std::cout << "TRACE_TN_FIDELITY_END\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 2;
    }

    return 0;
}
