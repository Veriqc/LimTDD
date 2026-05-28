#include "QuantumComputation.hpp"
#include "Cir_import.h"
#include "dd/Export.hpp"
#include "dd/Tensor.hpp"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>

using namespace dd;
xt::xarray<dd::ComplexValue> stateToArray(const BasisStates& state){
    switch (state) {
                    case BasisStates::zero:
                        return {complex_one,complex_zero};
                    case BasisStates::one:
                        return {complex_zero,complex_one};
                    case BasisStates::plus:
                        return {complex_SQRT2_2,complex_SQRT2_2};
                    case BasisStates::minus:
                        return {complex_SQRT2_2,complex_mSQRT2_2};
                    case BasisStates::right:
                        return {complex_SQRT2_2,complex_iSQRT2_2};
                    case BasisStates::left:
                        return {complex_SQRT2_2,complex_miSQRT2_2};
    }
}
TDD makezero(int n, dd::Package<>* ddpackage, std::vector<BasisStates> states) {
    TensorNetwork tn;
    if(n> states.size()){
        throw std::invalid_argument("wrong qubit number");
    }
    for(int i=0; i < n; i++){
        xt::xarray<dd::ComplexValue> array = stateToArray(states[i]);
        Tensor temp = Tensor(array,{{"x"+std::to_string(i)+"_0",0}});
        tn.add_ts(temp);
    }
    return tn.cont(ddpackage);
}

bool stepTraceEnabled();
bool stepInTraceWindow(std::size_t step);
void printRegressionDeltaStep(
    const char* phase,
    std::size_t step,
    std::size_t nodesBefore,
    std::size_t nodesAfter,
    const dd::Package<>::RegressionDiagnostics& before,
    const dd::Package<>::RegressionDiagnostics& after);

TDD cont(dd::TensorNetwork* tn,dd::Package<>* ddpackage, int n,bool simulate,const std::vector<BasisStates>& states,bool release = true) {
    if (!ddpackage) {
        throw std::runtime_error("ddpackage is null");
    }
    if (tn->tensors.size() == 0) {
        throw std::runtime_error("null tensor network");
    }

    clock_t start,end;
    start = clock();
    TDD res_dd = simulate ? makezero(n, ddpackage, states) : tn->tensors[0].to_tdd(ddpackage);
    ddpackage->incRef(res_dd.e);
    unsigned int MAX_NODE = ddpackage->size(res_dd.e);

    // The loop starts from 0 if simulating, 1 otherwise.
    for (size_t i = simulate ? 0 : 1; i < tn->tensors.size(); ++i) {
        try {
            const auto traceThisStep = ddpackage->enableRegressionDiagnostics && stepTraceEnabled() && stepInTraceWindow(i);
            const auto stateNodesBefore = traceThisStep ? ddpackage->size(res_dd.e) : 0U;
            const auto diagnosticsBeforeTensor = traceThisStep ? ddpackage->regressionDiagnostics : dd::Package<>::RegressionDiagnostics{};
            TDD current_dd = tn->tensors[i].to_tdd(ddpackage);
            const auto diagnosticsAfterTensor = traceThisStep ? ddpackage->regressionDiagnostics : dd::Package<>::RegressionDiagnostics{};
            ddpackage->setContStageTrace(traceThisStep, i);
            TDD temp_dd = ddpackage->cont(res_dd, current_dd);
            ddpackage->setContStageTrace(false);
            if (traceThisStep) {
                printRegressionDeltaStep("tensor", i, 0U, ddpackage->size(current_dd.e), diagnosticsBeforeTensor, diagnosticsAfterTensor);
                printRegressionDeltaStep("cont", i, stateNodesBefore, ddpackage->size(temp_dd.e), diagnosticsAfterTensor, ddpackage->regressionDiagnostics);
            }
            if (release) {
                ddpackage->incRef(temp_dd.e);
                ddpackage->decRef(res_dd.e);
                ddpackage->garbageCollect();
            }
            res_dd = temp_dd;
            MAX_NODE = std::max(MAX_NODE, ddpackage->size(res_dd.e));
        } catch (...) {
            std::exception_ptr p = std::current_exception();
            // std::clog << (p ? p.__cxa_exception_type()->name() : "null ") << std::endl;
        }
    }
    end = clock();
    std::cout<<"time: " << double(end-start)/CLOCKS_PER_SEC << "s" <<std::endl;
    std::cout<<"MAX node: " << MAX_NODE  <<std::endl;

    return res_dd;
};


BasisStates charToBasisState(char c) {
    static const std::map<char, BasisStates> stateMap = {
        {'0', BasisStates::zero},
        {'1', BasisStates::one},
        {'+', BasisStates::plus},
        {'-', BasisStates::minus},
        {'>', BasisStates::right},
        {'<', BasisStates::left}
    };

    auto it = stateMap.find(c);
    if (it != stateMap.end()) {
        return it->second;
    } else {
        throw std::invalid_argument("Invalid basis state character: " + std::string(1, c));
    }
}
std::vector<BasisStates> stringToBasisStates(const std::string& states) {
    std::vector<BasisStates> basisStates;
    for (char c : states) {
        basisStates.push_back(charToBasisState(c));
    }
    return basisStates;
}

bool envFlagEnabled(const char* name) {
    const auto* value = std::getenv(name);
    return value != nullptr && std::string_view(value) == "1";
}

std::optional<std::size_t> envSizeValue(const char* name) {
    const auto* value = std::getenv(name);
    if (value == nullptr || *value == '\0') {
        return std::nullopt;
    }
    return static_cast<std::size_t>(std::stoull(value));
}

std::string formatGateControls(const qc::Controls& controls) {
    std::ostringstream stream;
    bool first = true;
    for (const auto& control : controls) {
        if (!first) {
            stream << ",";
        }
        first = false;
        stream << (control.type == qc::Control::Type::Neg ? "!" : "") << control.qubit;
    }
    return stream.str();
}

std::string formatGateTargets(const qc::Targets& targets) {
    std::ostringstream stream;
    for (std::size_t i = 0; i < targets.size(); ++i) {
        if (i != 0) {
            stream << ",";
        }
        stream << targets[i];
    }
    return stream.str();
}

bool stepTraceEnabled() {
    return envSizeValue("LIMTDD_STEP_TRACE_START").has_value() ||
           envSizeValue("LIMTDD_STEP_TRACE_LEN").has_value();
}

bool stepInTraceWindow(const std::size_t step) {
    const auto start = envSizeValue("LIMTDD_STEP_TRACE_START").value_or(0);
    const auto length = envSizeValue("LIMTDD_STEP_TRACE_LEN").value_or(16);
    return step >= start && step < start + length;
}

void appendRegressionDelta(std::ostream& output, const char* label, const std::size_t before, const std::size_t after) {
    if (after != before) {
        output << " " << label << "+=" << (after - before);
    }
}

void printRegressionDeltaStep(
    const char* phase,
    const std::size_t step,
    const std::size_t nodesBefore,
    const std::size_t nodesAfter,
    const dd::Package<>::RegressionDiagnostics& before,
    const dd::Package<>::RegressionDiagnostics& after) {
    std::ostringstream output;
    output << "step_" << phase << "[" << step << "]:"
           << " nodes=" << nodesBefore << "->" << nodesAfter;
    appendRegressionDelta(output, "normalize.zero_children", before.normalizeZeroChildren, after.normalizeZeroChildren);
    appendRegressionDelta(output, "normalize.child_phase_adds", before.normalizeChildPhaseAdds, after.normalizeChildPhaseAdds);
    appendRegressionDelta(output, "normalize.root_phase_promotions", before.normalizeRootPhasePromotions, after.normalizeRootPhasePromotions);
    appendRegressionDelta(output, "normalize.child0_from_input0", before.normalizeChild0FromInput0, after.normalizeChild0FromInput0);
    appendRegressionDelta(output, "normalize.child0_from_input1", before.normalizeChild0FromInput1, after.normalizeChild0FromInput1);
    appendRegressionDelta(output, "normalize.child1_from_input0", before.normalizeChild1FromInput0, after.normalizeChild1FromInput0);
    appendRegressionDelta(output, "normalize.child1_from_input1", before.normalizeChild1FromInput1, after.normalizeChild1FromInput1);
    appendRegressionDelta(output, "normalize.child1_zeroed", before.normalizeChild1Zeroed, after.normalizeChild1Zeroed);
    appendRegressionDelta(output, "normalize.child1_same_weight_diff_map", before.normalizeChild1SameWeightDifferentMap, after.normalizeChild1SameWeightDifferentMap);
    appendRegressionDelta(output, "normalize.child1_same_map_diff_weight", before.normalizeChild1SameMapDifferentWeight, after.normalizeChild1SameMapDifferentWeight);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight", before.normalizeChild1DifferentMapAndWeight, after.normalizeChild1DifferentMapAndWeight);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_same_node", before.normalizeChild1DiffMapAndWeightSameNode, after.normalizeChild1DiffMapAndWeightSameNode);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_different_node", before.normalizeChild1DiffMapAndWeightDifferentNode, after.normalizeChild1DiffMapAndWeightDifferentNode);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_angle_snapped", before.normalizeChild1DiffMapAndWeightAngleSnapped, after.normalizeChild1DiffMapAndWeightAngleSnapped);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_residual_weight", before.normalizeChild1DiffMapAndWeightResidualWeight, after.normalizeChild1DiffMapAndWeightResidualWeight);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_rot_nonzero", before.normalizeChild1DiffMapAndWeightRotNonZero, after.normalizeChild1DiffMapAndWeightRotNonZero);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_near_snap_boundary", before.normalizeChild1DiffMapAndWeightNearSnapBoundary, after.normalizeChild1DiffMapAndWeightNearSnapBoundary);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_snap_interior", before.normalizeChild1DiffMapAndWeightSnapInterior, after.normalizeChild1DiffMapAndWeightSnapInterior);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_snap_boundary_inside", before.normalizeChild1DiffMapAndWeightSnapBoundaryInside, after.normalizeChild1DiffMapAndWeightSnapBoundaryInside);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_residual_boundary_outside", before.normalizeChild1DiffMapAndWeightResidualBoundaryOutside, after.normalizeChild1DiffMapAndWeightResidualBoundaryOutside);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_residual_far_outside", before.normalizeChild1DiffMapAndWeightResidualFarOutside, after.normalizeChild1DiffMapAndWeightResidualFarOutside);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_rot_positive", before.normalizeChild1DiffMapAndWeightRotPositive, after.normalizeChild1DiffMapAndWeightRotPositive);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_rot_negative", before.normalizeChild1DiffMapAndWeightRotNegative, after.normalizeChild1DiffMapAndWeightRotNegative);
    appendRegressionDelta(output, "normalize.child1_diff_map_and_weight_rot_zero", before.normalizeChild1DiffMapAndWeightRotZero, after.normalizeChild1DiffMapAndWeightRotZero);
    appendRegressionDelta(output, "normalize.child1_legacy_abs_snap_disagreement", before.normalizeChild1LegacyAbsSnapDisagreement, after.normalizeChild1LegacyAbsSnapDisagreement);
    appendRegressionDelta(output, "normalize.child1_phaseful_after_mapdiv", before.normalizeChild1PhasefulAfterMapdiv, after.normalizeChild1PhasefulAfterMapdiv);
    appendRegressionDelta(output, "normalize.root_phase_promotions_from_input0", before.normalizeRootPhasePromotionsFromInput0, after.normalizeRootPhasePromotionsFromInput0);
    appendRegressionDelta(output, "normalize.root_phase_promotions_from_input1", before.normalizeRootPhasePromotionsFromInput1, after.normalizeRootPhasePromotionsFromInput1);
    appendRegressionDelta(output, "unique.semantic_alternatives", before.uniqueSemanticAlternatives, after.uniqueSemanticAlternatives);
    appendRegressionDelta(output, "unique.semantic_alternative_misses", before.uniqueSemanticAlternativeMisses, after.uniqueSemanticAlternativeMisses);
    appendRegressionDelta(output, "unique.semantic_global_alternatives", before.uniqueSemanticGlobalAlternatives, after.uniqueSemanticGlobalAlternatives);
    appendRegressionDelta(output, "unique.semantic_global_alternative_misses", before.uniqueSemanticGlobalAlternativeMisses, after.uniqueSemanticGlobalAlternativeMisses);
    appendRegressionDelta(output, "unique.structural_alternatives", before.uniqueStructuralAlternatives, after.uniqueStructuralAlternatives);
    appendRegressionDelta(output, "unique.structural_alternative_misses", before.uniqueStructuralAlternativeMisses, after.uniqueStructuralAlternativeMisses);
    appendRegressionDelta(output, "unique.structural_same_weight_diff_map", before.uniqueStructuralSameWeightDifferentMap, after.uniqueStructuralSameWeightDifferentMap);
    appendRegressionDelta(output, "unique.structural_same_map_diff_weight", before.uniqueStructuralSameMapDifferentWeight, after.uniqueStructuralSameMapDifferentWeight);
    appendRegressionDelta(output, "unique.structural_diff_map_and_weight", before.uniqueStructuralDifferentMapAndWeight, after.uniqueStructuralDifferentMapAndWeight);
    appendRegressionDelta(output, "unique.pre_normalize_structural_alternatives", before.uniquePreNormalizeStructuralAlternatives, after.uniquePreNormalizeStructuralAlternatives);
    appendRegressionDelta(output, "unique.pre_normalize_same_weight_diff_map", before.uniquePreNormalizeSameWeightDifferentMap, after.uniquePreNormalizeSameWeightDifferentMap);
    appendRegressionDelta(output, "unique.pre_normalize_same_map_diff_weight", before.uniquePreNormalizeSameMapDifferentWeight, after.uniquePreNormalizeSameMapDifferentWeight);
    appendRegressionDelta(output, "unique.pre_normalize_diff_map_and_weight", before.uniquePreNormalizeDifferentMapAndWeight, after.uniquePreNormalizeDifferentMapAndWeight);
    appendRegressionDelta(output, "unique.normalize_transition_diff_both_to_map_only", before.uniqueNormalizeTransitionDiffBothToMapOnly, after.uniqueNormalizeTransitionDiffBothToMapOnly);
    appendRegressionDelta(output, "unique.normalize_transition_diff_both_to_weight_only", before.uniqueNormalizeTransitionDiffBothToWeightOnly, after.uniqueNormalizeTransitionDiffBothToWeightOnly);
    appendRegressionDelta(output, "unique.normalize_transition_diff_both_to_diff_both", before.uniqueNormalizeTransitionDiffBothToDiffBoth, after.uniqueNormalizeTransitionDiffBothToDiffBoth);
    appendRegressionDelta(output, "tadd.same_pointer_map_mismatch", before.taddSamePointerMapMismatch, after.taddSamePointerMapMismatch);
    appendRegressionDelta(output, "tadd.mismatch_residual_header", before.taddMismatchResidualHeader, after.taddMismatchResidualHeader);
    appendRegressionDelta(output, "tadd.mismatch_residual_non_header", before.taddMismatchResidualNonHeader, after.taddMismatchResidualNonHeader);
    appendRegressionDelta(output, "tadd.mismatch_residual_phaseful", before.taddMismatchResidualPhaseful, after.taddMismatchResidualPhaseful);
    appendRegressionDelta(output, "tadd.mismatch_add_hits", before.taddMismatchAddHits, after.taddMismatchAddHits);
    appendRegressionDelta(output, "tadd.mismatch_add_misses", before.taddMismatchAddMisses, after.taddMismatchAddMisses);
    appendRegressionDelta(output, "tadd.mismatch_add_miss_empty", before.taddMismatchAddMissEmpty, after.taddMismatchAddMissEmpty);
    appendRegressionDelta(output, "tadd.mismatch_add_miss_map_only", before.taddMismatchAddMissMapOnly, after.taddMismatchAddMissMapOnly);
    appendRegressionDelta(output, "tadd.mismatch_add_miss_weight_only", before.taddMismatchAddMissWeightOnly, after.taddMismatchAddMissWeightOnly);
    appendRegressionDelta(output, "tadd.mismatch_add_miss_map_and_weight", before.taddMismatchAddMissMapAndWeight, after.taddMismatchAddMissMapAndWeight);
    appendRegressionDelta(output, "mapmul.base_reset_self", before.mapmulBaseResetSelf, after.mapmulBaseResetSelf);
    appendRegressionDelta(output, "mapmul.base_reset_other", before.mapmulBaseResetOther, after.mapmulBaseResetOther);
    appendRegressionDelta(output, "mapmul.lookup_hits", before.mapmulLookupHits, after.mapmulLookupHits);
    appendRegressionDelta(output, "mapmul.lookup_phaseful", before.mapmulLookupPhaseful, after.mapmulLookupPhaseful);
    appendRegressionDelta(output, "mapmul.result_phaseful", before.mapmulResultPhaseful, after.mapmulResultPhaseful);
    appendRegressionDelta(output, "mapdiv.base_reset_self", before.mapdivBaseResetSelf, after.mapdivBaseResetSelf);
    appendRegressionDelta(output, "mapdiv.header_reset", before.mapdivHeaderReset, after.mapdivHeaderReset);
    appendRegressionDelta(output, "mapdiv.lookup_hits", before.mapdivLookupHits, after.mapdivLookupHits);
    appendRegressionDelta(output, "mapdiv.lookup_phaseful", before.mapdivLookupPhaseful, after.mapdivLookupPhaseful);
    appendRegressionDelta(output, "mapdiv.lookup_phase_overwrite", before.mapdivLookupPhaseOverwrite, after.mapdivLookupPhaseOverwrite);
    appendRegressionDelta(output, "mapdiv.lookup_phase_overwrite_non_header", before.mapdivLookupPhaseOverwriteNonHeader, after.mapdivLookupPhaseOverwriteNonHeader);
    appendRegressionDelta(output, "mapdiv.result_phaseful", before.mapdivResultPhaseful, after.mapdivResultPhaseful);
    appendRegressionDelta(output, "find_remain.phase_carries", before.findRemainPhaseCarries, after.findRemainPhaseCarries);
    std::cout << output.str() << std::endl;
}

void printGateWindow(const qc::QuantumComputation& qc) {
    const auto windowStart = envSizeValue("LIMTDD_GATE_WINDOW_START");
    const auto windowLength = envSizeValue("LIMTDD_GATE_WINDOW_LEN");
    if (!windowStart.has_value() && !windowLength.has_value()) {
        return;
    }

    const auto start = std::min(windowStart.value_or(0), qc.getNops());
    const auto length = windowLength.value_or(32);
    const auto end = std::min(start + length, qc.getNops());

    std::cout << "gate_window: [" << start << ", " << end << ")/" << qc.getNops() << std::endl;
    for (std::size_t index = start; index < end; ++index) {
        const auto& op = qc.at(index);
        std::cout << "gate[" << index << "]:"
                  << " name=" << op->getName()
                  << " type=" << qc::toString(op->getType())
                  << " controls=[" << formatGateControls(op->getControls()) << "]"
                  << " targets=[" << formatGateTargets(op->getTargets()) << "]"
                  << std::endl;
    }
}

int runXarraySelftest() {
    dd::ComplexValue one = {1, 0};
    dd::ComplexValue zero = {0, 0};
    dd::ComplexValue two = {2, 0};
    dd::ComplexValue three = {3, 0};

    auto ddPack = std::make_shared<dd::Package<>>(10);
    ddPack->varOrder = {{"x0", 0}, {"y0", 1}, {"x1", 2}, {"y1", 3}};
    xt::xarray<dd::ComplexValue> tensorCnot = {
        {{{one, zero}, {zero, one}}, {{zero, zero}, {zero, zero}}},
        {{{zero, zero}, {zero, zero}}, {{zero, one}, {one, zero}}},
    };
    xt::xarray<dd::ComplexValue> permutedTensorCnot = {
        {{{zero, zero}, {zero, zero}}, {{zero, zero}, {zero, zero}}},
        {{{zero, zero}, {zero, zero}}, {{zero, zero}, {zero, zero}}},
    };
    for (std::size_t a = 0; a < 2; ++a) {
        for (std::size_t b = 0; b < 2; ++b) {
            for (std::size_t c = 0; c < 2; ++c) {
                for (std::size_t d = 0; d < 2; ++d) {
                    permutedTensorCnot(a, b, c, d) = tensorCnot(c, d, a, b);
                }
            }
        }
    }

    std::vector<dd::Index> tensorIndices = {{"x0", 0}, {"y0", 0}, {"x1", 0}, {"y1", 0}};
    std::vector<dd::Index> permutedIndices = {{"x1", 0}, {"y1", 0}, {"x0", 0}, {"y0", 0}};

    auto tensorTdd = dd::Tensor(tensorCnot, tensorIndices, "tensor_cnot").to_tdd(ddPack.get());
    auto permutedTdd = dd::Tensor(permutedTensorCnot, permutedIndices, "permuted_tensor_cnot").to_tdd(ddPack.get());

    const bool orderEqual = tensorTdd.e == permutedTdd.e;
    std::cout << "xarray_selftest.equal: " << orderEqual << std::endl;
    std::cout << "xarray_selftest.tensor_nodes: " << ddPack->size(tensorTdd.e) << std::endl;
    std::cout << "xarray_selftest.permuted_nodes: " << ddPack->size(permutedTdd.e) << std::endl;

    auto hyperPack = std::make_shared<dd::Package<>>(10);
    hyperPack->varOrder = {{"x0", 0}, {"x1", 2}, {"y1", 3}};
    xt::xarray<dd::ComplexValue> repeatedIndexTensor = tensorCnot;
    xt::xarray<dd::ComplexValue> offDiagonalPerturbed = tensorCnot;
    offDiagonalPerturbed(0, 1, 0, 0) = two;
    offDiagonalPerturbed(0, 1, 1, 1) = three;
    offDiagonalPerturbed(1, 0, 0, 1) = three;
    offDiagonalPerturbed(1, 0, 1, 0) = two;

    std::vector<dd::Index> repeatedIndices = {{"x0", 0}, {"x0", 1}, {"x1", 0}, {"y1", 0}};
    auto repeatedTdd = dd::Tensor(repeatedIndexTensor, repeatedIndices, "repeated_index_base").to_tdd(hyperPack.get());
    auto perturbedTdd = dd::Tensor(offDiagonalPerturbed, repeatedIndices, "repeated_index_perturbed").to_tdd(hyperPack.get());
    const bool repeatedIndexEqual = repeatedTdd.e == perturbedTdd.e;
    std::cout << "xarray_selftest.repeated_index_equal: " << repeatedIndexEqual << std::endl;
    std::cout << "xarray_selftest.repeated_index_nodes: " << hyperPack->size(repeatedTdd.e) << std::endl;
    std::cout << "xarray_selftest.repeated_index_perturbed_nodes: " << hyperPack->size(perturbedTdd.e) << std::endl;

    return (orderEqual && repeatedIndexEqual) ? 0 : 2;
}

int main(int argc, char *argv[]) {
    if (envFlagEnabled("LIMTDD_XARRAY_SELFTEST")) {
        return runXarraySelftest();
    }

    // filename, initial state
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <number>\n";
        return 1;
    }
    
    
    std::string filename = argv[1];
	std::cout << filename << std::endl;
    
    std::ifstream fileStream(filename);
    std::stringstream buffer;
    buffer << fileStream.rdbuf();
    std::string fileContent = buffer.str();
    fileStream.close();

    // Use the file content with QuantumComputation::fromQASM
    const auto qc = qc::QuantumComputation::fromQASM(fileContent);
    std::shared_ptr<qc::QuantumComputation> QC = std::make_shared<qc::QuantumComputation>(std::move(qc));
    printGateWindow(*QC);
    auto ddPack = std::make_shared<dd::Package<>>(3*QC->getNqubits());
    ddPack->enableRegressionDiagnostics = envFlagEnabled("LIMTDD_REGRESSION_DIAG");
    ddPack->disableMapdivLookupWriteback = envFlagEnabled("LIMTDD_DISABLE_MAPDIV_LOOKUP_WRITEBACK");
    ddPack->enableUniqueSemanticProbe = envFlagEnabled("LIMTDD_UNIQUE_SEMANTIC_PROBE");
    auto tn = cir_2_tn(QC,ddPack);
    if (const auto prefixLimit = envSizeValue("LIMTDD_TN_PREFIX"); prefixLimit.has_value()) {
        const auto originalSize = tn.tensors.size();
        const auto effectiveSize = std::min(prefixLimit.value(), originalSize);
        tn.tensors.erase(tn.tensors.begin() + static_cast<std::ptrdiff_t>(effectiveSize), tn.tensors.end());
        std::cout << "tensor_prefix: " << effectiveSize << "/" << originalSize << std::endl;
    }

    bool simulate = false;
    std::vector<BasisStates> initialStates;
    if(argc > 2){
        simulate = true;
         try {
            initialStates = stringToBasisStates(argv[2]);
            std::cout << "Initial states vector size: " << initialStates.size() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }
    std::cout <<"simulate:" << simulate << std::endl;

	dd::TDD tdd = cont(&tn,ddPack.get(),QC->getNqubits(),simulate, initialStates);
    // dd::export2Dot(tdd.e,"test",true,true);
    
    std::cout<<"final node: " << ddPack->size(tdd.e) <<std::endl;
    if (ddPack->enableRegressionDiagnostics) {
        ddPack->printRegressionDiagnostics();
    }
    return 0;
}
