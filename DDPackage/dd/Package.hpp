#pragma once

#include "Complex.hpp"
#include "ComplexCache.hpp"
#include "ComplexNumbers.hpp"
#include "ComplexTable.hpp"
#include "ComplexValue.hpp"
#include "ComputeTable.hpp"
#include "Control.hpp"
#include "Definitions.hpp"
#include "Edge.hpp"
#include "GateMatrixDefinitions.hpp"
#include "Package_fwd.hpp"

#include "UniqueTable.hpp"

#include "Tdd.hpp"
#include "Maps.hpp"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <numeric>

#include <xtensor/containers/xarray.hpp>
#include <xtensor/core/xshape.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/views/xslice.hpp>
#include <xtensor/containers/xfixed.hpp>
#include <xtensor/views/xview.hpp>

namespace dd {

	bool mapCompare(const the_maps* map1, const the_maps* map2) {
		const auto maxArg = dd::PI;

		while (map1->level >= 0 || map2->level >= 0) {
			if (map1->level > map2->level) {
				// if (-std::pow(-1, map1->x)*ComplexNumbers::arg(map1->rotate) > maxArg) return true;
				const int signedRotate = map1->x ? map1->rotate : -map1->rotate;
				if (signedRotate % root_of_unit > int(root_of_unit / 2)) return true;
				
				map1 = (map1->level >= 0) ? map1->father : map1;
			} else if (map2->level > map1->level) {
				// if (ComplexNumbers::arg(map2->rotate) > maxArg) return true;
				if (map2->rotate > root_of_unit/2) return true;
				map2 = (map2->level >= 0) ? map2->father : map2;
			} else { // map1->level == map2->level
				// fp theta1 = ComplexNumbers::arg(map1->rotate);
				// fp theta2 = ComplexNumbers::arg(map2->rotate);
				// fp phaseDiff = theta2 - theta1 * std::pow(-1, map1->x ^ map2->x);
				// if (phaseDiff > maxArg) return true;
				const int phaseDiff = ((map1->x ^ map2->x) ? (map2->rotate + map1->rotate) : (map2->rotate - map1->rotate)) % root_of_unit;
                if (phaseDiff > root_of_unit/2) return true;
				map1 = (map1->level >= 0) ? map1->father : map1;
				map2 = (map2->level >= 0) ? map2->father : map2;
			}
		}

		return false;
	}






	template <class Config> class Package {
		static_assert(std::is_base_of_v<DDPackageConfig, Config>, "Config must be derived from DDPackageConfig");

		///
		/// Complex number handling
		///
	public:
		ComplexNumbers cn{};

		struct RegressionDiagnostics {
			std::size_t normalizeZeroChildren = 0;
			std::size_t normalizeChildPhaseAdds = 0;
			std::size_t normalizeRootPhasePromotions = 0;
			std::size_t normalizeChild0FromInput0 = 0;
			std::size_t normalizeChild0FromInput1 = 0;
			std::size_t normalizeChild1FromInput0 = 0;
			std::size_t normalizeChild1FromInput1 = 0;
			std::size_t normalizeChild1Zeroed = 0;
			std::size_t normalizeChild1SameWeightDifferentMap = 0;
			std::size_t normalizeChild1SameMapDifferentWeight = 0;
			std::size_t normalizeChild1DifferentMapAndWeight = 0;
			std::size_t normalizeChild1DiffMapAndWeightSameNode = 0;
			std::size_t normalizeChild1DiffMapAndWeightDifferentNode = 0;
			std::size_t normalizeChild1DiffMapAndWeightAngleSnapped = 0;
			std::size_t normalizeChild1DiffMapAndWeightResidualWeight = 0;
			std::size_t normalizeChild1DiffMapAndWeightRotNonZero = 0;
			std::size_t normalizeChild1DiffMapAndWeightNearSnapBoundary = 0;
			std::size_t normalizeChild1DiffMapAndWeightSnapInterior = 0;
			std::size_t normalizeChild1DiffMapAndWeightSnapBoundaryInside = 0;
			std::size_t normalizeChild1DiffMapAndWeightResidualBoundaryOutside = 0;
			std::size_t normalizeChild1DiffMapAndWeightResidualFarOutside = 0;
			std::size_t normalizeChild1DiffMapAndWeightRotPositive = 0;
			std::size_t normalizeChild1DiffMapAndWeightRotNegative = 0;
			std::size_t normalizeChild1DiffMapAndWeightRotZero = 0;
			std::size_t normalizeChild1LegacyAbsSnapDisagreement = 0;
			std::size_t normalizeChild1PhasefulAfterMapdiv = 0;
			std::size_t normalizeRootPhasePromotionsFromInput0 = 0;
			std::size_t normalizeRootPhasePromotionsFromInput1 = 0;
			std::size_t uniqueSemanticAlternatives = 0;
			std::size_t uniqueSemanticAlternativeMisses = 0;
			std::size_t uniqueSemanticGlobalAlternatives = 0;
			std::size_t uniqueSemanticGlobalAlternativeMisses = 0;
			std::size_t uniqueStructuralAlternatives = 0;
			std::size_t uniqueStructuralAlternativeMisses = 0;
			std::size_t uniqueStructuralSameWeightDifferentMap = 0;
			std::size_t uniqueStructuralSameMapDifferentWeight = 0;
			std::size_t uniqueStructuralDifferentMapAndWeight = 0;
			std::size_t uniquePreNormalizeStructuralAlternatives = 0;
			std::size_t uniquePreNormalizeSameWeightDifferentMap = 0;
			std::size_t uniquePreNormalizeSameMapDifferentWeight = 0;
			std::size_t uniquePreNormalizeDifferentMapAndWeight = 0;
			std::size_t uniqueNormalizeTransitionDiffBothToMapOnly = 0;
			std::size_t uniqueNormalizeTransitionDiffBothToWeightOnly = 0;
			std::size_t uniqueNormalizeTransitionDiffBothToDiffBoth = 0;
			std::size_t taddSamePointerMapMismatch = 0;
			std::size_t taddMismatchResidualHeader = 0;
			std::size_t taddMismatchResidualNonHeader = 0;
			std::size_t taddMismatchResidualPhaseful = 0;
			std::size_t taddMismatchAddHits = 0;
			std::size_t taddMismatchAddMisses = 0;
			std::size_t taddMismatchAddMissEmpty = 0;
			std::size_t taddMismatchAddMissMapOnly = 0;
			std::size_t taddMismatchAddMissWeightOnly = 0;
			std::size_t taddMismatchAddMissMapAndWeight = 0;
			std::size_t mapmulBaseResetSelf = 0;
			std::size_t mapmulBaseResetOther = 0;
			std::size_t mapmulLookupHits = 0;
			std::size_t mapmulLookupPhaseful = 0;
			std::size_t mapmulResultPhaseful = 0;
			std::size_t mapdivBaseResetSelf = 0;
			std::size_t mapdivHeaderReset = 0;
			std::size_t mapdivLookupHits = 0;
			std::size_t mapdivLookupPhaseful = 0;
			std::size_t mapdivLookupPhaseOverwrite = 0;
			std::size_t mapdivLookupPhaseOverwriteNonHeader = 0;
			std::size_t mapdivResultPhaseful = 0;
			std::size_t findRemainPhaseCarries = 0;
		};

		///
		/// Construction, destruction, information and reset
		///

		static constexpr std::size_t MAX_POSSIBLE_QUBITS = static_cast<std::make_unsigned_t<Qubit>>(std::numeric_limits<Qubit>::max()) + 1U;
		static constexpr std::size_t DEFAULT_QUBITS = 300;


		//==========================================我写的========================================
		bool to_test = false;
		bool enableRegressionDiagnostics = false;
		bool disableMapdivLookupWriteback = false;
		bool enableTailCxRenormExperiment = false;
		bool enableUniqueSemanticProbe = false;
		RegressionDiagnostics regressionDiagnostics{};
		bool enableContStageTrace = false;
		std::size_t contStageTraceStep = 0;
		std::size_t contStageTraceDepth = 0;
		bool contStageTraceNextChild = false;
		bool contStageTraceFocusedChildActive = false;
		std::size_t contStageTraceNextChildParentDepth = 0;
		const char* contStageTraceNextChildBranch = nullptr;
		int contStageTraceNextChildK = -1;
		std::size_t contStageTaddDepth = 0;
		std::size_t contStageMakeNodeDepth = 0;
		RegressionDiagnostics contStageTaddTotals{};
		RegressionDiagnostics contStageMakeNodeTotals{};

		int mode = 1;//设置提取的对角门的形式，mode=1,提取的只是Rz旋转门，mode=2,提取的是任意对角门；

		std::map<std::string, int> varOrder;

		Edge<mNode> identity;
		//==========================================我写的========================================
		explicit Package(std::size_t nq = DEFAULT_QUBITS) : nqubits(nq) {
			resize(nq);
			this->identity = this->xarray_2_edge({{{1,0},{0,0}},{{0,0},{1,0}}},{0,1});
		};
		~Package() = default;
		Package(const Package& package) = delete;

		Package& operator=(const Package& package) = delete;

		// resize the package instance
		void resize(std::size_t nq) {
			if (nq > MAX_POSSIBLE_QUBITS) {
				throw std::invalid_argument("Requested too many qubits from package. "
					"Qubit datatype only allows up to " +
					std::to_string(MAX_POSSIBLE_QUBITS) +
					" qubits, while " + std::to_string(nq) +
					" were requested. Please recompile the "
					"package with a wider Qubit type!");
			}
			nqubits = nq;
			nodeUniqueTable.resize(nqubits);
		}

		// reset package state
		void reset() {
			clearUniqueTables();
			clearComputeTables();
			cn.clear();
		}

		// getter for qubits
		[[nodiscard]] auto qubits() const { return nqubits; }

		void resetRegressionDiagnostics() {
			regressionDiagnostics = {};
		}

		void setContStageTrace(const bool enabled, const std::size_t step = 0) {
			enableContStageTrace = enabled;
			contStageTraceStep = step;
			contStageTaddTotals = {};
			contStageMakeNodeTotals = {};
			contStageTaddDepth = 0;
			contStageMakeNodeDepth = 0;
			contStageTraceNextChild = false;
			contStageTraceFocusedChildActive = false;
			contStageTraceNextChildParentDepth = 0;
			contStageTraceNextChildBranch = nullptr;
			contStageTraceNextChildK = -1;
			if (!enabled) {
				contStageTraceDepth = 0;
			}
		}

		static void accumulateRegressionDelta(RegressionDiagnostics& total, const RegressionDiagnostics& before, const RegressionDiagnostics& after) {
			total.normalizeZeroChildren += after.normalizeZeroChildren - before.normalizeZeroChildren;
			total.normalizeChildPhaseAdds += after.normalizeChildPhaseAdds - before.normalizeChildPhaseAdds;
			total.normalizeRootPhasePromotions += after.normalizeRootPhasePromotions - before.normalizeRootPhasePromotions;
			total.normalizeChild0FromInput0 += after.normalizeChild0FromInput0 - before.normalizeChild0FromInput0;
			total.normalizeChild0FromInput1 += after.normalizeChild0FromInput1 - before.normalizeChild0FromInput1;
			total.normalizeChild1FromInput0 += after.normalizeChild1FromInput0 - before.normalizeChild1FromInput0;
			total.normalizeChild1FromInput1 += after.normalizeChild1FromInput1 - before.normalizeChild1FromInput1;
			total.normalizeChild1Zeroed += after.normalizeChild1Zeroed - before.normalizeChild1Zeroed;
			total.normalizeChild1SameWeightDifferentMap += after.normalizeChild1SameWeightDifferentMap - before.normalizeChild1SameWeightDifferentMap;
			total.normalizeChild1SameMapDifferentWeight += after.normalizeChild1SameMapDifferentWeight - before.normalizeChild1SameMapDifferentWeight;
			total.normalizeChild1DifferentMapAndWeight += after.normalizeChild1DifferentMapAndWeight - before.normalizeChild1DifferentMapAndWeight;
			total.normalizeChild1DiffMapAndWeightSameNode += after.normalizeChild1DiffMapAndWeightSameNode - before.normalizeChild1DiffMapAndWeightSameNode;
			total.normalizeChild1DiffMapAndWeightDifferentNode += after.normalizeChild1DiffMapAndWeightDifferentNode - before.normalizeChild1DiffMapAndWeightDifferentNode;
			total.normalizeChild1DiffMapAndWeightAngleSnapped += after.normalizeChild1DiffMapAndWeightAngleSnapped - before.normalizeChild1DiffMapAndWeightAngleSnapped;
			total.normalizeChild1DiffMapAndWeightResidualWeight += after.normalizeChild1DiffMapAndWeightResidualWeight - before.normalizeChild1DiffMapAndWeightResidualWeight;
			total.normalizeChild1DiffMapAndWeightRotNonZero += after.normalizeChild1DiffMapAndWeightRotNonZero - before.normalizeChild1DiffMapAndWeightRotNonZero;
			total.normalizeChild1DiffMapAndWeightNearSnapBoundary += after.normalizeChild1DiffMapAndWeightNearSnapBoundary - before.normalizeChild1DiffMapAndWeightNearSnapBoundary;
			total.normalizeChild1DiffMapAndWeightSnapInterior += after.normalizeChild1DiffMapAndWeightSnapInterior - before.normalizeChild1DiffMapAndWeightSnapInterior;
			total.normalizeChild1DiffMapAndWeightSnapBoundaryInside += after.normalizeChild1DiffMapAndWeightSnapBoundaryInside - before.normalizeChild1DiffMapAndWeightSnapBoundaryInside;
			total.normalizeChild1DiffMapAndWeightResidualBoundaryOutside += after.normalizeChild1DiffMapAndWeightResidualBoundaryOutside - before.normalizeChild1DiffMapAndWeightResidualBoundaryOutside;
			total.normalizeChild1DiffMapAndWeightResidualFarOutside += after.normalizeChild1DiffMapAndWeightResidualFarOutside - before.normalizeChild1DiffMapAndWeightResidualFarOutside;
			total.normalizeChild1DiffMapAndWeightRotPositive += after.normalizeChild1DiffMapAndWeightRotPositive - before.normalizeChild1DiffMapAndWeightRotPositive;
			total.normalizeChild1DiffMapAndWeightRotNegative += after.normalizeChild1DiffMapAndWeightRotNegative - before.normalizeChild1DiffMapAndWeightRotNegative;
			total.normalizeChild1DiffMapAndWeightRotZero += after.normalizeChild1DiffMapAndWeightRotZero - before.normalizeChild1DiffMapAndWeightRotZero;
			total.normalizeChild1LegacyAbsSnapDisagreement += after.normalizeChild1LegacyAbsSnapDisagreement - before.normalizeChild1LegacyAbsSnapDisagreement;
			total.normalizeChild1PhasefulAfterMapdiv += after.normalizeChild1PhasefulAfterMapdiv - before.normalizeChild1PhasefulAfterMapdiv;
			total.normalizeRootPhasePromotionsFromInput0 += after.normalizeRootPhasePromotionsFromInput0 - before.normalizeRootPhasePromotionsFromInput0;
			total.normalizeRootPhasePromotionsFromInput1 += after.normalizeRootPhasePromotionsFromInput1 - before.normalizeRootPhasePromotionsFromInput1;
			total.uniqueSemanticAlternatives += after.uniqueSemanticAlternatives - before.uniqueSemanticAlternatives;
			total.uniqueSemanticAlternativeMisses += after.uniqueSemanticAlternativeMisses - before.uniqueSemanticAlternativeMisses;
			total.uniqueSemanticGlobalAlternatives += after.uniqueSemanticGlobalAlternatives - before.uniqueSemanticGlobalAlternatives;
			total.uniqueSemanticGlobalAlternativeMisses += after.uniqueSemanticGlobalAlternativeMisses - before.uniqueSemanticGlobalAlternativeMisses;
			total.uniqueStructuralAlternatives += after.uniqueStructuralAlternatives - before.uniqueStructuralAlternatives;
			total.uniqueStructuralAlternativeMisses += after.uniqueStructuralAlternativeMisses - before.uniqueStructuralAlternativeMisses;
			total.uniqueStructuralSameWeightDifferentMap += after.uniqueStructuralSameWeightDifferentMap - before.uniqueStructuralSameWeightDifferentMap;
			total.uniqueStructuralSameMapDifferentWeight += after.uniqueStructuralSameMapDifferentWeight - before.uniqueStructuralSameMapDifferentWeight;
			total.uniqueStructuralDifferentMapAndWeight += after.uniqueStructuralDifferentMapAndWeight - before.uniqueStructuralDifferentMapAndWeight;
			total.uniquePreNormalizeStructuralAlternatives += after.uniquePreNormalizeStructuralAlternatives - before.uniquePreNormalizeStructuralAlternatives;
			total.uniquePreNormalizeSameWeightDifferentMap += after.uniquePreNormalizeSameWeightDifferentMap - before.uniquePreNormalizeSameWeightDifferentMap;
			total.uniquePreNormalizeSameMapDifferentWeight += after.uniquePreNormalizeSameMapDifferentWeight - before.uniquePreNormalizeSameMapDifferentWeight;
			total.uniquePreNormalizeDifferentMapAndWeight += after.uniquePreNormalizeDifferentMapAndWeight - before.uniquePreNormalizeDifferentMapAndWeight;
			total.uniqueNormalizeTransitionDiffBothToMapOnly += after.uniqueNormalizeTransitionDiffBothToMapOnly - before.uniqueNormalizeTransitionDiffBothToMapOnly;
			total.uniqueNormalizeTransitionDiffBothToWeightOnly += after.uniqueNormalizeTransitionDiffBothToWeightOnly - before.uniqueNormalizeTransitionDiffBothToWeightOnly;
			total.uniqueNormalizeTransitionDiffBothToDiffBoth += after.uniqueNormalizeTransitionDiffBothToDiffBoth - before.uniqueNormalizeTransitionDiffBothToDiffBoth;
			total.taddSamePointerMapMismatch += after.taddSamePointerMapMismatch - before.taddSamePointerMapMismatch;
			total.taddMismatchResidualHeader += after.taddMismatchResidualHeader - before.taddMismatchResidualHeader;
			total.taddMismatchResidualNonHeader += after.taddMismatchResidualNonHeader - before.taddMismatchResidualNonHeader;
			total.taddMismatchResidualPhaseful += after.taddMismatchResidualPhaseful - before.taddMismatchResidualPhaseful;
			total.taddMismatchAddHits += after.taddMismatchAddHits - before.taddMismatchAddHits;
			total.taddMismatchAddMisses += after.taddMismatchAddMisses - before.taddMismatchAddMisses;
			total.taddMismatchAddMissEmpty += after.taddMismatchAddMissEmpty - before.taddMismatchAddMissEmpty;
			total.taddMismatchAddMissMapOnly += after.taddMismatchAddMissMapOnly - before.taddMismatchAddMissMapOnly;
			total.taddMismatchAddMissWeightOnly += after.taddMismatchAddMissWeightOnly - before.taddMismatchAddMissWeightOnly;
			total.taddMismatchAddMissMapAndWeight += after.taddMismatchAddMissMapAndWeight - before.taddMismatchAddMissMapAndWeight;
			total.mapmulBaseResetSelf += after.mapmulBaseResetSelf - before.mapmulBaseResetSelf;
			total.mapmulBaseResetOther += after.mapmulBaseResetOther - before.mapmulBaseResetOther;
			total.mapmulLookupHits += after.mapmulLookupHits - before.mapmulLookupHits;
			total.mapmulLookupPhaseful += after.mapmulLookupPhaseful - before.mapmulLookupPhaseful;
			total.mapmulResultPhaseful += after.mapmulResultPhaseful - before.mapmulResultPhaseful;
			total.mapdivBaseResetSelf += after.mapdivBaseResetSelf - before.mapdivBaseResetSelf;
			total.mapdivHeaderReset += after.mapdivHeaderReset - before.mapdivHeaderReset;
			total.mapdivLookupHits += after.mapdivLookupHits - before.mapdivLookupHits;
			total.mapdivLookupPhaseful += after.mapdivLookupPhaseful - before.mapdivLookupPhaseful;
			total.mapdivLookupPhaseOverwrite += after.mapdivLookupPhaseOverwrite - before.mapdivLookupPhaseOverwrite;
			total.mapdivLookupPhaseOverwriteNonHeader += after.mapdivLookupPhaseOverwriteNonHeader - before.mapdivLookupPhaseOverwriteNonHeader;
			total.mapdivResultPhaseful += after.mapdivResultPhaseful - before.mapdivResultPhaseful;
			total.findRemainPhaseCarries += after.findRemainPhaseCarries - before.findRemainPhaseCarries;
		}

		std::ostream& printRegressionDiagnostics(std::ostream& os = std::cout) {
			os << "regression.normalize.zero_children=" << regressionDiagnostics.normalizeZeroChildren << std::endl;
			os << "regression.normalize.child_phase_adds=" << regressionDiagnostics.normalizeChildPhaseAdds << std::endl;
			os << "regression.normalize.root_phase_promotions=" << regressionDiagnostics.normalizeRootPhasePromotions << std::endl;
			os << "regression.normalize.child0_from_input0=" << regressionDiagnostics.normalizeChild0FromInput0 << std::endl;
			os << "regression.normalize.child0_from_input1=" << regressionDiagnostics.normalizeChild0FromInput1 << std::endl;
			os << "regression.normalize.child1_from_input0=" << regressionDiagnostics.normalizeChild1FromInput0 << std::endl;
			os << "regression.normalize.child1_from_input1=" << regressionDiagnostics.normalizeChild1FromInput1 << std::endl;
			os << "regression.normalize.child1_zeroed=" << regressionDiagnostics.normalizeChild1Zeroed << std::endl;
			os << "regression.normalize.child1_same_weight_diff_map=" << regressionDiagnostics.normalizeChild1SameWeightDifferentMap << std::endl;
			os << "regression.normalize.child1_same_map_diff_weight=" << regressionDiagnostics.normalizeChild1SameMapDifferentWeight << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight=" << regressionDiagnostics.normalizeChild1DifferentMapAndWeight << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_same_node=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightSameNode << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_different_node=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightDifferentNode << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_angle_snapped=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightAngleSnapped << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_residual_weight=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightResidualWeight << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_rot_nonzero=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightRotNonZero << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_near_snap_boundary=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightNearSnapBoundary << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_snap_interior=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightSnapInterior << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_snap_boundary_inside=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightSnapBoundaryInside << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_residual_boundary_outside=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightResidualBoundaryOutside << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_residual_far_outside=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightResidualFarOutside << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_rot_positive=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightRotPositive << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_rot_negative=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightRotNegative << std::endl;
			os << "regression.normalize.child1_diff_map_and_weight_rot_zero=" << regressionDiagnostics.normalizeChild1DiffMapAndWeightRotZero << std::endl;
			os << "regression.normalize.child1_legacy_abs_snap_disagreement=" << regressionDiagnostics.normalizeChild1LegacyAbsSnapDisagreement << std::endl;
			os << "regression.normalize.child1_phaseful_after_mapdiv=" << regressionDiagnostics.normalizeChild1PhasefulAfterMapdiv << std::endl;
			os << "regression.normalize.root_phase_promotions_from_input0=" << regressionDiagnostics.normalizeRootPhasePromotionsFromInput0 << std::endl;
			os << "regression.normalize.root_phase_promotions_from_input1=" << regressionDiagnostics.normalizeRootPhasePromotionsFromInput1 << std::endl;
			os << "regression.unique.semantic_alternatives=" << regressionDiagnostics.uniqueSemanticAlternatives << std::endl;
			os << "regression.unique.semantic_alternative_misses=" << regressionDiagnostics.uniqueSemanticAlternativeMisses << std::endl;
			os << "regression.unique.semantic_global_alternatives=" << regressionDiagnostics.uniqueSemanticGlobalAlternatives << std::endl;
			os << "regression.unique.semantic_global_alternative_misses=" << regressionDiagnostics.uniqueSemanticGlobalAlternativeMisses << std::endl;
			os << "regression.unique.structural_alternatives=" << regressionDiagnostics.uniqueStructuralAlternatives << std::endl;
			os << "regression.unique.structural_alternative_misses=" << regressionDiagnostics.uniqueStructuralAlternativeMisses << std::endl;
			os << "regression.unique.structural_same_weight_diff_map=" << regressionDiagnostics.uniqueStructuralSameWeightDifferentMap << std::endl;
			os << "regression.unique.structural_same_map_diff_weight=" << regressionDiagnostics.uniqueStructuralSameMapDifferentWeight << std::endl;
			os << "regression.unique.structural_diff_map_and_weight=" << regressionDiagnostics.uniqueStructuralDifferentMapAndWeight << std::endl;
			os << "regression.unique.pre_normalize_structural_alternatives=" << regressionDiagnostics.uniquePreNormalizeStructuralAlternatives << std::endl;
			os << "regression.unique.pre_normalize_same_weight_diff_map=" << regressionDiagnostics.uniquePreNormalizeSameWeightDifferentMap << std::endl;
			os << "regression.unique.pre_normalize_same_map_diff_weight=" << regressionDiagnostics.uniquePreNormalizeSameMapDifferentWeight << std::endl;
			os << "regression.unique.pre_normalize_diff_map_and_weight=" << regressionDiagnostics.uniquePreNormalizeDifferentMapAndWeight << std::endl;
			os << "regression.unique.normalize_transition_diff_both_to_map_only=" << regressionDiagnostics.uniqueNormalizeTransitionDiffBothToMapOnly << std::endl;
			os << "regression.unique.normalize_transition_diff_both_to_weight_only=" << regressionDiagnostics.uniqueNormalizeTransitionDiffBothToWeightOnly << std::endl;
			os << "regression.unique.normalize_transition_diff_both_to_diff_both=" << regressionDiagnostics.uniqueNormalizeTransitionDiffBothToDiffBoth << std::endl;
			os << "regression.tadd.same_pointer_map_mismatch=" << regressionDiagnostics.taddSamePointerMapMismatch << std::endl;
			os << "regression.tadd.mismatch_residual_header=" << regressionDiagnostics.taddMismatchResidualHeader << std::endl;
			os << "regression.tadd.mismatch_residual_non_header=" << regressionDiagnostics.taddMismatchResidualNonHeader << std::endl;
			os << "regression.tadd.mismatch_residual_phaseful=" << regressionDiagnostics.taddMismatchResidualPhaseful << std::endl;
			os << "regression.tadd.mismatch_add_hits=" << regressionDiagnostics.taddMismatchAddHits << std::endl;
			os << "regression.tadd.mismatch_add_misses=" << regressionDiagnostics.taddMismatchAddMisses << std::endl;
			os << "regression.tadd.mismatch_add_miss_empty=" << regressionDiagnostics.taddMismatchAddMissEmpty << std::endl;
			os << "regression.tadd.mismatch_add_miss_map_only=" << regressionDiagnostics.taddMismatchAddMissMapOnly << std::endl;
			os << "regression.tadd.mismatch_add_miss_weight_only=" << regressionDiagnostics.taddMismatchAddMissWeightOnly << std::endl;
			os << "regression.tadd.mismatch_add_miss_map_and_weight=" << regressionDiagnostics.taddMismatchAddMissMapAndWeight << std::endl;
			os << "regression.mapmul.base_reset_self=" << regressionDiagnostics.mapmulBaseResetSelf << std::endl;
			os << "regression.mapmul.base_reset_other=" << regressionDiagnostics.mapmulBaseResetOther << std::endl;
			os << "regression.mapmul.lookup_hits=" << regressionDiagnostics.mapmulLookupHits << std::endl;
			os << "regression.mapmul.lookup_phaseful=" << regressionDiagnostics.mapmulLookupPhaseful << std::endl;
			os << "regression.mapmul.result_phaseful=" << regressionDiagnostics.mapmulResultPhaseful << std::endl;
			os << "regression.mapdiv.base_reset_self=" << regressionDiagnostics.mapdivBaseResetSelf << std::endl;
			os << "regression.mapdiv.header_reset=" << regressionDiagnostics.mapdivHeaderReset << std::endl;
			os << "regression.mapdiv.lookup_hits=" << regressionDiagnostics.mapdivLookupHits << std::endl;
			os << "regression.mapdiv.lookup_phaseful=" << regressionDiagnostics.mapdivLookupPhaseful << std::endl;
			os << "regression.mapdiv.lookup_phase_overwrite=" << regressionDiagnostics.mapdivLookupPhaseOverwrite << std::endl;
			os << "regression.mapdiv.lookup_phase_overwrite_non_header=" << regressionDiagnostics.mapdivLookupPhaseOverwriteNonHeader << std::endl;
			os << "regression.mapdiv.result_phaseful=" << regressionDiagnostics.mapdivResultPhaseful << std::endl;
			os << "regression.find_remain.phase_carries=" << regressionDiagnostics.findRemainPhaseCarries << std::endl;
			return os;
		}

	private:
		std::size_t nqubits;

		static bool inline ifContract(float k){
			//when k is half, means this index will be contracted
			return std::abs(k - std::round(k)) > 0.4999999 && std::abs(k - std::round(k)) < 0.5000001;
		}

		///
		/// Vector nodes, edges and quantum states
		///
	public:


		//==========================================我写的========================================
		Edge<mNode> xarray_2_edge(
			const xt::xarray<ComplexValue>& array,
			std::vector<int> order){
						if (array.size() == 1) {
				 if (array[0].approximatelyZero()) {
				 	return Edge<mNode>::zero;
				 }
				 else if (array[0].approximatelyOne()) {
				 	return Edge<mNode>::one;
				 }
				 else {
				 	return Edge<mNode>::terminal(cn.lookup(array[0]));
				 }
			}

			auto split_pos = std::distance(order.begin(), std::max_element(order.begin(), order.end()));
			int16_t x = order[split_pos];
			order[split_pos] = -1;
			
			std::vector<xt::xarray<ComplexValue>> split_U;
			for (const auto& u : xt::split(array, array.shape(split_pos), split_pos)) {
				split_U.push_back(u);
			}
			// xt::split return a const type data, it is not good.

			while (std::find(order.begin(), order.end(), x) != order.end()) {
				auto split_pos = std::distance(order.begin(), std::find(order.begin(), order.end(), x));
				order[split_pos] = -1;
				std::vector<xt::xarray<ComplexValue>> temp_U;
				for (int i = 0; i < split_U.size(); ++i) {
					auto temp = xt::split(split_U.at(i), split_U.at(i).shape(split_pos), split_pos).at(i);
					temp_U.push_back(temp);
				}
				split_U = temp_U;
			}

			std::vector<Edge<mNode>> edges;
			for (const auto u : split_U) {
				edges.push_back(xarray_2_edge(u, order));
			}

			return makeDDNode((int16_t)x, edges, false);

		}

		int update_map_value() {
			if (varOrder.empty()) {
				return 0;
			}
			auto compare_function = [](const auto& a, const auto& b) { return a.second < b.second; };
			auto max_pair = std::max_element(varOrder.begin(), varOrder.end(), compare_function);
			return (max_pair->second) + 1;
		}


		//==========================================我写的========================================

		///
		/// Matrix nodes, edges and quantum gates
		///
		template <class Node> Edge<Node> normalize(const Edge<Node>& e, bool cached) {
			const auto originalEdges = e.p->e;
			const auto traceNanProbe = enableContStageTrace && contStageTraceDepth > 0;
			auto formatComplex = [](const Complex& value) {
				std::ostringstream output;
				output << "(" << CTEntry::val(value.r) << ", " << CTEntry::val(value.i) << ")";
				return output.str();
			};
			auto throwOnNonFinite = [&](const char* stage,
										 const std::size_t edgeIndex,
										 const Complex& current,
										 const Complex& currentMaxValue,
										 const fp angle,
										 const fp deltaAngle,
										 const fp magnitude,
										 const int rot,
										 const the_maps* currentMap,
										 const the_maps* currentResMap) {
				if (!traceNanProbe) {
					return;
				}
				const auto real = CTEntry::val(current.r);
				const auto imag = CTEntry::val(current.i);
				const auto maxReal = CTEntry::val(currentMaxValue.r);
				const auto maxImag = CTEntry::val(currentMaxValue.i);
				if (std::isfinite(real) && std::isfinite(imag) && std::isfinite(maxReal) && std::isfinite(maxImag) &&
					std::isfinite(angle) && std::isfinite(deltaAngle) && std::isfinite(magnitude)) {
					return;
				}
				std::ostringstream message;
				message << "normalize non-finite"
					<< " step=" << contStageTraceStep
					<< " stage=" << stage
					<< " var=" << static_cast<int>(e.p->v)
					<< " edge_index=" << edgeIndex
					<< " current=" << formatComplex(current)
					<< " max_value=" << formatComplex(currentMaxValue)
					<< " angle=" << angle
					<< " delta_angle=" << deltaAngle
					<< " magnitude=" << magnitude
					<< " rot=" << rot
					<< " edge_map_extra_phase=" << (currentMap ? currentMap->extra_phase : -999999)
					<< " res_map_extra_phase=" << (currentResMap ? currentResMap->extra_phase : -999999);
				throw std::runtime_error(message.str());
			};

			auto maxArgIndex = -1;
			// v0 = e.p->e[0].p
			// v1 = e.p->e[1].p
			// w0 = e.p->e[0].w
			// w0 = e.p->e[1].w
			//auto zero = std::array{	e.p->e[0].w.approximatelyZero(), e.p->e[1].w.approximatelyZero()};
			int nodeCount = e.p->e.size();
			// during now, nodeCount should be  2
			// if(to_test){
			// 	std::cout << "224: " << nodeCount <<std::endl;
			// }
			assert(nodeCount == 2);

			std::vector<bool> isZero(nodeCount, false);

			// Check if any values are approximately isZero
			for (int k = 0; k < nodeCount; k++) {
				isZero[k] = e.p->e[k].w.approximatelyZero();
			}

			// Release cached numbers approximately zero, but not exactly zero
			if(cached) {
				for (auto i = 0U; i < nodeCount; i++) {
					if (isZero[i] && e.p->e[i].w != Complex::zero) {
						// cn.returnToCache(e.p->e[i].w);
						e.p->e[i] = Edge<Node>::zero;
					}
				}
			}

			fp max_mag2 = 0;
			auto max_value = Complex::one;
			// determine max amplitude
			for (auto i = 0U; i < nodeCount; ++i) {
				//std::cout << 299 << " " << isZero[i]<<" "<< e.p->e[i].w << std::endl;
				if (isZero[i]) {
					continue;
				}
				if (maxArgIndex == -1) {
					maxArgIndex = static_cast<decltype(maxArgIndex)>(i);
					max_mag2 = ComplexNumbers::mag2(e.p->e[i].w);
					throwOnNonFinite("max_init", i, e.p->e[i].w, e.p->e[i].w, 0.0, 0.0, max_mag2, 0, e.p->e[i].map, e.map);
					max_value = e.p->e[i].w;
				}
				else {
					auto mag = ComplexNumbers::mag2(e.p->e[i].w);
					throwOnNonFinite("max_compare", i, e.p->e[i].w, max_value, 0.0, 0.0, mag, 0, e.p->e[i].map, e.map);
					if (mag - max_mag2 > ComplexTable<>::tolerance()/2) {
						maxArgIndex = static_cast<decltype(maxArgIndex)>(i);
						max_mag2 = mag;
						max_value = e.p->e[i].w;
					}
				}
				//std::cout << 315 << " " << i << " " << maxArgIndex << " " << max << " " << max_value << std::endl;
			}

			// all equal to zero
			if (maxArgIndex == -1) {
				if (!cached && !e.isTerminal()) {
					// If it is not a cached computation, the node has to be put back into
					// the chain
					getUniqueTable<Node>().returnNode(e.p);
				}
				return Edge<Node>::zero;
			}

			bool add_x = 0;
			auto res = e;

			if (maxArgIndex > 0) {
				add_x = true;
			}else{
				add_x = false;
			}
			// else if(ComplexNumbers::mag2(res.p->e[0].w)-ComplexNumbers::mag2(res.p->e[1].w) > ComplexTable<>::tolerance()){
			// 	add_x = false;
			// }
			// else if(res.p->e[0].p >  res.p->e[1].p){
			// 	add_x = true ;
			// }
			// else if(res.p->e[0].p <  res.p->e[1].p){
			// 	add_x = false;
			// }
			// else{
			// 	add_x = mapCompare(res.p->e[0].map,res.p->e[1].map);
			// }
				// if( && ComplexNumbers::arg(res.p->e[0].w) >  ComplexNumbers::arg(res.p->e[1].w)){
				// 	res.p->e = {res.p->e[1],res.p->e[0]};
				// 	isZero = { isZero[1],isZero[0] };
				// 	add_x = 1;
				// }

			if(add_x){
				res.p->e = {res.p->e[1],res.p->e[0]};
				isZero = { isZero[1],isZero[0] };
			}
			const auto child0SourceIndex = add_x ? 1U : 0U;
			const auto child1SourceIndex = add_x ? 0U : 1U;
			const auto& child1SourceEdge = originalEdges[child1SourceIndex];
			const auto child1SourceWasZero = child1SourceEdge.w.approximatelyZero();
			bool child1AngleSnapped = false;
			bool child1ResidualWeight = false;
			bool child1RotNonZero = false;
			bool child1NearSnapBoundary = false;
			bool child1SnapInterior = false;
			bool child1SnapBoundaryInside = false;
			bool child1ResidualBoundaryOutside = false;
			bool child1ResidualFarOutside = false;
			int child1RotSign = 0;
			bool child1LegacyAbsSnapDisagreement = false;
			if (enableRegressionDiagnostics) {
				if (child0SourceIndex == 0U) {
					regressionDiagnostics.normalizeChild0FromInput0++;
					regressionDiagnostics.normalizeChild1FromInput1++;
				} else {
					regressionDiagnostics.normalizeChild0FromInput1++;
					regressionDiagnostics.normalizeChild1FromInput0++;
				}
			}
			maxArgIndex = 0;

			//std::cout << "aaa" << std::endl;
			//std::cout << r.p->e[0].w << " " << r.p->e[1].w << std::endl;
			//the_maps::print_maps(r.p->e[0].map);
			//the_maps::print_maps(r.p->e[1].map);

			//std::cout << argmax<<" " << max_value << std::endl;
			// divide each entry by max
			for (auto i = 0U; i < nodeCount; ++i) {
				if (static_cast<decltype(maxArgIndex)>(i) == maxArgIndex) {
					if (cached) {
						if (res.w.exactlyOne()) {
							res.w = max_value;
						}
						else {
							assert(res.w != Complex::zero);
							ComplexNumbers::mul(res.w, res.w, max_value);
						}
					}
					else {
						if (res.w.exactlyOne()) {
							res.w = max_value;
						}
						else {
							auto c = cn.getTemporary();
							assert(c != Complex::zero);
							ComplexNumbers::mul(c, res.w, max_value);
							res.w = cn.lookup(c);
						}
					}
					res.map = res.p->e[i].map;
					res.p->e[i].w = Complex::one;
					res.p->e[i].map = the_maps::the_maps_header();
				}
				else {
					if (isZero[i]) {
						if (cached && res.p->e[i].w != Complex::zero) {
							//assert(r.p->e[i].w != Complex::zero);
							// cn.returnToCache(res.p->e[i].w);
						}
						//r.p->e[i] = Edge<Node>::zero;
						res.p->e[i] = { res.p->e[0].p,Complex::zero, the_maps::the_maps_header() };
						if (enableRegressionDiagnostics) {
							regressionDiagnostics.normalizeZeroChildren++;
							regressionDiagnostics.normalizeChild1Zeroed++;
						}
						res.p->e[i].map->extra_phase = 0;
						continue;
					}
					if (cached && !isZero[i] && !res.p->e[i].w.exactlyOne()) {
						//assert(r.p->e[i].w != Complex::zero);
						// cn.returnToCache(res.p->e[i].w);
					}
					if (res.p->e[i].w.approximatelyOne()) {
						res.p->e[i].w = Complex::one;
					}

					if (mode == 2) {
						auto c = cn.getCached();
						ComplexNumbers::div(c, res.p->e[i].w, max_value);
						res.p->e[i].map = mapdiv(res.p->e[i].map, res.map);
						// cn.mul(res.p->e[i].map->extra_phase, res.p->e[i].map->extra_phase, c);
						// cn.mul(res.p->e[i].map->extra_phase, res.p->e[i].map->extra_phase, c);
						if (enableRegressionDiagnostics) {
							regressionDiagnostics.normalizeChildPhaseAdds++;
						}
						res.p->e[i].map->extra_phase=res.p->e[i].map->extra_phase+int(ComplexNumbers::arg(c)/rotate_angle);
						// cn.returnToCache(c);
						res.p->e[i].w = Complex::one;

					}
					else {
						auto c = cn.getTemporary();
						ComplexNumbers::div(c, res.p->e[i].w, max_value);
						auto angle = ComplexNumbers::arg(c);
						int rot = round(angle / rotate_angle);
						double detla_angle = angle - rot * rotate_angle;
						throwOnNonFinite("post_div", i, c, max_value, angle, detla_angle, ComplexNumbers::mag2(c), rot, res.p->e[i].map, res.map);
						if (i == 1U) {
							child1RotNonZero = (rot != 0);
							const auto snapThreshold = ComplexTable<>::tolerance() * rotate_angle;
							const auto absDeltaAngle = std::abs(detla_angle);
							const auto boundaryDistance = std::abs(absDeltaAngle - snapThreshold);
							child1NearSnapBoundary = boundaryDistance <= snapThreshold;
							child1SnapInterior = absDeltaAngle <= snapThreshold * 0.5;
							child1SnapBoundaryInside = absDeltaAngle > snapThreshold * 0.5 && absDeltaAngle < snapThreshold;
							child1ResidualBoundaryOutside = absDeltaAngle >= snapThreshold && absDeltaAngle <= snapThreshold * 1.5;
							child1ResidualFarOutside = absDeltaAngle > snapThreshold * 1.5;
							child1RotSign = (rot > 0) ? 1 : ((rot < 0) ? -1 : 0);
							const auto legacyAbsSnapped = abs(detla_angle) < ComplexTable<>::tolerance() * rotate_angle;
							const auto explicitAbsSnapped = absDeltaAngle < snapThreshold;
							child1LegacyAbsSnapDisagreement = legacyAbsSnapped != explicitAbsSnapped;
						}
						if (std::abs(detla_angle) < ComplexTable<>::tolerance()* rotate_angle) {
							if (i == 1U) {
								child1AngleSnapped = true;
							}
							if (std::getenv("LIMTDD_TRACE_NORMALIZE_SNAP") != nullptr && enableContStageTrace && contStageTraceDepth > 0) {
								std::cerr << "normalize_snap"
									<< "\tstep\t" << contStageTraceStep
									<< "\tdepth\t" << contStageTraceDepth
									<< "\tvar\t" << static_cast<int>(e.p->v)
									<< "\tedge\t" << i
									<< "\tc_re_before\t" << CTEntry::val(c.r)
									<< "\tc_im_before\t" << CTEntry::val(c.i)
									<< "\tmag2\t" << ComplexNumbers::mag2(c)
									<< "\trot\t" << rot
									<< "\tdelta\t" << detla_angle
									<< "\n";
							}
							c.r->value = sqrt(ComplexNumbers::mag2(c));
							c.i->value = 0;
							throwOnNonFinite("snap_lookup", i, c, max_value, angle, detla_angle, ComplexNumbers::mag2(c), rot, res.p->e[i].map, res.map);
							//std::cout << c << " a " << ComplexNumbers::mag2(c) << std::endl;
							res.p->e[i].w = cn.lookup(c);
						}
						else {
							if (i == 1U) {
								child1ResidualWeight = true;
							}
							//c.r->value = sqrt(ComplexNumbers::mag2(c))*cos(angle- rot * rotate_angle);
							//c.i->value = sqrt(ComplexNumbers::mag2(c))*sin(angle - rot * rotate_angle);
							double mags = sqrt(ComplexNumbers::mag2(c));
							c.r->value = mags * cos(detla_angle);
							c.i->value = mags * sin(detla_angle);
							throwOnNonFinite("residual_lookup", i, c, max_value, angle, detla_angle, mags, rot, res.p->e[i].map, res.map);
							//std::cout << mags * cos(angle - rot * rotate_angle) * mags * cos(angle - rot * rotate_angle) + mags * sin(angle - rot * rotate_angle) * mags * sin(angle - rot * rotate_angle) << std::endl;
							//std::cout << c.r->value * c.r->value + c.i->value * c.i->value << std::endl;
							res.p->e[i].w = cn.lookup(c);
						}

						res.p->e[i].map = mapdiv(res.p->e[i].map, res.map);
						if (enableRegressionDiagnostics) {
							regressionDiagnostics.normalizeChildPhaseAdds++;
						}
						res.p->e[i].map->extra_phase = res.p->e[i].map->extra_phase + rot;
						
						// cn.mul(res.p->e[i].map->extra_phase, res.p->e[i].map->extra_phase, cn.getTemporary(cos(angle), sin(angle)));
						
						//std::cout << angle << " " << rotate_angle << " " << angle / rotate_angle << " " << round(angle / rotate_angle);
					}
				}
			}

			if (enableRegressionDiagnostics && !child1SourceWasZero && !res.p->e[1].w.approximatelyZero()) {
				const auto sameWeight = res.p->e[1].w.approximatelyEquals(child1SourceEdge.w);
				const auto sameMap = res.p->e[1].map == child1SourceEdge.map;
				if (sameWeight && !sameMap) {
					regressionDiagnostics.normalizeChild1SameWeightDifferentMap++;
				} else if (!sameWeight && sameMap) {
					regressionDiagnostics.normalizeChild1SameMapDifferentWeight++;
				} else if (!sameWeight && !sameMap) {
					regressionDiagnostics.normalizeChild1DifferentMapAndWeight++;
					if (res.p->e[1].p == child1SourceEdge.p) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightSameNode++;
					} else {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightDifferentNode++;
					}
					if (child1AngleSnapped) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightAngleSnapped++;
					}
					if (child1ResidualWeight) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightResidualWeight++;
					}
					if (child1RotNonZero) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightRotNonZero++;
					}
					if (child1NearSnapBoundary) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightNearSnapBoundary++;
					}
					if (child1SnapInterior) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightSnapInterior++;
					}
					if (child1SnapBoundaryInside) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightSnapBoundaryInside++;
					}
					if (child1ResidualBoundaryOutside) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightResidualBoundaryOutside++;
					}
					if (child1ResidualFarOutside) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightResidualFarOutside++;
					}
					if (child1RotSign > 0) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightRotPositive++;
					} else if (child1RotSign < 0) {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightRotNegative++;
					} else {
						regressionDiagnostics.normalizeChild1DiffMapAndWeightRotZero++;
					}
					if (child1LegacyAbsSnapDisagreement) {
						regressionDiagnostics.normalizeChild1LegacyAbsSnapDisagreement++;
					}
				}
				if (res.p->e[1].map->extra_phase != 0) {
					regressionDiagnostics.normalizeChild1PhasefulAfterMapdiv++;
				}
			}

			if (enableRegressionDiagnostics && res.p->e[1].map->extra_phase != 0) {
				regressionDiagnostics.normalizeRootPhasePromotions++;
				if (child1SourceIndex == 0U) {
					regressionDiagnostics.normalizeRootPhasePromotionsFromInput0++;
				} else {
					regressionDiagnostics.normalizeRootPhasePromotionsFromInput1++;
				}
			}
			res.map = append_new_map(res.map, res.p->v, add_x, res.p->e[1].map->extra_phase);
			if (!isZero[1]) {
				// cn.returnToCache(res.p->e[1].map->extra_phase);
			}
			//std::cout << r.w << std::endl;
			//the_maps::print_maps(r.map);
			//std::cout << "bbb" << std::endl;
			return res;

		}


	private:

		///
		/// Unique tables, Reference counting and garbage collection
		///
	public:
		// unique tables
		template <class Node> [[nodiscard]] auto& getUniqueTable() {
			return nodeUniqueTable;
		}

		template <class Node> void incRef(const Edge<Node>& e) {
			getUniqueTable<Node>().incRef(e);
		}
		template <class Node> void decRef(const Edge<Node>& e) {
			getUniqueTable<Node>().decRef(e);
		}

		UniqueTable<mNode, Config::UT_MAT_NBUCKET, Config::UT_MAT_INITIAL_ALLOCATION_SIZE>	nodeUniqueTable{nqubits};

		bool garbageCollect(bool force = false) {
			// return immediately if no table needs collection
			if (!force &&
				!nodeUniqueTable.possiblyNeedsCollection() &&
				!cn.complexTable.possiblyNeedsCollection()) {
				return false;
			}

			auto cCollect = cn.garbageCollect(force);
			if (cCollect > 0) {
				// Collecting garbage in the complex numbers table requires collecting the
				// node tables as well
				force = true;
			}

			auto mCollect = nodeUniqueTable.garbageCollect(force);

			// invalidate all compute tables where any component of the entry contains
			// numbers from the complex table if any complex numbers were collected
			if (mCollect > 0) {

				addTable.clear();
				contTable.clear();
			}
			return  mCollect > 0;
		}

		void clearUniqueTables() {
			nodeUniqueTable.clear();

		}

		// create a normalized DD node and return an edge pointing to it. The node is
		// not recreated if it already exists.
		template <class Node>
		Edge<Node> makeDDNode(
			Qubit var,
			const std::vector<Edge<Node>>& edges,
			bool cached = false) {
			const auto traceMakeNode = enableContStageTrace && contStageTraceDepth > 0;
			const auto traceBefore = traceMakeNode ? regressionDiagnostics : RegressionDiagnostics{};
			if (traceMakeNode) {
				contStageMakeNodeDepth++;
			}
			struct MakeNodeTraceGuard {
				Package* pkg;
				bool active;
				RegressionDiagnostics before;
				~MakeNodeTraceGuard() {
					if (!active) {
						return;
					}
					if (pkg->contStageMakeNodeDepth == 1) {
						Package::accumulateRegressionDelta(pkg->contStageMakeNodeTotals, before, pkg->regressionDiagnostics);
					}
					pkg->contStageMakeNodeDepth--;
				}
			} makeNodeTraceGuard{this, traceMakeNode, traceBefore};
				if(to_test){
					std::cout << "var: " << var << std::endl;
					std::cout << "edge.node.key: " << std::endl;
					for(auto edge:edges){
						std::cout << edge.p->v << " " ;
						// if(edge.p->v == var) throw std::runtime_error("bug here");
					}
					std::cout << std::endl;

				}

			auto& uniqueTable = getUniqueTable<Node>();
			Edge<Node> e{uniqueTable.getNode(), Complex::one};
			e.p->v = var;
			e.p->e = edges;

			assert(e.p->ref == 0);
			if (traceMakeNode) {
				for (std::size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
					const auto real = CTEntry::val(edges[edgeIndex].w.r);
					const auto imag = CTEntry::val(edges[edgeIndex].w.i);
					if (std::isfinite(real) && std::isfinite(imag)) {
						continue;
					}
					std::ostringstream message;
					message << "makeDDNode non-finite input"
						<< " step=" << contStageTraceStep
						<< " var=" << static_cast<int>(var)
						<< " edge_index=" << edgeIndex
						<< " weight=(" << real << ", " << imag << ")"
						<< " child_var=" << static_cast<int>(edges[edgeIndex].p->v)
						<< " map_extra_phase=" << (edges[edgeIndex].map ? edges[edgeIndex].map->extra_phase : -999999);
					throw std::runtime_error(message.str());
				}
			}
			const auto* focusedMakeNodeVarEnv = std::getenv("LIMTDD_FOCUSED_MAKENODE_VAR");
			const auto traceFocusedMakeNode = contStageTraceDepth > 1 &&
				(focusedMakeNodeVarEnv == nullptr || static_cast<int>(var) == std::stoi(focusedMakeNodeVarEnv));
			const auto traceThisMakeNode = traceMakeNode && (contStageTraceDepth == 1 || traceFocusedMakeNode);
			if (traceThisMakeNode) {
				std::cerr << "makeDDNode_trace_begin\tstep\t" << contStageTraceStep
					<< "\tvar\t" << static_cast<int>(var)
					<< "\tcached\t" << cached
					<< "\tedge_count\t" << edges.size()
					<< "\n";
				for (std::size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
					std::cerr << "makeDDNode_trace_input\tstep\t" << contStageTraceStep
						<< "\tvar\t" << static_cast<int>(var)
						<< "\tedge\t" << edgeIndex
						<< "\tw_re\t" << CTEntry::val(edges[edgeIndex].w.r)
						<< "\tw_im\t" << CTEntry::val(edges[edgeIndex].w.i)
						<< "\tchild_var\t" << static_cast<int>(edges[edgeIndex].p->v)
						<< "\tmap_level\t" << (edges[edgeIndex].map ? static_cast<int>(edges[edgeIndex].map->level) : -999)
						<< "\tmap_x\t" << (edges[edgeIndex].map ? static_cast<int>(edges[edgeIndex].map->x) : -1)
						<< "\tmap_rot\t" << (edges[edgeIndex].map ? edges[edgeIndex].map->rotate : -999)
						<< "\tmap_ep\t" << (edges[edgeIndex].map ? edges[edgeIndex].map->extra_phase : -999)
						<< "\n";
				}
			}


			//if (edges[0].p == edges[1].p && e.p->e[0].w.approximatelyEquals(e.p->e[1].w)) {
			//	if (cached) {
			//		if (e.p->e[1].w != Complex::zero) {
			//			cn.returnToCache(e.p->e[1].w);
			//			return edges[0];
			//		}

			//		return edges[0];

			//	}
			//	return edges[0];
			//}
			//std::cout << "--" << std::endl;
			//std::cout << 486 << "   " << edges[0].w << std::endl;
			//std::cout << 486 << "   " << edges[1].w << std::endl;

			const auto preNormalizeStructuralProfile = enableUniqueSemanticProbe
				? uniqueTable.profileStructuralAlternatives(e)
				: decltype(uniqueTable.profileStructuralAlternatives(e)){};
			if (enableRegressionDiagnostics && preNormalizeStructuralProfile.found) {
				regressionDiagnostics.uniquePreNormalizeStructuralAlternatives++;
			}
			if (enableRegressionDiagnostics && preNormalizeStructuralProfile.sameWeightDifferentMap) {
				regressionDiagnostics.uniquePreNormalizeSameWeightDifferentMap++;
			}
			if (enableRegressionDiagnostics && preNormalizeStructuralProfile.sameMapDifferentWeight) {
				regressionDiagnostics.uniquePreNormalizeSameMapDifferentWeight++;
			}
			if (enableRegressionDiagnostics && preNormalizeStructuralProfile.differentMapAndWeight) {
				regressionDiagnostics.uniquePreNormalizeDifferentMapAndWeight++;
			}

			e = normalize(e, cached);

			if (traceThisMakeNode) {
				std::cerr << "makeDDNode_trace_after_normalize\tstep\t" << contStageTraceStep
					<< "\tvar\t" << static_cast<int>(var)
					<< "\tret_w_re\t" << CTEntry::val(e.w.r)
					<< "\tret_w_im\t" << CTEntry::val(e.w.i)
					<< "\tret_var\t" << static_cast<int>(e.p->v)
					<< "\tret_map_level\t" << (e.map ? static_cast<int>(e.map->level) : -999)
					<< "\tret_map_x\t" << (e.map ? static_cast<int>(e.map->x) : -1)
					<< "\tret_map_rot\t" << (e.map ? e.map->rotate : -999)
					<< "\tret_map_ep\t" << (e.map ? e.map->extra_phase : -999)
					<< "\n";
				for (std::size_t edgeIndex = 0; edgeIndex < e.p->e.size(); ++edgeIndex) {
					std::cerr << "makeDDNode_trace_norm_edge\tstep\t" << contStageTraceStep
						<< "\tvar\t" << static_cast<int>(var)
						<< "\tedge\t" << edgeIndex
						<< "\tw_re\t" << CTEntry::val(e.p->e[edgeIndex].w.r)
						<< "\tw_im\t" << CTEntry::val(e.p->e[edgeIndex].w.i)
						<< "\tchild_var\t" << static_cast<int>(e.p->e[edgeIndex].p->v)
						<< "\tmap_level\t" << (e.p->e[edgeIndex].map ? static_cast<int>(e.p->e[edgeIndex].map->level) : -999)
						<< "\tmap_x\t" << (e.p->e[edgeIndex].map ? static_cast<int>(e.p->e[edgeIndex].map->x) : -1)
						<< "\tmap_rot\t" << (e.p->e[edgeIndex].map ? e.p->e[edgeIndex].map->rotate : -999)
						<< "\tmap_ep\t" << (e.p->e[edgeIndex].map ? e.p->e[edgeIndex].map->extra_phase : -999)
						<< "\n";
				}
			}

			assert(e.p->v == var || e.isTerminal());

			// look it up in the unique tables
			const auto semanticAlternative = enableUniqueSemanticProbe && uniqueTable.hasSemanticAlternative(e);
			const auto semanticGlobalAlternativeBucket = enableUniqueSemanticProbe && !semanticAlternative
				? uniqueTable.findSemanticAlternativeBucket(e)
				: decltype(uniqueTable.findSemanticAlternativeBucket(e)){uniqueTable.getTables().front().size()};
			const auto structuralAlternativeBucket = enableUniqueSemanticProbe
				? uniqueTable.findStructuralAlternativeBucket(e)
				: decltype(uniqueTable.findStructuralAlternativeBucket(e)){uniqueTable.getTables().front().size()};
			const auto structuralAlternativeProfile = enableUniqueSemanticProbe
				? uniqueTable.profileStructuralAlternatives(e)
				: decltype(uniqueTable.profileStructuralAlternatives(e)){};
			if (enableRegressionDiagnostics && semanticAlternative) {
				regressionDiagnostics.uniqueSemanticAlternatives++;
			}
			if (enableRegressionDiagnostics && semanticGlobalAlternativeBucket != uniqueTable.getTables().front().size()) {
				regressionDiagnostics.uniqueSemanticGlobalAlternatives++;
			}
			if (enableRegressionDiagnostics && structuralAlternativeBucket != uniqueTable.getTables().front().size()) {
				regressionDiagnostics.uniqueStructuralAlternatives++;
			}
			if (enableRegressionDiagnostics && structuralAlternativeProfile.sameWeightDifferentMap) {
				regressionDiagnostics.uniqueStructuralSameWeightDifferentMap++;
			}
			if (enableRegressionDiagnostics && structuralAlternativeProfile.sameMapDifferentWeight) {
				regressionDiagnostics.uniqueStructuralSameMapDifferentWeight++;
			}
			if (enableRegressionDiagnostics && structuralAlternativeProfile.differentMapAndWeight) {
				regressionDiagnostics.uniqueStructuralDifferentMapAndWeight++;
			}
			if (enableRegressionDiagnostics && preNormalizeStructuralProfile.differentMapAndWeight) {
				if (structuralAlternativeProfile.sameWeightDifferentMap) {
					regressionDiagnostics.uniqueNormalizeTransitionDiffBothToMapOnly++;
				}
				if (structuralAlternativeProfile.sameMapDifferentWeight) {
					regressionDiagnostics.uniqueNormalizeTransitionDiffBothToWeightOnly++;
				}
				if (structuralAlternativeProfile.differentMapAndWeight) {
					regressionDiagnostics.uniqueNormalizeTransitionDiffBothToDiffBoth++;
				}
			}
			auto l = uniqueTable.lookup(e, false);
			if (enableRegressionDiagnostics && semanticAlternative && l.p == e.p) {
				regressionDiagnostics.uniqueSemanticAlternativeMisses++;
			}
			if (enableRegressionDiagnostics && semanticGlobalAlternativeBucket != uniqueTable.getTables().front().size() && l.p == e.p) {
				regressionDiagnostics.uniqueSemanticGlobalAlternativeMisses++;
			}
			if (enableRegressionDiagnostics && structuralAlternativeBucket != uniqueTable.getTables().front().size() && l.p == e.p) {
				regressionDiagnostics.uniqueStructuralAlternativeMisses++;
			}

			assert(l.p->v == var || l.isTerminal());
			//std::cout << 486 << "   " << l.w << std::endl;
			//std::cout << "--" << std::endl;

			return l;
		}


		///
		/// Compute table definitions
		///
	public:
		void clearComputeTables() {

		}


	public:

		void traceMapChain(const char* label, const the_maps* map) const {
			if (!enableContStageTrace) {
				return;
			}
			std::cerr << label << "\tstep\t" << contStageTraceStep;
			int depth = 0;
			for (const auto* cur = map; cur != nullptr; cur = cur->father) {
				std::cerr << "\tnode" << depth
					<< "_addr\t" << cur
					<< "\tnode" << depth << "_level\t" << static_cast<int>(cur->level)
					<< "\tnode" << depth << "_x\t" << static_cast<int>(cur->x)
					<< "\tnode" << depth << "_rot\t" << cur->rotate
					<< "\tnode" << depth << "_ep\t" << cur->extra_phase;
				++depth;
				if (cur->level == -1) {
					break;
				}
			}
			std::cerr << "\n";
		}

		the_maps* append_new_map(the_maps* self, short level, bool x, int rotate) {

			rotate = (rotate % root_of_unit + root_of_unit) % root_of_unit;

			if (x == 0 && rotate==0) {
				return self;
			}

			std::string new_key = std::to_string(level) + "_" + std::to_string(x) + "_" + std::to_string(rotate);

			auto it = self->next.find(new_key);

			if (it != self->next.end()) {
				return it->second;
			}
			else {
				self->next[new_key] = new the_maps{ level, x, rotate,0,{}, self };
				//std::cout << 570 << " " << x << " " << rotate<< " " << rotate % root_of_unit << std::endl;
				self->next[new_key]->extra_phase = 0;
				return self->next[new_key];
			}
		}


		ComputeTable3 <the_maps*, the_maps*, the_maps*>  mapmulTable{};

		the_maps* mapmul(the_maps* self, the_maps* other) {

			const auto traceMapOps = enableContStageTrace;
			if (traceMapOps && ((self == the_maps::the_maps_header() && self->extra_phase != 0) || (other == the_maps::the_maps_header() && other->extra_phase != 0))) {
				std::cerr << "mapmul_header_phase_input\tstep\t" << contStageTraceStep
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tself_header\t" << (self == the_maps::the_maps_header())
					<< "\tself_ep\t" << self->extra_phase
					<< "\tother_header\t" << (other == the_maps::the_maps_header())
					<< "\tother_ep\t" << other->extra_phase
					<< "\n";
				traceMapChain("mapmul_input_self_chain", self);
				traceMapChain("mapmul_input_other_chain", other);
			}

			if (self->level == -1) {
				if (traceMapOps && other->extra_phase != 0) {
					std::cerr << "mapmul_base_reset_other\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tother_header\t" << (other == the_maps::the_maps_header())
						<< "\tother_level\t" << static_cast<int>(other->level)
						<< "\told_ep\t" << other->extra_phase
						<< "\n";
				}
				if (enableRegressionDiagnostics) {
					regressionDiagnostics.mapmulBaseResetOther++;
				}
				other->extra_phase = 0;
				return other;
			}

			if (other->level == -1) {
				if (traceMapOps && self->extra_phase != 0) {
					std::cerr << "mapmul_base_reset_self\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tself_header\t" << (self == the_maps::the_maps_header())
						<< "\tself_level\t" << static_cast<int>(self->level)
						<< "\told_ep\t" << self->extra_phase
						<< "\n";
				}
				if (enableRegressionDiagnostics) {
					regressionDiagnostics.mapmulBaseResetSelf++;
				}
				self->extra_phase = 0;
				return self;
			}

				auto r = mapmulTable.lookup(self, other);
				if (r != nullptr) {
					if (traceMapOps && (r == the_maps::the_maps_header() || r->extra_phase != 0)) {
						std::cerr << "mapmul_lookup_return\tstep\t" << contStageTraceStep
							<< "\tdepth\t" << contStageTraceDepth
							<< "\tresult_header\t" << (r == the_maps::the_maps_header())
							<< "\tresult_level\t" << static_cast<int>(r->level)
							<< "\tresult_ep\t" << r->extra_phase
							<< "\n";
						traceMapChain("mapmul_lookup_result_chain", r);
					}
				if (enableRegressionDiagnostics) {
					regressionDiagnostics.mapmulLookupHits++;
					if (r->extra_phase != 0) {
						regressionDiagnostics.mapmulLookupPhaseful++;
					}
				}
				if (disableMapdivLookupWriteback) {
					r->extra_phase = 0;
				}
				return r;
			}
			the_maps* res;
			if (self->level > other->level) {
				auto r = mapmul(self->father, other);
				res = append_new_map(r, self->level, self->x, self->rotate);
				res->extra_phase = r->extra_phase;
			}
			else if (self->level < other->level) {
				auto r = mapmul(self, other->father);
				res = append_new_map(r, other->level, other->x, other->rotate);
				res->extra_phase = r->extra_phase;
			}
			else {
				auto r = mapmul(self->father, other->father);
				//long int rotate = other->rotate + self->rotate * pow(-1, other->x);
				auto rotate = 0;
				if (other->x == 0) {
					rotate=other->rotate+self->rotate;
				}
				else {
					rotate=other->rotate-self->rotate;
				}

				res = append_new_map(r, self->level, (self->x + other->x) % 2, rotate%root_of_unit);
				res->extra_phase = r->extra_phase;
				if (other->x) {
					res->extra_phase = res->extra_phase+self->rotate;
				}
			}

				if (traceMapOps && (res == the_maps::the_maps_header() || res->extra_phase != 0)) {
					std::cerr << "mapmul_insert_result\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tresult_header\t" << (res == the_maps::the_maps_header())
						<< "\tresult_level\t" << static_cast<int>(res->level)
						<< "\tresult_ep\t" << res->extra_phase
						<< "\tstored_ep\t" << (res->extra_phase % root_of_unit)
						<< "\n";
					traceMapChain("mapmul_insert_result_chain", res);
				}
				mapmulTable.insert(self, other, res, res->extra_phase%root_of_unit);
			if (enableRegressionDiagnostics && res->extra_phase != 0) {
				regressionDiagnostics.mapmulResultPhaseful++;
			}

			return res;
		}

		ComputeTable3 <the_maps*, the_maps*, the_maps*>  mapdivTable{};

		the_maps* mapdiv(the_maps* self, the_maps* other) {

				const auto traceMapOps = enableContStageTrace;
				if (traceMapOps && ((self == the_maps::the_maps_header() && self->extra_phase != 0) || (other == the_maps::the_maps_header() && other->extra_phase != 0))) {
					std::cerr << "mapdiv_header_phase_input\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tself_header\t" << (self == the_maps::the_maps_header())
						<< "\tself_ep\t" << self->extra_phase
						<< "\tother_header\t" << (other == the_maps::the_maps_header())
						<< "\tother_ep\t" << other->extra_phase
						<< "\n";
					traceMapChain("mapdiv_input_self_chain", self);
					traceMapChain("mapdiv_input_other_chain", other);
				}

			if (other->level == -1) {
				if (traceMapOps && self->extra_phase != 0) {
					std::cerr << "mapdiv_base_reset_self\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tself_header\t" << (self == the_maps::the_maps_header())
						<< "\tself_level\t" << static_cast<int>(self->level)
						<< "\told_ep\t" << self->extra_phase
						<< "\n";
				}
				if (enableRegressionDiagnostics) {
					regressionDiagnostics.mapdivBaseResetSelf++;
				}
				self->extra_phase = 0;
				return self;
			}
			if (self == other) {
				auto the_maps_header = the_maps::the_maps_header();
				if (traceMapOps && the_maps_header->extra_phase != 0) {
					std::cerr << "mapdiv_same_return_header_reset\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\told_ep\t" << the_maps_header->extra_phase
						<< "\n";
				}
				if (enableRegressionDiagnostics) {
					regressionDiagnostics.mapdivHeaderReset++;
				}
				the_maps_header->extra_phase = 0;
				return the_maps_header;
			}
			
			if (const auto* entry = mapdivTable.findEntry(self, other); entry != nullptr) {
				if (traceMapOps && (entry->result == the_maps::the_maps_header() || entry->extra_phase != 0 || entry->result->extra_phase != entry->extra_phase)) {
					std::cerr << "mapdiv_lookup_return\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tresult_header\t" << (entry->result == the_maps::the_maps_header())
						<< "\tresult_level\t" << static_cast<int>(entry->result->level)
						<< "\tresult_ep_before\t" << entry->result->extra_phase
						<< "\tentry_ep\t" << entry->extra_phase
						<< "\n";
					traceMapChain("mapdiv_lookup_result_chain", entry->result);
				}
				if (enableRegressionDiagnostics) {
					regressionDiagnostics.mapdivLookupHits++;
					if (entry->extra_phase != 0) {
						regressionDiagnostics.mapdivLookupPhaseful++;
					}
					if (entry->extra_phase != entry->result->extra_phase) {
						regressionDiagnostics.mapdivLookupPhaseOverwrite++;
						if (entry->result != the_maps::the_maps_header()) {
							regressionDiagnostics.mapdivLookupPhaseOverwriteNonHeader++;
						}
					}
				}
				if (!disableMapdivLookupWriteback) {
					entry->result->extra_phase = entry->extra_phase;
				}
				return entry->result;
			}
			
			the_maps* res;
			if (self->level > other->level) {
				auto r = mapdiv(self->father, other);
				res = append_new_map(r, self->level, self->x, self->rotate);
				res->extra_phase = r->extra_phase;
			}
			else if (self->level < other->level) {
				auto r = mapdiv(self, other->father);

				if (other->x == 0) {
					if (mode == 2) {
						// auto temp = cn.getTemporary();
						// cn.div(temp, Complex::one, other->rotate);
						// res = append_new_map(r, other->level, other->x, cn.lookup(temp));
						res = append_new_map(r, other->level, other->x, (-other->rotate)%root_of_unit);
					}
					else {
						res = append_new_map(r, other->level, other->x, (-other->rotate)%root_of_unit);
					}
					
					res->extra_phase = r->extra_phase;
				}
				else {
					res = append_new_map(r, other->level, other->x, other->rotate);
					res->extra_phase= r->extra_phase;
					res->extra_phase=res->extra_phase-other->rotate;
				}
			}
			else {
				auto r = mapdiv(self->father, other->father);

				bool x = (self->x + other->x) % 2;

				auto rotate = 0;

				if (x==1) {
					rotate = self->rotate + other->rotate;
				}
				else {
					rotate = self->rotate-other->rotate;
				}
				res = append_new_map(r, self->level, x, rotate%root_of_unit);
				res->extra_phase = r->extra_phase;
				if (x == 1) {
					res->extra_phase = res->extra_phase - other->rotate;
				}
			}
			if (traceMapOps && (res == the_maps::the_maps_header() || res->extra_phase != 0)) {
				std::cerr << "mapdiv_insert_result\tstep\t" << contStageTraceStep
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tresult_header\t" << (res == the_maps::the_maps_header())
					<< "\tresult_level\t" << static_cast<int>(res->level)
					<< "\tresult_ep\t" << res->extra_phase
					<< "\tstored_ep\t" << (res->extra_phase % root_of_unit)
					<< "\n";
				traceMapChain("mapdiv_insert_result_chain", res);
			}
			mapdivTable.insert(self, other, res, res->extra_phase%root_of_unit);
			if (enableRegressionDiagnostics && res->extra_phase != 0) {
				regressionDiagnostics.mapdivResultPhaseful++;
			}
			return res;
		}


	public:
		//==========================================我写的========================================


		ComputeTable<mCachedEdge, mCachedEdge, mCachedEdge, Config::CT_VEC_ADD_NBUCKET>	addTable{};
		ComputeTable2<mEdge, mEdge, mCachedEdge, Config::CT_MAT_MAT_MULT_NBUCKET>  contTable{};

		key_2_new_key_node key_2_new_key_tree_header_element = { -1,-1,{},nullptr };
		key_2_new_key_node* key_2_new_key_tree_header = &key_2_new_key_tree_header_element;

		key_2_new_key_node* append_new_key(key_2_new_key_node* self, float new_key) {

			auto it = self->next.find(new_key);
			
			if (it != self->next.end()) {
				return self->next[new_key];
			}
			else {
				self->next[new_key] = new key_2_new_key_node{ short(self->level + 1), new_key, {}, self };


				return self->next[new_key];
			}
		}

		template <class Edge> Edge T_add(const Edge& x, const Edge& y) {

			return y;
		}


		//template <class LeftOperand, class RightOperand>
		TDD cont(TDD tdd1, TDD tdd2) {

			TDD res;

			std::vector<Index> var_out;
			std::vector<std::string> var_cont_temp;
			std::vector<std::string> var_cont;
			std::vector<std::string> var_out_key;

			int k;
			int k1;
			for (k = 0; k < tdd1.index_set.size(); ++k) {
				bool flag = true;

				for (k1 = 0; k1 < tdd2.index_set.size(); ++k1) {
					if (tdd2.index_set[k1].idx == tdd1.index_set[k].idx && tdd2.index_set[k1].key == tdd1.index_set[k].key) {
						var_cont_temp.push_back(tdd1.index_set[k].key);
						flag = false;
						break;
					}
				}
				if (flag) {
					var_out.push_back(tdd1.index_set[k]);
					var_out_key.push_back(tdd1.index_set[k].key);
				}
			}

			for (k = 0; k < tdd2.index_set.size(); ++k) {
				bool flag = true;
				for (k1 = 0; k1 < tdd1.index_set.size(); ++k1) {
					if (tdd1.index_set[k1].idx == tdd2.index_set[k].idx && tdd1.index_set[k1].key == tdd2.index_set[k].key) {
						flag = false;
						break;
					}
				}
				if (flag) {
					var_out.push_back(tdd2.index_set[k]);
					var_out_key.push_back(tdd2.index_set[k].key);
				}
			}
			for (k = 0; k < var_cont_temp.size(); ++k) {
				if (find(var_out_key.begin(), var_out_key.end(), var_cont_temp[k]) == var_out_key.end()) {
					if (find(var_cont.begin(), var_cont.end(), var_cont_temp[k]) == var_cont.end()) {
						var_cont.push_back(var_cont_temp[k]);
					}
				}
			}


			key_2_new_key_node* key_2_new_key1 = key_2_new_key_tree_header;
			key_2_new_key_node* key_2_new_key2 = key_2_new_key_tree_header;

			std::vector<std::string> new_key_2_index;
			k1 = 0;
			int k2 = 0;
			int new_key = 0;
			int m1 = tdd1.key_2_index.size();
			int m2 = tdd2.key_2_index.size();
			int repeat_time = 1;
			float last_cont_idx = -2;
			const auto contracted_key_step = 1.0F / static_cast<float>(3 * nqubits);

			while (k1 < m1 || k2 < m2) {

				if (k1 == m1) {
					for (k2; k2 < m2; ++k2)
					{
						key_2_new_key2 = append_new_key(key_2_new_key2, new_key);
						new_key_2_index.push_back(tdd2.key_2_index[k2]);
						new_key++;
					}
					break;
				}
				if (k2 == m2) {
					for (k1; k1 < m1; ++k1)
					{
						key_2_new_key1 = append_new_key(key_2_new_key1, new_key);
						new_key_2_index.push_back(tdd1.key_2_index[k1]);
						new_key++;
					}
					break;
				}

				if (varOrder[tdd1.key_2_index[k1]] < varOrder[tdd2.key_2_index[k2]]) {
					key_2_new_key1 = append_new_key(key_2_new_key1, new_key);
					new_key_2_index.push_back(tdd1.key_2_index[k1]);
					new_key++;
					k1++;
				}
				else if (varOrder[tdd1.key_2_index[k1]] > varOrder[tdd2.key_2_index[k2]]) {
					key_2_new_key2 = append_new_key(key_2_new_key2, new_key);
					new_key_2_index.push_back(tdd2.key_2_index[k2]);
					new_key++;
					k2++;
				}
				else if (find(var_out_key.begin(), var_out_key.end(), tdd1.key_2_index[k1]) == var_out_key.end()) {
					if (new_key - last_cont_idx <= 0.5) {
						last_cont_idx = last_cont_idx + contracted_key_step * repeat_time;
						repeat_time += 1;
						key_2_new_key1 = append_new_key(key_2_new_key1, last_cont_idx);
						key_2_new_key2 = append_new_key(key_2_new_key2, last_cont_idx);
						k1++;
						k2++;
					}
					else {
						key_2_new_key1 = append_new_key(key_2_new_key1, new_key - 0.5);
						key_2_new_key2 = append_new_key(key_2_new_key2, new_key - 0.5);
						last_cont_idx = new_key - 0.5;
						repeat_time = 1;
						k1++;
						k2++;
					}

				}
				else {
					key_2_new_key1 = append_new_key(key_2_new_key1, new_key);
					key_2_new_key2 = append_new_key(key_2_new_key2, new_key);
					new_key_2_index.push_back(tdd1.key_2_index[k1]);
					new_key++;
					k1++;
					k2++;
				}
			}

			res.index_set = var_out;
			res.key_2_index = new_key_2_index;

			if (to_test) {
				std::cout << "dd order: " << std::endl;
				for(auto it: this->varOrder){
					std::cout << it.first << " " << it.second ;
				}
				std::cout << std::endl;
				std::cout << "TDD1: ";
				for (const auto& element : tdd1.key_2_index) {
					std::cout << element << " ";
				}
				std::cout << std::endl;

				std::cout << "TDD2: ";
				for (const auto& element : tdd2.key_2_index) {
					std::cout << element << " ";
				}
				std::cout << std::endl;
			}



			[[maybe_unused]] const auto before = cn.cacheCount();
			// std::cout << "-----------" << std::endl;
			// std::cout << tdd1.e.w<<" "<<tdd2.e.w << std::endl;
			res.e = cont2(tdd1.e, tdd2.e, key_2_new_key1, key_2_new_key2, var_cont.size());

			if (to_test) {
				std::cout << "TDD: ";
				for (const auto& element : res.key_2_index) {
					std::cout << element << " ";
				}
				std::cout << std::endl;

			}

			var_out.clear();
			var_cont_temp.clear();
			var_out_key.clear();
			var_cont.clear();

			if (!res.e.w.exactlyZero() && !res.e.w.exactlyOne()) {
				//assert(res.e.w != Complex::zero);
				// cn.returnToCache(res.e.w);
				res.e.w = cn.lookup(res.e.w);
			}

			[[maybe_unused]] const auto after = cn.cacheCount();

			// assert(before == after);

			

			/*std::cout << tdd1.e.w.r->value << " " << tdd1.e.w.i->value << " " << tdd2.e.w.r->value << " " << tdd2.e.w.i->value << " " << res.e.w.r->value << " " << res.e.w.i->value << std::endl;*/
			//std::cout << "-----------" << std::endl;
			//the_maps::print_maps(tdd1.e.map);
			//the_maps::print_maps(tdd2.e.map);
			//the_maps::print_maps(res.e.map);
			//std::cout << tdd1.e.w << " " << tdd2.e.w <<" "<<res.e.w << std::endl;
			//std::cout << size(res.e)<<" "<< res.e.p->v << std::endl;
			//std::cout << "-----------" << std::endl;
			return res;
		}

		template <class Node>
		Edge<Node> renormalize(const Edge<Node>& e) {

			if (e.p->v == -1) {
				return e;
			}
			std::vector<Edge<Node>> edge(2);
			edge[0] = renormalize(Slicing(e, e.p->v, 0));
			edge[1] = renormalize(Slicing(e, e.p->v, 1));
			
			//std::cout << "-----aa----"<< ComplexTable<>::tolerance() << std::endl;

			auto e2 = makeDDNode(e.p->v, edge, true);

			return e2;
		}

	private:

		template <class Node>
		Edge<Node> Slicing(const Edge<Node>& e, int x, int c) {
			auto throwOnNonFiniteSlicing = [&](const char* stage, const Edge<Node>& edge) {
				if (!(enableContStageTrace && contStageTraceDepth > 0)) {
					return;
				}
				const auto real = CTEntry::val(edge.w.r);
				const auto imag = CTEntry::val(edge.w.i);
				if (std::isfinite(real) && std::isfinite(imag)) {
					return;
				}
				std::ostringstream message;
				message << "Slicing non-finite edge"
					<< " step=" << contStageTraceStep
					<< " stage=" << stage
					<< " slice_var=" << x
					<< " branch=" << c
					<< " edge_var=" << static_cast<int>(edge.p->v)
					<< " weight=(" << real << ", " << imag << ")"
					<< " map_extra_phase=" << (edge.map ? edge.map->extra_phase : -999999)
					<< " input_var=" << static_cast<int>(e.p->v)
					<< " input_map_extra_phase=" << (e.map ? e.map->extra_phase : -999999);
				throw std::runtime_error(message.str());
			};
			// used for add
			assert(e.w != Complex::zero);
			throwOnNonFiniteSlicing("entry", e);
			if (e.p->v == -1) {
				return e;
			}
			if (e.p->v < x) {
				return e;
			}
			if (e.p->v == x) {
				if (e.p->v != e.map->level) {
					auto temp = e.p->e[c];
					throwOnNonFiniteSlicing("child_raw_direct", temp);
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = mapmul(e.map, temp.map);

						assert(temp.w != Complex::zero);
						// cn.mul(temp.w, temp.w, temp.map->extra_phase);
						cn.mul(temp.w, temp.w, cn.getTemporary(cos(temp.map->extra_phase*rotate_angle),sin(temp.map->extra_phase*rotate_angle)));
						throwOnNonFiniteSlicing("child_after_direct", temp);
						// cn.returnToCache(temp.map->extra_phase);
					}
					//std::cout << "Slicing " << temp.w << std::endl;
					return temp;
				}
				else if (e.map->x == 0) {
					auto temp = e.p->e[c];
					throwOnNonFiniteSlicing("child_raw_x0", temp);
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = mapmul(e.map->father, temp.map);

						assert(temp.w != Complex::zero);
						// cn.mul(temp.w, temp.w, temp.map->extra_phase);
						cn.mul(temp.w, temp.w, cn.getTemporary(cos(temp.map->extra_phase*rotate_angle),sin(temp.map->extra_phase*rotate_angle)));
						// cn.returnToCache(temp.map->extra_phase);

						if (c == 1) {
								assert(temp.w != Complex::zero);
								// cn.mul(temp.w, temp.w, e.map->rotate);
								cn.mul(temp.w, temp.w, cn.getTemporary(cos(e.map->rotate*rotate_angle),sin(e.map->rotate*rotate_angle)));
						}
						throwOnNonFiniteSlicing("child_after_x0", temp);

					}
					//std::cout << "Slicing " << temp.w << std::endl;
					return temp;
				}
				else {
					auto temp = e.p->e[1 - c];
					throwOnNonFiniteSlicing("child_raw_x1", temp);
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = mapmul(e.map->father, temp.map);

						assert(temp.w != Complex::zero);
						// cn.mul(temp.w, temp.w, temp.map->extra_phase);
						cn.mul(temp.w, temp.w, cn.getTemporary(cos(temp.map->extra_phase*rotate_angle),sin(temp.map->extra_phase*rotate_angle)));
						if (c == 0) {

							assert(temp.w != Complex::zero);
							// cn.mul(temp.w, temp.w, e.map->rotate);
                            cn.mul(temp.w, temp.w, cn.getTemporary(cos(e.map->rotate*rotate_angle),sin(e.map->rotate*rotate_angle)));
						}
						throwOnNonFiniteSlicing("child_after_x1", temp);
					}
					//std::cout << "Slicing " << temp.w << std::endl;
					return temp;
				}

			}
			else {
				std::cout << "Slicing not support yet" << std::endl;
				return e;
			}

		}

		template<class Node>
		Edge<Node>& copyEdge(const Edge<Node>& edge) {
			auto temp = edge;
			// if (cn.inCache(edge.w)) {
			// 	std::cout << "959: complex number not in cache" << std::endl;
			// 	temp.w = cn.getCached();
			// 	temp.w.r->value = edge.w.r->value;
			// 	temp.w.i->value = edge.w.i->value;
			// 	// temp = deepCopyEdge(edge);
			// }
			// else{
			// 	temp = edge;
			// }
			// temp.w = cn.getCached(edge.w.r,edge.w.i);
			//std::cout << temp.w << " value in: " << &(temp.w) << std::endl;
			//std::cout << "temp in table? " << cn.inTable(temp.w) << std::endl;
			//std::cout << "temp in cache? " << cn.inCache(temp.w) << std::endl;
			return temp;
		}
		void returnToCache(Complex& c){
			if(!c.exactlyZero() && !c.exactlyOne()){
				//std::cout << "975: return to cache: " << c.r->value << ","<< c.i->value << std::endl;
				//std::cout << Complex::zero <<" zero?: "<< (c==Complex::zero) << " " << Complex::one << " one?: "<< (c == Complex::one)  << std::endl;
				cn.returnToCache(c);
			}
		}


		template <class Node>
		Edge<Node> Slicing2(Edge<Node>& e, int x, int c) {

			assert(e.w != Complex::zero);
			const auto traceSlicing2 = std::getenv("LIMTDD_SLICING2_TRACE") != nullptr && enableContStageTrace;
			const auto traceSlicing2ForDepth = traceSlicing2 &&
				(std::getenv("LIMTDD_FOCUSED_CONT_DEPTH") == nullptr ||
				 contStageTraceDepth == static_cast<std::size_t>(std::stoull(std::getenv("LIMTDD_FOCUSED_CONT_DEPTH"))));
			if (traceSlicing2ForDepth) {
				std::cerr << "slicing2_entry\tx\t" << x << "\tc\t" << c
					<< "\tdepth\t" << contStageTraceDepth
					<< "\te_v\t" << static_cast<int>(e.p->v)
					<< "\te_w_re\t" << CTEntry::val(e.w.r)
					<< "\te_w_im\t" << CTEntry::val(e.w.i)
					<< "\te_map_level\t" << (e.map ? static_cast<int>(e.map->level) : -999)
					<< "\te_map_x\t" << (e.map ? static_cast<int>(e.map->x) : -1)
					<< "\te_map_rot\t" << (e.map ? e.map->rotate : -999)
					<< "\te_map_ep\t" << (e.map ? e.map->extra_phase : -999)
					<< "\n";
			}
		// used for contract
			if (e.p->v == -1) {
				if (traceSlicing2ForDepth) {
					std::cerr << "slicing2_branch\tterminal\tx\t" << x << "\tc\t" << c << "\tdepth\t" << contStageTraceDepth << "\n";
				}
				return e;
			}
			if (e.p->v < x) {
				if (traceSlicing2ForDepth) {
					std::cerr << "slicing2_branch\tpass_through\tx\t" << x << "\tc\t" << c
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tout_v\t" << static_cast<int>(e.p->v)
						<< "\tout_w_re\t" << CTEntry::val(e.w.r)
						<< "\tout_map_level\t" << (e.map ? static_cast<int>(e.map->level) : -999)
						<< "\n";
				}
				return e;
			}
			if (e.p->v == x) {
				if (e.p->v != e.map->level) {
					Edge<Node> temp = e.p->e[c];
					if (temp.w != Complex::zero) {
						temp.map = mapmul(e.map, temp.map);
						temp.w = cn.mulCached(temp.w, cn.getTemporary(cos(temp.map->extra_phase * rotate_angle), sin(temp.map->extra_phase * rotate_angle)));
					}
					if (traceSlicing2ForDepth) {
						std::cerr << "slicing2_branch\tv_eq_x_level_mismatch\tx\t" << x << "\tc\t" << c
							<< "\tdepth\t" << contStageTraceDepth
							<< "\tout_w_re\t" << CTEntry::val(temp.w.r)
							<< "\tout_w_im\t" << CTEntry::val(temp.w.i)
							<< "\tout_v\t" << static_cast<int>(temp.p ? temp.p->v : -999)
							<< "\tout_map_level\t" << (temp.map ? static_cast<int>(temp.map->level) : -999)
							<< "\tout_map_x\t" << (temp.map ? static_cast<int>(temp.map->x) : -1)
							<< "\tout_map_rot\t" << (temp.map ? temp.map->rotate : -999)
							<< "\n";
					}
					return temp;
				}
				else if (e.map->x == 0) {
					Edge<Node> temp = e.p->e[c];
					temp.w = cn.lookup(e.p->e[c].w);
					//std::cout << "1011 temp w: " << temp->w << " " << temp->w.i << " " << temp->w.r << " " << temp->p << std::endl;
					//std::cout << "979 w: " << e.p->e[c].w.i << " " << e.p->e[c].w.r << " " << e.p->e[c].p << std::endl;
					// std::cout << "1012 ref count:" << temp->w.i->refCount << " " << temp->w.r->refCount << std::endl;
					// std::cout << "979: " <<  & (e.p->e[c]) << " "<<& (temp) << std::endl;
					if (temp.w != Complex::zero) {
						temp.map = mapmul(e.map->father, temp.map);
						// if(temp->w == Complex::one) {
						// 	temp->w = cn.getCached(1., 0.)
						// }
						// temp->w = cn.mulCached(temp->w, temp->map->extra_phase);
						//std::cout << "Scling2 2 " << temp->w << std::endl;
						temp.w = cn.mulCached(temp.w, cn.getTemporary(cos(temp.map->extra_phase * rotate_angle), sin(temp.map->extra_phase * rotate_angle)));
						//std::cout << "1018 temp w: " << temp->w << " " << temp->w.i << " " << temp->w.r << std::endl;
						// cn.returnToCache(temp->map->extra_phase);
						//std::cout << "Scling2 2 " << temp->w << std::endl;
						if (c == 1) {
							assert(temp.w != Complex::zero);
							// cn.mul(temp->w, temp->w, e.map->rotate);
							cn.mul(temp.w, temp.w, cn.getTemporary(cos(e.map->rotate * rotate_angle), sin(e.map->rotate * rotate_angle)));
							//std::cout << "Scling2 2 " << temp->w << std::endl;
							//std::cout<< e.map->rotate<<" "<< e.map->rotate * rotate_angle << " " << cos(e.map->rotate * rotate_angle) << " " << sin(e.map->rotate * rotate_angle) <<std::endl;
						}
						// std::cout << "1021 ref count:" << temp->w.i->refCount << " " << temp->w.r->refCount << std::endl;
						// std::cout << "1020 temp w: " << temp->w << " in:" << temp->w.i << " " << temp->w.r << std::endl;
						// return {temp->p,  temp_w, temp_map};
					}
					//std::cout << "Scling2 2 " << temp->w << std::endl;
					if (traceSlicing2ForDepth) {
						std::cerr << "slicing2_branch\tv_eq_x_level_match_x0\tx\t" << x << "\tc\t" << c
							<< "\tdepth\t" << contStageTraceDepth
							<< "\te_map_rot\t" << (e.map ? e.map->rotate : -999)
							<< "\tout_w_re\t" << CTEntry::val(temp.w.r)
							<< "\tout_w_im\t" << CTEntry::val(temp.w.i)
							<< "\tout_v\t" << static_cast<int>(temp.p ? temp.p->v : -999)
							<< "\tout_map_level\t" << (temp.map ? static_cast<int>(temp.map->level) : -999)
							<< "\tout_map_x\t" << (temp.map ? static_cast<int>(temp.map->x) : -1)
							<< "\tout_map_rot\t" << (temp.map ? temp.map->rotate : -999)
							<< "\n";
					}
					return temp;
				}
				else {
					Edge<Node> temp = e.p->e[1-c];
					// std::cout << "1026: " << & (temp->w) << std::endl;
					// std::cout << "1029 ref count:" << temp->w.i->refCount << " " << temp->w.r->refCount << std::endl;
					if (temp.w != Complex::zero) {
						temp.map = mapmul(e.map->father, temp.map);
						// temp->w = cn.mulCached(temp->w, temp->map->extra_phase);
						temp.w = cn.mulCached(temp.w, cn.getTemporary(cos(temp.map->extra_phase * rotate_angle), sin(temp.map->extra_phase * rotate_angle)));
						// cn.returnToCache(temp->map->extra_phase);
						if (c == 0) {
							assert(temp.w != Complex::zero);
							// cn.mul(temp->w, temp->w, e.map->rotate);
							cn.mul(temp.w, temp.w, cn.getTemporary(cos(e.map->rotate * rotate_angle), sin(e.map->rotate * rotate_angle)));

						}
						// std::cout << "1038 ref count:" << temp->w.i->refCount << " " << temp->w.r->refCount << std::endl;
						// std::cout << "1037 temp w: " << temp->w << " in:" << &(temp->w.i) << std::endl;
						// return {temp.p,  temp_w, temp_map};
					}
					//std::cout << "Scling2 3 " << temp->w << std::endl;
					if (traceSlicing2ForDepth) {
						std::cerr << "slicing2_branch\tv_eq_x_level_match_x1\tx\t" << x << "\tc\t" << c
							<< "\tdepth\t" << contStageTraceDepth
							<< "\te_map_rot\t" << (e.map ? e.map->rotate : -999)
							<< "\tout_w_re\t" << CTEntry::val(temp.w.r)
							<< "\tout_w_im\t" << CTEntry::val(temp.w.i)
							<< "\tout_v\t" << static_cast<int>(temp.p ? temp.p->v : -999)
							<< "\tout_map_level\t" << (temp.map ? static_cast<int>(temp.map->level) : -999)
							<< "\tout_map_x\t" << (temp.map ? static_cast<int>(temp.map->x) : -1)
							<< "\tout_map_rot\t" << (temp.map ? temp.map->rotate : -999)
							<< "\n";
					}
					return temp;
				}

			}
			else {
				std::cout << "Slicing2 not support yet" << std::endl;
				return e;
			}

		}



		template <class Node>
		Edge<Node> T_add2(const Edge<Node>& x, const Edge<Node>& y) {
			const auto traceTadd = enableContStageTrace && contStageTraceDepth > 0;
			const auto traceBefore = traceTadd ? regressionDiagnostics : RegressionDiagnostics{};
			auto throwOnNonFiniteEdge = [&](const char* stage, const Edge<Node>& edge) {
				if (!traceTadd) {
					return;
				}
				const auto real = CTEntry::val(edge.w.r);
				const auto imag = CTEntry::val(edge.w.i);
				if (std::isfinite(real) && std::isfinite(imag)) {
					return;
				}
				std::ostringstream message;
				message << "T_add2 non-finite edge"
					<< " step=" << contStageTraceStep
					<< " stage=" << stage
					<< " var=" << static_cast<int>(edge.p->v)
					<< " weight=(" << real << ", " << imag << ")"
					<< " map_extra_phase=" << (edge.map ? edge.map->extra_phase : -999999);
				throw std::runtime_error(message.str());
			};
			auto throwOnNonFiniteNamedEdge = [&](const char* stage,
											const Edge<Node>& edge,
											const char* source,
											const std::size_t childIndex,
											const Edge<Node>& lhs,
											const Edge<Node>& rhs) {
				if (!traceTadd) {
					return;
				}
				const auto real = CTEntry::val(edge.w.r);
				const auto imag = CTEntry::val(edge.w.i);
				if (std::isfinite(real) && std::isfinite(imag)) {
					return;
				}
				std::ostringstream message;
				message << "T_add2 non-finite edge"
					<< " step=" << contStageTraceStep
					<< " stage=" << stage
					<< " source=" << source
					<< " child_index=" << childIndex
					<< " edge_var=" << static_cast<int>(edge.p->v)
					<< " weight=(" << real << ", " << imag << ")"
					<< " map_extra_phase=" << (edge.map ? edge.map->extra_phase : -999999)
					<< " lhs_var=" << static_cast<int>(lhs.p->v)
					<< " rhs_var=" << static_cast<int>(rhs.p->v);
				throw std::runtime_error(message.str());
			};
			auto throwOnBadDivInputs = [&](const Complex& numerator,
									   const Complex& denominator,
									   const Edge<Node>& lhs,
									   const Edge<Node>& rhs) {
				if (!traceTadd) {
					return;
				}
				const auto nr = CTEntry::val(numerator.r);
				const auto ni = CTEntry::val(numerator.i);
				const auto dr = CTEntry::val(denominator.r);
				const auto di = CTEntry::val(denominator.i);
				const auto denomMag2 = dr * dr + di * di;
				const auto numerFinite = std::isfinite(nr) && std::isfinite(ni);
				const auto denomFinite = std::isfinite(dr) && std::isfinite(di) && std::isfinite(denomMag2);
				const auto denomApproxZero = std::abs(denomMag2) < ComplexTable<>::tolerance();
				if (numerFinite && denomFinite) {
					return;
				}
				std::ostringstream message;
				message << "T_add2 bad div inputs"
					<< " step=" << contStageTraceStep
					<< " numerator=(" << nr << ", " << ni << ")"
					<< " denominator=(" << dr << ", " << di << ")"
					<< " denominator_mag2=" << denomMag2
					<< " denominator_approx_zero=" << denomApproxZero
					<< " lhs_var=" << static_cast<int>(lhs.p->v)
					<< " rhs_var=" << static_cast<int>(rhs.p->v)
					<< " lhs_map_extra_phase=" << (lhs.map ? lhs.map->extra_phase : -999999)
					<< " rhs_map_extra_phase=" << (rhs.map ? rhs.map->extra_phase : -999999);
				throw std::runtime_error(message.str());
			};
			if (traceTadd) {
				contStageTaddDepth++;
			}
			struct TAddTraceGuard {
				Package* pkg;
				bool active;
				RegressionDiagnostics before;
				~TAddTraceGuard() {
					if (!active) {
						return;
					}
					if (pkg->contStageTaddDepth == 1) {
						Package::accumulateRegressionDelta(pkg->contStageTaddTotals, before, pkg->regressionDiagnostics);
					}
					pkg->contStageTaddDepth--;
				}
			} taddTraceGuard{this, traceTadd, traceBefore};

			if (x.p > y.p) {
				return T_add2(y, x);
			}

			if (x.w.approximatelyZero()) {
				if (y.w.approximatelyZero()) {
					return Edge<Node>::zero;
				}
				auto r = y;
				r.w = cn.getCached(CTEntry::val(y.w.r), CTEntry::val(y.w.i));
				if (r.w.approximatelyZero()) {
					return Edge<Node>::zero;
				}
				return r;
			}
			if (y.w.approximatelyZero()) {
				auto r = x;
				r.w = cn.getCached(CTEntry::val(x.w.r), CTEntry::val(x.w.i));
				if (r.w.approximatelyZero()) {
					return Edge<Node>::zero;
				}
				return r;
			}
			if (x.p == y.p && x.map==y.map) {
				//std::cout << "Case 0" << std::endl;
				auto r = y;
				r.w = cn.addCached(x.w, y.w);
				if (r.w.approximatelyZero()) {
					//assert(r.w != Complex::zero);
					// cn.returnToCache(r.w);
					return Edge<Node>::zero;
				}
				r.map = x.map;
				
				return r;
			}
			const auto samePointerMapMismatch = x.p == y.p && x.map != y.map;
			if (enableRegressionDiagnostics && samePointerMapMismatch) {
				regressionDiagnostics.taddSamePointerMapMismatch++;
			}

			auto xCopy = x;
			auto yCopy = y;


			xCopy.w = Complex::one;
			xCopy.map = the_maps::the_maps_header();
			throwOnBadDivInputs(y.w, x.w, x, y);
			auto divResult = cn.getCached();
			if (traceTadd) {
				const auto aliasesYR = divResult.r == y.w.r;
				const auto aliasesYI = divResult.i == y.w.i;
				const auto aliasesXR = divResult.r == x.w.r;
				const auto aliasesXI = divResult.i == x.w.i;
				if (aliasesYR || aliasesYI || aliasesXR || aliasesXI) {
					std::ostringstream message;
					message << "T_add2 div cache alias"
						<< " step=" << contStageTraceStep
						<< " result_r=" << divResult.r
						<< " result_i=" << divResult.i
						<< " y_r=" << y.w.r
						<< " y_i=" << y.w.i
						<< " x_r=" << x.w.r
						<< " x_i=" << x.w.i;
					throw std::runtime_error(message.str());
				}
			}
			const auto manualNr = CTEntry::val(y.w.r);
			const auto manualNi = CTEntry::val(y.w.i);
			const auto manualDr = CTEntry::val(x.w.r);
			const auto manualDi = CTEntry::val(x.w.i);
			const auto manualCmag = manualDr * manualDr + manualDi * manualDi;
			const auto manualReal = (manualNr * manualDr + manualNi * manualDi) / manualCmag;
			const auto manualImag = (manualNi * manualDr - manualNr * manualDi) / manualCmag;
			ComplexNumbers::div(divResult, y.w, x.w);
			if (traceTadd) {
				const auto actualReal = CTEntry::val(divResult.r);
				const auto actualImag = CTEntry::val(divResult.i);
				if (!std::isfinite(actualReal) || !std::isfinite(actualImag)) {
					std::ostringstream message;
					message << "T_add2 div result non-finite"
						<< " step=" << contStageTraceStep
						<< " numerator=(" << manualNr << ", " << manualNi << ")"
						<< " denominator=(" << manualDr << ", " << manualDi << ")"
						<< " numerator_exact_zero=" << y.w.exactlyZero()
						<< " numerator_approx_zero=" << y.w.approximatelyZero()
						<< " denominator_exact_zero=" << x.w.exactlyZero()
						<< " denominator_approx_zero=" << x.w.approximatelyZero()
						<< " manual_cmag=" << manualCmag
						<< " manual_real=" << manualReal
						<< " manual_imag=" << manualImag
						<< " actual_real=" << actualReal
						<< " actual_imag=" << actualImag
						<< " result_r_ptr=" << divResult.r
						<< " result_i_ptr=" << divResult.i
						<< " numerator_r_ptr=" << y.w.r
						<< " numerator_i_ptr=" << y.w.i
						<< " denominator_r_ptr=" << x.w.r
						<< " denominator_i_ptr=" << x.w.i;
					throw std::runtime_error(message.str());
				}
			}
			yCopy.w = divResult;
			throwOnNonFiniteNamedEdge("ycopy_after_div_only", yCopy, "yCopy", 0, x, y);
			yCopy.map = mapdiv(y.map, x.map);
			throwOnNonFiniteNamedEdge("ycopy_after_div", yCopy, "yCopy", 0, x, y);
			if (enableRegressionDiagnostics && samePointerMapMismatch) {
				if (yCopy.map == the_maps::the_maps_header()) {
					regressionDiagnostics.taddMismatchResidualHeader++;
				} else {
					regressionDiagnostics.taddMismatchResidualNonHeader++;
				}
				if (yCopy.map->extra_phase != 0) {
					regressionDiagnostics.taddMismatchResidualPhaseful++;
				}
			}
			if (yCopy.w != Complex::zero) {
				// cn.mul(yCopy.w, yCopy.w, yCopy.map->extra_phase);
				cn.mul(yCopy.w, yCopy.w, cn.getTemporary(cos(yCopy.map->extra_phase*rotate_angle),sin(yCopy.map->extra_phase*rotate_angle)));
				
			}
			throwOnNonFiniteNamedEdge("ycopy_after_phase", yCopy, "yCopy", 0, x, y);
			// cn.returnToCache(yCopy.map->extra_phase);


			const auto leftLookup = CachedEdge<Node>{ xCopy.p, xCopy.w, xCopy.map };
			const auto rightLookup = CachedEdge<Node>{ yCopy.p, yCopy.w, yCopy.map };
			auto r = addTable.lookup(leftLookup, rightLookup);
			if (enableRegressionDiagnostics && samePointerMapMismatch) {
				if (r.p != nullptr) {
					regressionDiagnostics.taddMismatchAddHits++;
				} else {
					regressionDiagnostics.taddMismatchAddMisses++;
					if (const auto* entry = addTable.findEntry(leftLookup, rightLookup); entry == nullptr) {
						regressionDiagnostics.taddMismatchAddMissEmpty++;
					} else {
						const auto sameMap = entry->leftOperand.map == leftLookup.map && entry->rightOperand.map == rightLookup.map;
						const auto sameWeight = entry->leftOperand.w.approximatelyEquals(leftLookup.w) && entry->rightOperand.w.approximatelyEquals(rightLookup.w);
						if (!sameMap && sameWeight) {
							regressionDiagnostics.taddMismatchAddMissMapOnly++;
						} else if (sameMap && !sameWeight) {
							regressionDiagnostics.taddMismatchAddMissWeightOnly++;
						} else {
							regressionDiagnostics.taddMismatchAddMissMapAndWeight++;
						}
					}
				}
			}

			if (r.p != nullptr) {
				//std::cout << "Case 1" << std::endl;
				//assert(yCopy.w != Complex::zero);
				// cn.returnToCache(yCopy.w);
				if (r.w.approximatelyZero()) {
					return Edge<Node>::zero;
				}
				auto c = cn.getCached(r.w);

				if (c != Complex::zero) {
					cn.mul(c, c, x.w);
				}

				auto temp_map = mapmul(x.map, r.map);
				if (c != Complex::zero) {
					// cn.mul(c, c, temp_map->extra_phase);
					cn.mul(c, c, cn.getTemporary(cos(temp_map->extra_phase*rotate_angle),sin(temp_map->extra_phase*rotate_angle)));
				}
				// cn.returnToCache(temp_map->extra_phase);
				auto result = Edge<Node>{ r.p, c,temp_map };
				throwOnNonFiniteEdge("lookup_return", result);
				return result;
			}

			const Qubit w = (x.isTerminal() || (!y.isTerminal() && y.p->v > x.p->v))
				? y.p->v
				: x.p->v;

			int n = (x.p->v != w)? y.p->e.size() : x.p->e.size();

			std::vector<Edge<Node>> edge(n);
			for (std::size_t i = 0U; i < n; i++) {
				Edge<Node> e1{};
				if (!x.isTerminal() && x.p->v == w) {
					//e1 = x.p->e[i];
					//if (e1.w != Complex::zero) {
					//	e1.w = cn.mulCached(e1.w, xCopy.w);
					//}
					e1 = Slicing(xCopy, xCopy.p->v, i);
				}
				else {
					e1 = xCopy;
					if (y.p->e[i].p == nullptr) {
						e1 = { nullptr, Complex::zero };
					}
				}
				Edge<Node> e2{};
				if (!y.isTerminal() && y.p->v == w) {
					//e2 = y.p->e[i];
					//if (e2.w != Complex::zero) {
					//	e2.w = cn.mulCached(e2.w, yCopy.w);
					//}
					e2 = Slicing(yCopy, yCopy.p->v, i);
					throwOnNonFiniteNamedEdge("rhs_after_slicing", e2, "Slicing(yCopy)", i, e1, yCopy);
				}
				else {
					e2 = yCopy;
					if (x.p->e[i].p == nullptr) {
						e2 = { nullptr, Complex::zero };
					}
					throwOnNonFiniteNamedEdge("rhs_passthrough", e2, "yCopy", i, e1, yCopy);
				}
				throwOnNonFiniteEdge("child_input_lhs", e1);
				throwOnNonFiniteEdge("child_input_rhs", e2);
				edge[i] = T_add2(e1, e2);
				throwOnNonFiniteEdge("child_sum", edge[i]);


				if (!x.isTerminal() && x.p->v == w && e1.w != Complex::zero) {
					//assert(e1.w != Complex::zero);
					// cn.returnToCache(e1.w);
				}

				if (!y.isTerminal() && y.p->v == w && e2.w != Complex::zero) {
					//assert(e2.w != Complex::zero);
					// cn.returnToCache(e2.w);
				}
			}
			if(to_test){
					std::cout << "T_add2 function: " << 1083 << std::endl; 
					std::cout << "var: " << w << std::endl;
					std::cout << "edge.node.key: " << std::endl;
					for(auto e:edge){
						std::cout << e.p->v << " " ;
						if(e.p->v == w) throw std::runtime_error("bug here");
					}
					std::cout << std::endl;

				}
			auto e = makeDDNode(w, edge, true);

			addTable.insert({ xCopy.p,xCopy.w,xCopy.map }, { yCopy.p,yCopy.w,yCopy.map }, { e.p, e.w,e.map });
			//if (x.w != Complex::one) {
	
			//	assert(e.w != Complex::zero);
			//	//assert(yCopy.w != Complex::zero);
			//	
			//}
			if (e.w != Complex::zero) {
				assert(e.w != Complex::zero);
				cn.mul(e.w, e.w, x.w);
				e.map = mapmul(x.map, e.map);

				assert(e.w != Complex::zero);
				// cn.mul(e.w, e.w, e.map->extra_phase);
				cn.mul(e.w, e.w, cn.getTemporary(cos(e.map->extra_phase*rotate_angle),sin(e.map->extra_phase*rotate_angle)));
				// cn.returnToCache(e.map->extra_phase);
			}
			throwOnNonFiniteEdge("final_return", e);

			// cn.returnToCache(yCopy.w);
			//std::cout << "Case 2" << std::endl;
			return e;
		}


		comm_maps* find_remain_map(the_maps* map1, the_maps* map2, key_2_new_key_node* key_2_new_key1, key_2_new_key_node* key_2_new_key2) {

			//the_maps* res[3];
			//std::cout << 868 << "   " << map1->level << " " << map2->level << std::endl;

			//int to_tset2 = 2;
			//if (to_tset2 == 1) {
			//	comm_maps* res = new comm_maps{ the_maps::the_maps_header(),map1,map2 };
			//	res->remain_map->extra_phase = cn.getCached(1, 0);
			//	return res;
			//}


			key_2_new_key_node* temp_key_2_new_key1 = key_2_new_key1;
			while (temp_key_2_new_key1->level > map1->level) {
				temp_key_2_new_key1 = temp_key_2_new_key1->father;
			}

			key_2_new_key_node* temp_key_2_new_key2 = key_2_new_key2;
			while (temp_key_2_new_key2->level > map2->level) {
				temp_key_2_new_key2 = temp_key_2_new_key2->father;
			}

			float newk1 = temp_key_2_new_key1->new_key;
			float newk2 = temp_key_2_new_key2->new_key;

			if (newk1 > newk2 && !ifContract(newk1)) {
				auto res = find_remain_map(map1->father, map2, temp_key_2_new_key1, temp_key_2_new_key2);
				auto temp_pahse = res->remain_map->extra_phase;
				res->remain_map = append_new_map(res->remain_map, newk1, map1->x, map1->rotate);
				if (enableRegressionDiagnostics && temp_pahse != 0) {
					regressionDiagnostics.findRemainPhaseCarries++;
				}
				res->remain_map->extra_phase = temp_pahse;
				return res;
			}
			if (newk1 < newk2 && !ifContract(newk2)) {
				auto res = find_remain_map(map1, map2->father, temp_key_2_new_key1, temp_key_2_new_key2);
				auto temp_pahse = res->remain_map->extra_phase;
				res->remain_map = append_new_map(res->remain_map, newk2, map2->x, map2->rotate);
				if (enableRegressionDiagnostics && temp_pahse != 0) {
					regressionDiagnostics.findRemainPhaseCarries++;
				}
				res->remain_map->extra_phase = temp_pahse;
				return res;
			}
			if (map1->level == -1 && map2->level == -1) {
				if (enableContStageTrace && (map1->extra_phase != 0 || map2->extra_phase != 0)) {
					std::cerr << "find_remain_header_phase_input\tstep\t" << contStageTraceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tmap1_header\t" << (map1 == the_maps::the_maps_header())
						<< "\tmap1_ep\t" << map1->extra_phase
						<< "\tmap2_header\t" << (map2 == the_maps::the_maps_header())
						<< "\tmap2_ep\t" << map2->extra_phase
						<< "\n";
					traceMapChain("find_remain_input_map1_chain", map1);
					traceMapChain("find_remain_input_map2_chain", map2);
				}
				comm_maps* res=new comm_maps{ the_maps::the_maps_header(),the_maps::the_maps_header(),the_maps::the_maps_header() };
				res->remain_map->extra_phase = 0;
				return res;
			}
			if (newk1 > newk2) {
				auto res = find_remain_map(map1->father, map2, temp_key_2_new_key1, temp_key_2_new_key2);
				res->cont_map1 = append_new_map(res->cont_map1, map1->level, map1->x, map1->rotate);
				return res;
			}
			if (newk1 < newk2) {
				auto res = find_remain_map(map1, map2->father, temp_key_2_new_key1, temp_key_2_new_key2);
				res->cont_map2 = append_new_map(res->cont_map2, map2->level, map2->x, map2->rotate);
				return res;
			}
			auto res = find_remain_map(map1->father, map2->father, temp_key_2_new_key1, temp_key_2_new_key2);

			auto x = (map1->x + map2->x) % 2;
			if (x == 1) {
				// assert(res->remain_map->extra_phase != Complex::zero);
				if (enableRegressionDiagnostics) {
					regressionDiagnostics.findRemainPhaseCarries++;
				}
res->remain_map->extra_phase =  res->remain_map->extra_phase+ map2->rotate;
			}

			auto rotate = 0;
			if (x == 0) {
				// assert(rotate != Complex::zero);
				rotate = map1->rotate+map2->rotate;
			}
			else {
				rotate = map1->rotate- map2->rotate;
			}

			res->cont_map1 = append_new_map(res->cont_map1, map1->level, x, rotate%root_of_unit);

			return res;
		}

		//template <class LeftOperandNode, class RightOperandNode>
		Edge<mNode> cont2(const Edge<mNode>& x, const Edge<mNode>& y, key_2_new_key_node* key_2_new_key1, key_2_new_key_node* key_2_new_key2, const int var_num) {
			auto& id = this->identity;
			if (enableContStageTrace) {
				contStageTraceDepth++;
			}
			using ResultEdge = Edge<mNode>;
			const auto traceRootCall = enableContStageTrace && contStageTraceDepth == 1;
			const auto traceChildCall = enableContStageTrace && contStageTraceNextChild && contStageTraceDepth == contStageTraceNextChildParentDepth + 1;
			const auto traceStep = contStageTraceStep;
			const auto traceEdgeNodesBefore = traceRootCall ? size(x) : 0U;
			const auto traceStart = traceRootCall ? regressionDiagnostics : RegressionDiagnostics{};
			const auto traceChildBranch = traceChildCall && contStageTraceNextChildBranch ? contStageTraceNextChildBranch : "";
			const auto traceChildK = traceChildCall ? contStageTraceNextChildK : -1;
			if (traceChildCall) {
				contStageTraceNextChild = false;
				contStageTraceFocusedChildActive = true;
				contStageTraceNextChildParentDepth = 0;
				contStageTraceNextChildBranch = nullptr;
				contStageTraceNextChildK = -1;
			}
			struct ContStageTraceGuard {
				Package* pkg;
				bool focusedChild;
				~ContStageTraceGuard() {
					if (focusedChild) {
						pkg->contStageTraceFocusedChildActive = false;
					}
					if (pkg->enableContStageTrace && pkg->contStageTraceDepth > 0) {
						pkg->contStageTraceDepth--;
					}
				}
			} traceGuard{this, traceChildCall};
			auto printContStageDelta = [&](const char* stage,
									  const RegressionDiagnostics& before,
									  const RegressionDiagnostics& after,
									  const unsigned int nodesBefore,
									  const unsigned int nodesAfter) {
				if (!traceRootCall) {
					return;
				}
				auto appendDelta = [](std::ostream& os, const char* label, const std::size_t lhs, const std::size_t rhs) {
					if (rhs != lhs) {
						os << " " << label << "+=" << (rhs - lhs);
					}
				};
				std::ostringstream output;
				output << "cont_stage[" << traceStep << "]:" << stage << " nodes=" << nodesBefore << "->" << nodesAfter;
				appendDelta(output, "normalize.zero_children", before.normalizeZeroChildren, after.normalizeZeroChildren);
				appendDelta(output, "normalize.child_phase_adds", before.normalizeChildPhaseAdds, after.normalizeChildPhaseAdds);
				appendDelta(output, "normalize.root_phase_promotions", before.normalizeRootPhasePromotions, after.normalizeRootPhasePromotions);
				appendDelta(output, "tadd.same_pointer_map_mismatch", before.taddSamePointerMapMismatch, after.taddSamePointerMapMismatch);
				appendDelta(output, "tadd.mismatch_residual_header", before.taddMismatchResidualHeader, after.taddMismatchResidualHeader);
				appendDelta(output, "tadd.mismatch_residual_non_header", before.taddMismatchResidualNonHeader, after.taddMismatchResidualNonHeader);
				appendDelta(output, "tadd.mismatch_residual_phaseful", before.taddMismatchResidualPhaseful, after.taddMismatchResidualPhaseful);
				appendDelta(output, "tadd.mismatch_add_hits", before.taddMismatchAddHits, after.taddMismatchAddHits);
				appendDelta(output, "tadd.mismatch_add_misses", before.taddMismatchAddMisses, after.taddMismatchAddMisses);
				appendDelta(output, "tadd.mismatch_add_miss_empty", before.taddMismatchAddMissEmpty, after.taddMismatchAddMissEmpty);
				appendDelta(output, "tadd.mismatch_add_miss_map_only", before.taddMismatchAddMissMapOnly, after.taddMismatchAddMissMapOnly);
				appendDelta(output, "tadd.mismatch_add_miss_weight_only", before.taddMismatchAddMissWeightOnly, after.taddMismatchAddMissWeightOnly);
				appendDelta(output, "tadd.mismatch_add_miss_map_and_weight", before.taddMismatchAddMissMapAndWeight, after.taddMismatchAddMissMapAndWeight);
				appendDelta(output, "mapmul.lookup_phaseful", before.mapmulLookupPhaseful, after.mapmulLookupPhaseful);
				appendDelta(output, "mapmul.result_phaseful", before.mapmulResultPhaseful, after.mapmulResultPhaseful);
				appendDelta(output, "mapdiv.lookup_phaseful", before.mapdivLookupPhaseful, after.mapdivLookupPhaseful);
				appendDelta(output, "mapdiv.lookup_phase_overwrite", before.mapdivLookupPhaseOverwrite, after.mapdivLookupPhaseOverwrite);
				appendDelta(output, "mapdiv.lookup_phase_overwrite_non_header", before.mapdivLookupPhaseOverwriteNonHeader, after.mapdivLookupPhaseOverwriteNonHeader);
				appendDelta(output, "mapdiv.result_phaseful", before.mapdivResultPhaseful, after.mapdivResultPhaseful);
				appendDelta(output, "find_remain.phase_carries", before.findRemainPhaseCarries, after.findRemainPhaseCarries);
				std::cout << output.str() << std::endl;
			};
			auto printContCounterDelta = [&](const char* stage,
									 const RegressionDiagnostics& before,
									 const RegressionDiagnostics& after) {
				if (!traceRootCall) {
					return;
				}
				auto appendDelta = [](std::ostream& os, const char* label, const std::size_t lhs, const std::size_t rhs) {
					if (rhs != lhs) {
						os << " " << label << "+=" << (rhs - lhs);
					}
				};
				std::ostringstream output;
				output << "cont_stage[" << traceStep << "]:" << stage;
				appendDelta(output, "normalize.zero_children", before.normalizeZeroChildren, after.normalizeZeroChildren);
				appendDelta(output, "normalize.child_phase_adds", before.normalizeChildPhaseAdds, after.normalizeChildPhaseAdds);
				appendDelta(output, "normalize.root_phase_promotions", before.normalizeRootPhasePromotions, after.normalizeRootPhasePromotions);
				appendDelta(output, "tadd.same_pointer_map_mismatch", before.taddSamePointerMapMismatch, after.taddSamePointerMapMismatch);
				appendDelta(output, "tadd.mismatch_residual_header", before.taddMismatchResidualHeader, after.taddMismatchResidualHeader);
				appendDelta(output, "tadd.mismatch_residual_non_header", before.taddMismatchResidualNonHeader, after.taddMismatchResidualNonHeader);
				appendDelta(output, "tadd.mismatch_residual_phaseful", before.taddMismatchResidualPhaseful, after.taddMismatchResidualPhaseful);
				appendDelta(output, "tadd.mismatch_add_hits", before.taddMismatchAddHits, after.taddMismatchAddHits);
				appendDelta(output, "tadd.mismatch_add_misses", before.taddMismatchAddMisses, after.taddMismatchAddMisses);
				appendDelta(output, "tadd.mismatch_add_miss_empty", before.taddMismatchAddMissEmpty, after.taddMismatchAddMissEmpty);
				appendDelta(output, "tadd.mismatch_add_miss_map_only", before.taddMismatchAddMissMapOnly, after.taddMismatchAddMissMapOnly);
				appendDelta(output, "tadd.mismatch_add_miss_weight_only", before.taddMismatchAddMissWeightOnly, after.taddMismatchAddMissWeightOnly);
				appendDelta(output, "tadd.mismatch_add_miss_map_and_weight", before.taddMismatchAddMissMapAndWeight, after.taddMismatchAddMissMapAndWeight);
				appendDelta(output, "mapmul.lookup_phaseful", before.mapmulLookupPhaseful, after.mapmulLookupPhaseful);
				appendDelta(output, "mapmul.result_phaseful", before.mapmulResultPhaseful, after.mapmulResultPhaseful);
				appendDelta(output, "mapdiv.lookup_phaseful", before.mapdivLookupPhaseful, after.mapdivLookupPhaseful);
				appendDelta(output, "mapdiv.lookup_phase_overwrite", before.mapdivLookupPhaseOverwrite, after.mapdivLookupPhaseOverwrite);
				appendDelta(output, "mapdiv.lookup_phase_overwrite_non_header", before.mapdivLookupPhaseOverwriteNonHeader, after.mapdivLookupPhaseOverwriteNonHeader);
				appendDelta(output, "mapdiv.result_phaseful", before.mapdivResultPhaseful, after.mapdivResultPhaseful);
				appendDelta(output, "find_remain.phase_carries", before.findRemainPhaseCarries, after.findRemainPhaseCarries);
				std::cout << output.str() << std::endl;
			};
			auto throwOnNonFiniteResult = [&](const char* stage, const ResultEdge& edge) {
				if (!enableContStageTrace || contStageTraceDepth == 0) {
					return;
				}
				const auto real = CTEntry::val(edge.w.r);
				const auto imag = CTEntry::val(edge.w.i);
				if (std::isfinite(real) && std::isfinite(imag)) {
					return;
				}
				std::ostringstream message;
				message << "cont2 non-finite result"
					<< " step=" << traceStep
					<< " depth=" << contStageTraceDepth
					<< " stage=" << stage
					<< " var_num=" << var_num
					<< " node_var=" << static_cast<int>(edge.p->v)
					<< " weight=(" << real << ", " << imag << ")"
					<< " map_extra_phase=" << (edge.map ? edge.map->extra_phase : -999999);
				throw std::runtime_error(message.str());
			};
			bool traceSawRootMakeNode = false;
			auto traceBeforeRootMakeNode = RegressionDiagnostics{};
			auto traceAfterRootMakeNode = RegressionDiagnostics{};
				if (traceRootCall) {
					std::cerr << "cont_root_entry\tstep\t" << traceStep
						<< "\tvar_num\t" << var_num
						<< "\tx_w_re\t" << CTEntry::val(x.w.r)
					<< "\tx_w_im\t" << CTEntry::val(x.w.i)
					<< "\tx_var\t" << static_cast<int>(x.p->v)
					<< "\tx_map_level\t" << (x.map ? static_cast<int>(x.map->level) : -999)
					<< "\tx_map_x\t" << (x.map ? static_cast<int>(x.map->x) : -1)
					<< "\tx_map_rot\t" << (x.map ? x.map->rotate : -999)
					<< "\tx_map_ep\t" << (x.map ? x.map->extra_phase : -999)
					<< "\ty_w_re\t" << CTEntry::val(y.w.r)
					<< "\ty_w_im\t" << CTEntry::val(y.w.i)
					<< "\ty_var\t" << static_cast<int>(y.p->v)
					<< "\ty_map_level\t" << (y.map ? static_cast<int>(y.map->level) : -999)
					<< "\ty_map_x\t" << (y.map ? static_cast<int>(y.map->x) : -1)
					<< "\ty_map_rot\t" << (y.map ? y.map->rotate : -999)
					<< "\ty_map_ep\t" << (y.map ? y.map->extra_phase : -999)
					<< "\n";
			}
			if (traceChildCall) {
				std::cerr << "cont_child_entry\tstep\t" << traceStep
					<< "\tparent_branch\t" << traceChildBranch
					<< "\tparent_k\t" << traceChildK
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tvar_num\t" << var_num
					<< "\tx_w_re\t" << CTEntry::val(x.w.r)
					<< "\tx_w_im\t" << CTEntry::val(x.w.i)
					<< "\tx_var\t" << static_cast<int>(x.p->v)
					<< "\tx_map_level\t" << (x.map ? static_cast<int>(x.map->level) : -999)
					<< "\tx_map_x\t" << (x.map ? static_cast<int>(x.map->x) : -1)
					<< "\tx_map_rot\t" << (x.map ? x.map->rotate : -999)
					<< "\tx_map_ep\t" << (x.map ? x.map->extra_phase : -999)
					<< "\ty_w_re\t" << CTEntry::val(y.w.r)
					<< "\ty_w_im\t" << CTEntry::val(y.w.i)
					<< "\ty_var\t" << static_cast<int>(y.p->v)
					<< "\ty_map_level\t" << (y.map ? static_cast<int>(y.map->level) : -999)
					<< "\ty_map_x\t" << (y.map ? static_cast<int>(y.map->x) : -1)
					<< "\ty_map_rot\t" << (y.map ? y.map->rotate : -999)
					<< "\ty_map_ep\t" << (y.map ? y.map->extra_phase : -999)
					<< "\n";
			}
			//std::cout <<"838 " << x.w << " " << y.w.r->value<<" "<<y.w.i->value<< std::endl;
			//std::cout <<"838 " << x.w << " " << y.w << " " << int(x.p->v) << " " << int(y.p->v) << std::endl;
			//the_maps::print_maps(x.map);
			//the_maps::print_maps(y.map);

			if (x.p == nullptr) {
				return { nullptr, Complex::zero };
			}
			if (y.p == nullptr) {
				return y;
			}

			if (x.w.exactlyZero() || y.w.exactlyZero()) {
				return ResultEdge::zero;
			}


			if (x.p->v == -1 && y.p->v == -1)
			{
				auto c = cn.mulCached(x.w, y.w);
				const auto traceTerminalScale = enableContStageTrace && (traceStep == 519 || (contStageTraceDepth > 1 && var_num > 0));
				if (traceTerminalScale) {
					std::cerr << "cont_terminal_entry\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tvar_num\t" << var_num
						<< "\tx_w_re\t" << CTEntry::val(x.w.r)
						<< "\tx_w_im\t" << CTEntry::val(x.w.i)
						<< "\ty_w_re\t" << CTEntry::val(y.w.r)
						<< "\ty_w_im\t" << CTEntry::val(y.w.i)
						<< "\tc_pre_re\t" << CTEntry::val(c.r)
						<< "\tc_pre_im\t" << CTEntry::val(c.i)
						<< "\n";
				}

				if (var_num > 0 && std::getenv("LIMTDD_DISABLE_TERMINAL_VAR_SCALE") == nullptr) {
					assert(c != Complex::zero);
					ComplexNumbers::mul(c, c, cn.getTemporary(pow(2, var_num), 0));
				}
				if (traceTerminalScale) {
					std::cerr << "cont_terminal_return\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tvar_num\t" << var_num
						<< "\tpow2\t" << (var_num > 0 ? pow(2, var_num) : 1)
						<< "\tc_re\t" << CTEntry::val(c.r)
						<< "\tc_im\t" << CTEntry::val(c.i)
						<< "\n";
				}
				//std::cout << "Case 00" << std::endl;
				return ResultEdge::terminal(c);
			}

			key_2_new_key_node* temp_key_2_new_key2 = key_2_new_key2;
			while (temp_key_2_new_key2->level > y.p->v) {
				temp_key_2_new_key2 = temp_key_2_new_key2->father;
			}


			if (x.p->v == -1 && var_num == 0 && std::abs(temp_key_2_new_key2->new_key - y.p->v) < 1e-10) {
				//std::cout << "Case 01" << std::endl;
				return 	ResultEdge{ y.p, cn.mulCached(x.w, y.w) ,y.map};
			}

			key_2_new_key_node* temp_key_2_new_key1 = key_2_new_key1;
			while (temp_key_2_new_key1->level > x.p->v) {
				temp_key_2_new_key1 = temp_key_2_new_key1->father;
			}

			if (y.p->v == -1 && var_num == 0 && std::abs(temp_key_2_new_key1->new_key - x.p->v) < 1e-10) {
				//std::cout << "Case 02" << std::endl;
				return 	ResultEdge{ x.p, cn.mulCached(x.w, y.w) ,x.map};
			}


			auto xCopy = x;
			xCopy.w = Complex::one;
			auto yCopy = y;
			yCopy.w = Complex::one;

			
			auto r_maps = find_remain_map(x.map, y.map, key_2_new_key1, key_2_new_key2);
			if (traceRootCall) {
				traceMapChain("cont_root_x_map_chain", x.map);
				traceMapChain("cont_root_y_map_chain", y.map);
				traceMapChain("cont_root_cont1_chain", r_maps->cont_map1);
				traceMapChain("cont_root_cont2_chain", r_maps->cont_map2);
				traceMapChain("cont_root_remain_chain", r_maps->remain_map);
			}
			if (traceChildCall) {
				traceMapChain("cont_child_x_map_chain", x.map);
				traceMapChain("cont_child_y_map_chain", y.map);
				traceMapChain("cont_child_cont1_chain", r_maps->cont_map1);
				traceMapChain("cont_child_cont2_chain", r_maps->cont_map2);
				traceMapChain("cont_child_remain_chain", r_maps->remain_map);
			}
			const auto traceAfterFindRemain = traceRootCall ? regressionDiagnostics : RegressionDiagnostics{};
			if (traceRootCall) {
				printContStageDelta("find_remain", traceStart, traceAfterFindRemain, traceEdgeNodesBefore, traceEdgeNodesBefore);
			}

			xCopy.map = r_maps->cont_map1;
			yCopy.map = r_maps->cont_map2;
			// yCopy.map->print_maps(yCopy.map);
			//auto extra_phase = cn.getCached(r_maps->remain_map->extra_phase.r->value, r_maps->remain_map->extra_phase.i->value);
			auto extra_phase = r_maps->remain_map->extra_phase;
			if (traceChildCall) {
				std::cerr << "cont_child_maps\tstep\t" << traceStep
					<< "\tparent_branch\t" << traceChildBranch
					<< "\tparent_k\t" << traceChildK
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tcont1_level\t" << (r_maps->cont_map1 ? static_cast<int>(r_maps->cont_map1->level) : -999)
					<< "\tcont1_x\t" << (r_maps->cont_map1 ? static_cast<int>(r_maps->cont_map1->x) : -1)
					<< "\tcont1_rot\t" << (r_maps->cont_map1 ? r_maps->cont_map1->rotate : -999)
					<< "\tcont1_ep\t" << (r_maps->cont_map1 ? r_maps->cont_map1->extra_phase : -999)
					<< "\tcont2_level\t" << (r_maps->cont_map2 ? static_cast<int>(r_maps->cont_map2->level) : -999)
					<< "\tcont2_x\t" << (r_maps->cont_map2 ? static_cast<int>(r_maps->cont_map2->x) : -1)
					<< "\tcont2_rot\t" << (r_maps->cont_map2 ? r_maps->cont_map2->rotate : -999)
					<< "\tcont2_ep\t" << (r_maps->cont_map2 ? r_maps->cont_map2->extra_phase : -999)
					<< "\tremain_level\t" << (r_maps->remain_map ? static_cast<int>(r_maps->remain_map->level) : -999)
					<< "\tremain_x\t" << (r_maps->remain_map ? static_cast<int>(r_maps->remain_map->x) : -1)
					<< "\tremain_rot\t" << (r_maps->remain_map ? r_maps->remain_map->rotate : -999)
					<< "\tremain_ep\t" << (r_maps->remain_map ? r_maps->remain_map->extra_phase : -999)
					<< "\textra_phase\t" << extra_phase
					<< "\n";
			}
			if (traceRootCall) {
				std::cerr << "cont_root_maps\tstep\t" << traceStep
					<< "\tcont1_level\t" << (r_maps->cont_map1 ? static_cast<int>(r_maps->cont_map1->level) : -999)
					<< "\tcont1_x\t" << (r_maps->cont_map1 ? static_cast<int>(r_maps->cont_map1->x) : -1)
					<< "\tcont1_rot\t" << (r_maps->cont_map1 ? r_maps->cont_map1->rotate : -999)
					<< "\tcont1_ep\t" << (r_maps->cont_map1 ? r_maps->cont_map1->extra_phase : -999)
					<< "\tcont2_level\t" << (r_maps->cont_map2 ? static_cast<int>(r_maps->cont_map2->level) : -999)
					<< "\tcont2_x\t" << (r_maps->cont_map2 ? static_cast<int>(r_maps->cont_map2->x) : -1)
					<< "\tcont2_rot\t" << (r_maps->cont_map2 ? r_maps->cont_map2->rotate : -999)
					<< "\tcont2_ep\t" << (r_maps->cont_map2 ? r_maps->cont_map2->extra_phase : -999)
					<< "\tremain_level\t" << (r_maps->remain_map ? static_cast<int>(r_maps->remain_map->level) : -999)
					<< "\tremain_x\t" << (r_maps->remain_map ? static_cast<int>(r_maps->remain_map->x) : -1)
					<< "\tremain_rot\t" << (r_maps->remain_map ? r_maps->remain_map->rotate : -999)
					<< "\tremain_ep\t" << (r_maps->remain_map ? r_maps->remain_map->extra_phase : -999)
					<< "\textra_phase\t" << extra_phase
					<< "\n";
			}

			auto res = (std::getenv("LIMTDD_DISABLE_CONT_CACHE") == nullptr)
				? contTable.lookup(xCopy, yCopy, temp_key_2_new_key1, temp_key_2_new_key2)
				: decltype(contTable.lookup(xCopy, yCopy, temp_key_2_new_key1, temp_key_2_new_key2)){};
			if (std::getenv("LIMTDD_BYPASS_REMAIN_MAP_CONT_CACHE") != nullptr && res.e.p != nullptr &&
				((x.map != xCopy.map && xCopy.map != nullptr && xCopy.map->level == -1) ||
				 (y.map != yCopy.map && yCopy.map != nullptr && yCopy.map->level == -1))) {
				res = {};
			}
			if (res.e.p != nullptr) {
				if (enableContStageTrace && traceStep == 519) {
					std::cerr << "cont_cache_hit\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tcont_num\t" << res.cont_num
						<< "\tvar_num\t" << var_num
						<< "\tx_w_re\t" << CTEntry::val(x.w.r)
						<< "\tx_w_im\t" << CTEntry::val(x.w.i)
						<< "\tx_var\t" << static_cast<int>(x.p->v)
						<< "\tx_map_level\t" << (x.map ? static_cast<int>(x.map->level) : -999)
						<< "\tx_copy_map_level\t" << (xCopy.map ? static_cast<int>(xCopy.map->level) : -999)
						<< "\tx_copy_map_x\t" << (xCopy.map ? static_cast<int>(xCopy.map->x) : -1)
						<< "\tx_copy_map_rot\t" << (xCopy.map ? xCopy.map->rotate : -999)
						<< "\tx_copy_map_ep\t" << (xCopy.map ? xCopy.map->extra_phase : -999)
						<< "\ty_w_re\t" << CTEntry::val(y.w.r)
						<< "\ty_w_im\t" << CTEntry::val(y.w.i)
						<< "\ty_var\t" << static_cast<int>(y.p->v)
						<< "\ty_map_level\t" << (y.map ? static_cast<int>(y.map->level) : -999)
						<< "\tres_w_re\t" << res.e.w.r
						<< "\tres_w_im\t" << res.e.w.i
						<< "\tres_var\t" << static_cast<int>(res.e.p->v)
						<< "\tres_map_level\t" << (res.e.map ? static_cast<int>(res.e.map->level) : -999)
						<< "\n";
				}
				if (traceRootCall) {
					std::cerr << "cont_root_cache\tstep\t" << traceStep
						<< "\tcont_num\t" << res.cont_num
						<< "\tvar_num\t" << var_num
						<< "\tres_w_re\t" << res.e.w.r
						<< "\tres_w_im\t" << res.e.w.i
						<< "\n";
				}
				if (res.e.w.approximatelyZero()) {
					// cn.returnToCache(extra_phase);
					return ResultEdge::zero;
				}
				auto e = ResultEdge{ res.e.p, cn.getCached(res.e.w),res.e.map };
				assert(e.w != Complex::zero);
				if (enableContStageTrace && traceStep == 519) {
					std::cerr << "cont_cache_return_stage\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tstage\tloaded"
						<< "\tvar_num\t" << var_num
						<< "\tcont_num\t" << res.cont_num
						<< "\te_w_re\t" << CTEntry::val(e.w.r)
						<< "\te_w_im\t" << CTEntry::val(e.w.i)
						<< "\te_map_level\t" << (e.map ? static_cast<int>(e.map->level) : -999)
						<< "\n";
				}
				ComplexNumbers::mul(e.w, e.w, x.w);
				ComplexNumbers::mul(e.w, e.w, y.w);
				if (enableContStageTrace && traceStep == 519) {
					std::cerr << "cont_cache_return_stage\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tstage\tafter_xy"
						<< "\tvar_num\t" << var_num
						<< "\tcont_num\t" << res.cont_num
						<< "\te_w_re\t" << CTEntry::val(e.w.r)
						<< "\te_w_im\t" << CTEntry::val(e.w.i)
						<< "\n";
				}
				if (e.w.approximatelyZero()) {
					//assert(e.w != Complex::zero);
					// cn.returnToCache(e.w);
					// cn.returnToCache(extra_phase);
					return ResultEdge::zero;
				}
				//std::cout << "1160 " << var_num << " " << res.cont_num << std::endl;
				if (res.cont_num != var_num && std::getenv("LIMTDD_DISABLE_CACHE_CONT_NUM_SCALE") == nullptr) {
					assert(e.w != Complex::zero);
					ComplexNumbers::mul(e.w, e.w, cn.getTemporary(pow(2, var_num - res.cont_num), 0));//对于一般形状的tensor,以2为底数可能有问题
					// TODO: pow(2,n) can be optimized by 1<<n
				}
				if (enableContStageTrace && traceStep == 519) {
					std::cerr << "cont_cache_return_stage\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tstage\tafter_cont_num"
						<< "\tvar_num\t" << var_num
						<< "\tcont_num\t" << res.cont_num
						<< "\tpow2\t" << (res.cont_num != var_num ? pow(2, var_num - res.cont_num) : 1)
						<< "\te_w_re\t" << CTEntry::val(e.w.r)
						<< "\te_w_im\t" << CTEntry::val(e.w.i)
						<< "\n";
				}
				e.map = mapmul(r_maps->remain_map, e.map);
				assert(e.w != Complex::zero);
				// cn.mul(e.w, e.w, e.map->extra_phase);
				cn.mul(e.w, e.w, cn.getTemporary(cos(e.map->extra_phase*rotate_angle),sin(e.map->extra_phase*rotate_angle)));
				if (enableContStageTrace && traceStep == 519) {
					std::cerr << "cont_cache_return_stage\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tstage\tafter_map_phase"
						<< "\te_w_re\t" << CTEntry::val(e.w.r)
						<< "\te_w_im\t" << CTEntry::val(e.w.i)
						<< "\te_map_level\t" << (e.map ? static_cast<int>(e.map->level) : -999)
						<< "\te_map_ep\t" << (e.map ? e.map->extra_phase : -999)
						<< "\n";
				}
				// cn.returnToCache(e.map->extra_phase);
				assert(e.w != Complex::zero);
				// cn.mul(e.w, e.w, extra_phase);
				cn.mul(e.w, e.w, cn.getTemporary(cos(extra_phase*rotate_angle),sin(extra_phase*rotate_angle)));
				if (enableContStageTrace && traceStep == 519) {
					std::cerr << "cont_cache_return_stage\tstep\t" << traceStep
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tstage\tafter_extra_phase"
						<< "\textra_phase\t" << extra_phase
						<< "\te_w_re\t" << CTEntry::val(e.w.r)
						<< "\te_w_im\t" << CTEntry::val(e.w.i)
						<< "\n";
				}
				// cn.returnToCache(extra_phase);
				if (traceRootCall) {
					std::cerr << "cont_root_cache_return\tstep\t" << traceStep
						<< "\te_w_re\t" << CTEntry::val(e.w.r)
						<< "\te_w_im\t" << CTEntry::val(e.w.i)
						<< "\n";
				}
				return e;
			}
			// TODO: add if here
			/*
			check if y is identity.
			(identity == yCopy)
			{ auto y0 = varOrder(yCopy.first)
			auto y1 = varOrder(yCopy.second)
			auto x = varOrder(yCopy.first.father)
			auto z = varOrder(yCopy.first.next)
			if(x>y1 and y1 > z){
				return xCopy
			}
			else{
				continue ?
			}
			}
			*/
			// bool test_1314 = false;
			if(yCopy == this->identity){
				if(to_test){
					std::cout << "a identity" << std::endl;
				}
				auto temp_temp_key_2_new_key1 = temp_key_2_new_key1;
				bool flag = true ;
				int i = 0;
				while (temp_temp_key_2_new_key1->level > -1) {
					if (ifContract(temp_temp_key_2_new_key1->new_key)){
						temp_temp_key_2_new_key1 = temp_temp_key_2_new_key1->father;
						++ i;
					}
					else{
						if(temp_temp_key_2_new_key1->level != temp_temp_key_2_new_key1->new_key){
							flag = false;
							break;
						}
						else{
							temp_temp_key_2_new_key1 = temp_temp_key_2_new_key1->father;
						}
					}
				}

				// if(flag){
				// 	auto e = ResultEdge{ xCopy.p, cn.getCached(1,0),xCopy.map };
				// 	assert(e.w != Complex::zero);
				// 	ComplexNumbers::mul(e.w, e.w, x.w);
				// 	ComplexNumbers::mul(e.w, e.w, y.w);
				// 	if (e.w.approximatelyZero()) {
				// 		//assert(e.w != Complex::zero);
				// 		cn.returnToCache(e.w);
				// 		cn.returnToCache(extra_phase);
				// 		return ResultEdge::zero;
				// 	}
				// 	e.map = mapmul(r_maps->remain_map, e.map);
				// 	assert(e.w != Complex::zero);
				// 	cn.mul(e.w, e.w, e.map->extra_phase);
				// 	cn.returnToCache(e.map->extra_phase);
				// 	assert(e.w != Complex::zero);
				// 	cn.mul(e.w, e.w, extra_phase);
				// 	cn.returnToCache(extra_phase);
				// 	std::cout << i << std::endl;
				// 	return e;
				// }
				// test_1314 = flag;
			}


			float newk1 = temp_key_2_new_key1->new_key;

			float newk2 = temp_key_2_new_key2->new_key;
			const auto* focusedContNewkEnv = std::getenv("LIMTDD_FOCUSED_CONT_NEWK1");
			const auto* focusedContDepthEnv = std::getenv("LIMTDD_FOCUSED_CONT_DEPTH");
			const auto traceFocusedContByNewk = focusedContNewkEnv != nullptr &&
				(std::abs(newk1 - std::stof(focusedContNewkEnv)) < 1e-6 || std::abs(newk2 - std::stof(focusedContNewkEnv)) < 1e-6);
			const auto traceFocusedContByDepth = focusedContDepthEnv != nullptr &&
				contStageTraceDepth == static_cast<std::size_t>(std::stoull(focusedContDepthEnv));
			const auto traceFocusedContCall = enableContStageTrace && contStageTraceDepth > 1 &&
				(traceFocusedContByNewk || traceFocusedContByDepth);
			if (traceFocusedContCall) {
				std::cerr << "cont_focused_keys\tstep\t" << traceStep
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tvar_num\t" << var_num
					<< "\tnewk1\t" << newk1
					<< "\tnewk2\t" << newk2
					<< "\tifc1\t" << ifContract(newk1)
					<< "\tifc2\t" << ifContract(newk2)
					<< "\tx_var\t" << static_cast<int>(x.p->v)
					<< "\ty_var\t" << static_cast<int>(y.p->v)
					<< "\tx_w_re\t" << CTEntry::val(x.w.r)
					<< "\ty_w_re\t" << CTEntry::val(y.w.r)
					<< "\n";
			}
			if (traceRootCall) {
				std::cerr << "cont_root_keys\tstep\t" << traceStep
					<< "\tvar_num\t" << var_num
					<< "\tnewk1\t" << newk1
					<< "\tnewk2\t" << newk2
					<< "\tifc1\t" << ifContract(newk1)
					<< "\tifc2\t" << ifContract(newk2)
					<< "\tx_var\t" << static_cast<int>(x.p->v)
					<< "\ty_var\t" << static_cast<int>(y.p->v)
					<< "\n";
			}
			if (traceChildCall) {
				std::cerr << "cont_child_keys\tstep\t" << traceStep
					<< "\tparent_branch\t" << traceChildBranch
					<< "\tparent_k\t" << traceChildK
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tvar_num\t" << var_num
					<< "\tnewk1\t" << newk1
					<< "\tnewk2\t" << newk2
					<< "\tifc1\t" << ifContract(newk1)
					<< "\tifc2\t" << ifContract(newk2)
					<< "\tx_var\t" << static_cast<int>(x.p->v)
					<< "\ty_var\t" << static_cast<int>(y.p->v)
					<< "\n";
			}
			if(to_test){
				std::cout << 1298 << std::endl;
				std::cout << "newk1: " << newk1 << " newk2: " << newk2 << std::endl;
			}
			ResultEdge r;
			auto traceRootAggEdge = [&](const char* branch, const char* stage, const int k, const ResultEdge& edge) {
				if (!traceRootCall && !traceChildCall && !traceFocusedContCall) {
					return;
				}
				std::cerr << (traceChildCall ? "cont_child_agg" : (traceFocusedContCall ? "cont_focused_agg" : "cont_root_agg")) << "\tstep\t" << traceStep;
				if (traceChildCall) {
					std::cerr << "\tparent_branch\t" << traceChildBranch
						<< "\tparent_k\t" << traceChildK
						<< "\tdepth\t" << contStageTraceDepth;
				} else if (traceFocusedContCall) {
					std::cerr << "\tdepth\t" << contStageTraceDepth;
				}
				std::cerr << "\tbranch\t" << branch
					<< "\tstage\t" << stage
					<< "\tk\t" << k
					<< "\tvar_num\t" << var_num
					<< "\tnewk1\t" << newk1
					<< "\tnewk2\t" << newk2
					<< "\tw_re\t" << CTEntry::val(edge.w.r)
					<< "\tw_im\t" << CTEntry::val(edge.w.i)
					<< "\tchild_var\t" << (edge.p ? static_cast<int>(edge.p->v) : -999)
					<< "\tmap_level\t" << (edge.map ? static_cast<int>(edge.map->level) : -999)
					<< "\tmap_x\t" << (edge.map ? static_cast<int>(edge.map->x) : -1)
					<< "\tmap_rot\t" << (edge.map ? edge.map->rotate : -999)
					<< "\tmap_ep\t" << (edge.map ? edge.map->extra_phase : -999)
					<< "\n";
			};
			auto traceRootAggAction = [&](const char* branch, const char* action, const int k) {
				if (!traceRootCall && !traceChildCall && !traceFocusedContCall) {
					return;
				}
				std::cerr << (traceChildCall ? "cont_child_agg_action" : (traceFocusedContCall ? "cont_focused_agg_action" : "cont_root_agg_action")) << "\tstep\t" << traceStep;
				if (traceChildCall) {
					std::cerr << "\tparent_branch\t" << traceChildBranch
						<< "\tparent_k\t" << traceChildK
						<< "\tdepth\t" << contStageTraceDepth;
				} else if (traceFocusedContCall) {
					std::cerr << "\tdepth\t" << contStageTraceDepth;
				}
				std::cerr << "\tbranch\t" << branch
					<< "\taction\t" << action
					<< "\tk\t" << k
					<< "\tvar_num\t" << var_num
					<< "\tnewk1\t" << newk1
					<< "\tnewk2\t" << newk2
					<< "\n";
			};

			if (newk1 > newk2) {
				// TODO: half integer?
				// bool isHalfInteger = std::abs(newk1 - std::round(newk1)) > 0.4999999 && std::abs(newk1 - std::round(newk1)) < 0.5000001;
				if (ifContract(newk1)) {
					r = ResultEdge::zero;
					ResultEdge etemp;
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						auto e1 = Slicing2(xCopy, xCopy.p->v, k);
						auto& e2 = yCopy;
						etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
						traceRootAggEdge("gt_contract", "input_e1", k, e1);
						traceRootAggEdge("gt_contract", "input_e2", k, e2);
						traceRootAggEdge("gt_contract", "etemp", k, etemp);
						throwOnNonFiniteResult("recurse_gt", etemp);
						if (e1.w != Complex::zero) {
							// cn.returnToCache(e1.w);
						}
						if (etemp.w != Complex::zero) {
							traceRootAggEdge("gt_contract", "r_before", k, r);
							if (r != ResultEdge::zero) {
								traceRootAggAction("gt_contract", "T_add2", k);
								auto temp = r.w;
								r = T_add2(r, etemp);
								traceRootAggEdge("gt_contract", "r_after_add", k, r);
								//assert(temp != Complex::zero);
								//assert(etemp.w != Complex::zero);
								// cn.returnToCache(temp);
								// cn.returnToCache(etemp.w);
							}
							else {
								traceRootAggAction("gt_contract", "assign", k);
								r = etemp;
								traceRootAggEdge("gt_contract", "r_after_assign", k, r);
							}
						}
					}
				}
				else {
					std::vector<ResultEdge> e;
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						auto e1 = Slicing2(xCopy, xCopy.p->v, k);
						auto& e2 = yCopy;
						if (traceRootCall) {
							std::cerr << "cont_root_child_input\tstep\t" << traceStep
								<< "\tbranch\tgt"
								<< "\tk\t" << k
								<< "\tvar_num\t" << var_num
								<< "\te1_w_re\t" << CTEntry::val(e1.w.r)
								<< "\te1_w_im\t" << CTEntry::val(e1.w.i)
								<< "\te1_var\t" << static_cast<int>(e1.p->v)
								<< "\te1_map_level\t" << (e1.map ? static_cast<int>(e1.map->level) : -999)
								<< "\te1_map_x\t" << (e1.map ? static_cast<int>(e1.map->x) : -1)
								<< "\te1_map_rot\t" << (e1.map ? e1.map->rotate : -999)
								<< "\te1_map_ep\t" << (e1.map ? e1.map->extra_phase : -999)
								<< "\te2_w_re\t" << CTEntry::val(e2.w.r)
								<< "\te2_w_im\t" << CTEntry::val(e2.w.i)
								<< "\te2_var\t" << static_cast<int>(e2.p->v)
								<< "\te2_map_level\t" << (e2.map ? static_cast<int>(e2.map->level) : -999)
								<< "\te2_map_x\t" << (e2.map ? static_cast<int>(e2.map->x) : -1)
								<< "\te2_map_rot\t" << (e2.map ? e2.map->rotate : -999)
								<< "\te2_map_ep\t" << (e2.map ? e2.map->extra_phase : -999)
								<< "\n";
						}
						if (traceRootCall && k == 0) {
							contStageTraceNextChild = true;
							contStageTraceNextChildParentDepth = contStageTraceDepth;
							contStageTraceNextChildBranch = "gt";
							contStageTraceNextChildK = k;
						}
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						if (traceRootCall) {
							const auto& child = e.back();
							std::cerr << "cont_root_child_output\tstep\t" << traceStep
								<< "\tbranch\tgt"
								<< "\tk\t" << k
								<< "\tchild_w_re\t" << CTEntry::val(child.w.r)
								<< "\tchild_w_im\t" << CTEntry::val(child.w.i)
								<< "\tchild_var\t" << static_cast<int>(child.p->v)
								<< "\tchild_map_level\t" << (child.map ? static_cast<int>(child.map->level) : -999)
								<< "\tchild_map_x\t" << (child.map ? static_cast<int>(child.map->x) : -1)
								<< "\tchild_map_rot\t" << (child.map ? child.map->rotate : -999)
								<< "\tchild_map_ep\t" << (child.map ? child.map->extra_phase : -999)
								<< "\n";
						}
						throwOnNonFiniteResult("vector_gt", e.back());
						if (e1.w != Complex::zero) {
							// cn.returnToCache(e1.w);
						}
					}
					if(to_test){
					std::cout <<"cont2 function: " << 1342 << std::endl; 
					std::cout << "var: " << newk1 << std::endl;
					std::cout << "var(newk2): " << newk2 << std::endl;
					std::cout << "edge.node.key: " << std::endl;
					for(auto edge:e){
						std::cout << edge.p->v << " " ;
						if(edge.p->v == newk1) throw std::runtime_error("bug here");
					}
					std::cout << std::endl;

				}
					if (traceRootCall || traceFocusedContCall) {
						traceBeforeRootMakeNode = regressionDiagnostics;
						for (std::size_t edgeIndex = 0; edgeIndex < e.size(); ++edgeIndex) {
							std::cerr << (traceFocusedContCall ? "cont_focused_edges" : "cont_root_edges") << "\tstep\t" << traceStep
								<< "\tedge\t" << edgeIndex
								<< "\tw_re\t" << CTEntry::val(e[edgeIndex].w.r)
								<< "\tw_im\t" << CTEntry::val(e[edgeIndex].w.i)
								<< "\tchild_var\t" << static_cast<int>(e[edgeIndex].p->v)
								<< "\tmap_level\t" << (e[edgeIndex].map ? static_cast<int>(e[edgeIndex].map->level) : -999)
								<< "\tmap_x\t" << (e[edgeIndex].map ? static_cast<int>(e[edgeIndex].map->x) : -1)
								<< "\tmap_rot\t" << (e[edgeIndex].map ? e[edgeIndex].map->rotate : -999)
								<< "\tmap_ep\t" << (e[edgeIndex].map ? e[edgeIndex].map->extra_phase : -999)
								<< "\n";
						}
					}
					r = makeDDNode(Qubit(newk1), e, true);
					if (traceRootCall) {
						std::cerr << "cont_root_after_make\tstep\t" << traceStep
							<< "\tbranch\tgt"
							<< "\tr_w_re\t" << CTEntry::val(r.w.r)
							<< "\tr_w_im\t" << CTEntry::val(r.w.i)
							<< "\tr_var\t" << static_cast<int>(r.p->v)
							<< "\tr_map_level\t" << (r.map ? static_cast<int>(r.map->level) : -999)
							<< "\tr_map_x\t" << (r.map ? static_cast<int>(r.map->x) : -1)
							<< "\tr_map_rot\t" << (r.map ? r.map->rotate : -999)
							<< "\tr_map_ep\t" << (r.map ? r.map->extra_phase : -999)
							<< "\n";
						traceSawRootMakeNode = true;
						traceAfterRootMakeNode = regressionDiagnostics;
					}
				}
			}
			else if (newk1 < newk2) {
				if (ifContract(newk2)) {
					r = ResultEdge::zero;
					ResultEdge etemp;
					for (int k = 0; k < y.p->e.size(); ++k) {
						auto& e1 = xCopy;
						//e2 = y.p->e[k];
						auto e2 = Slicing2(yCopy, yCopy.p->v, k);
							etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
							traceRootAggEdge("lt_contract", "input_e1", k, e1);
							traceRootAggEdge("lt_contract", "input_e2", k, e2);
							traceRootAggEdge("lt_contract", "etemp", k, etemp);
							throwOnNonFiniteResult("recurse_lt", etemp);
						if (e2.w != Complex::zero) {
							// cn.returnToCache(e2.w);
						}
						if (etemp.w != Complex::zero) {
							traceRootAggEdge("lt_contract", "r_before", k, r);
							if (r != ResultEdge::zero) {
								traceRootAggAction("lt_contract", "T_add2", k);
								auto temp = r.w;
								r = T_add2(r, etemp);
								traceRootAggEdge("lt_contract", "r_after_add", k, r);
								//assert(temp != Complex::zero);
								//assert(etemp.w != Complex::zero);
								// cn.returnToCache(temp);
								// cn.returnToCache(etemp.w);
							}
							else {
								traceRootAggAction("lt_contract", "assign", k);
								r = etemp;
								traceRootAggEdge("lt_contract", "r_after_assign", k, r);
							}
						}
					}
				}
				else {
					std::vector<ResultEdge> e;
					for (int k = 0; k < y.p->e.size(); ++k) {
						auto& e1 = xCopy;
						//e2 = y.p->e[k];
						auto e2 = Slicing2(yCopy, yCopy.p->v, k);
						traceRootAggEdge("lt", "input_e1", k, e1);
						traceRootAggEdge("lt", "input_e2", k, e2);
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						traceRootAggEdge("lt", "child_output", k, e.back());
						throwOnNonFiniteResult("vector_lt", e.back());
						if (e2.w != Complex::zero) {
							// cn.returnToCache(e2.w);
						}
					}
					if(to_test){
					std::cout << "cont2 function: " << 1391 << std::endl; 
					std::cout << "var: " << newk2 << std::endl;
					std::cout << "edge.node.key: " << std::endl;
					for(auto edge:e){
						std::cout << edge.p->v << " " ;
						if(edge.p->v == newk2) throw std::runtime_error("bug here");
					}
					std::cout << std::endl;

				}
					if (traceRootCall || traceFocusedContCall) {
						traceBeforeRootMakeNode = regressionDiagnostics;
						for (std::size_t edgeIndex = 0; edgeIndex < e.size(); ++edgeIndex) {
							std::cerr << (traceFocusedContCall ? "cont_focused_edges" : "cont_root_edges") << "\tstep\t" << traceStep
								<< "\tedge\t" << edgeIndex
								<< "\tw_re\t" << CTEntry::val(e[edgeIndex].w.r)
								<< "\tw_im\t" << CTEntry::val(e[edgeIndex].w.i)
								<< "\tchild_var\t" << static_cast<int>(e[edgeIndex].p->v)
								<< "\tmap_level\t" << (e[edgeIndex].map ? static_cast<int>(e[edgeIndex].map->level) : -999)
								<< "\tmap_x\t" << (e[edgeIndex].map ? static_cast<int>(e[edgeIndex].map->x) : -1)
								<< "\tmap_rot\t" << (e[edgeIndex].map ? e[edgeIndex].map->rotate : -999)
								<< "\tmap_ep\t" << (e[edgeIndex].map ? e[edgeIndex].map->extra_phase : -999)
								<< "\n";
						}
					}
					r = makeDDNode(Qubit(newk2), e, true);
					if (traceRootCall || traceFocusedContCall) {
						std::cerr << (traceFocusedContCall ? "cont_focused_after_make" : "cont_root_after_make") << "\tstep\t" << traceStep
							<< "\tbranch\tlt"
							<< "\tr_w_re\t" << CTEntry::val(r.w.r)
							<< "\tr_w_im\t" << CTEntry::val(r.w.i)
							<< "\tr_var\t" << static_cast<int>(r.p->v)
							<< "\tr_map_level\t" << (r.map ? static_cast<int>(r.map->level) : -999)
							<< "\tr_map_x\t" << (r.map ? static_cast<int>(r.map->x) : -1)
							<< "\tr_map_rot\t" << (r.map ? r.map->rotate : -999)
							<< "\tr_map_ep\t" << (r.map ? r.map->extra_phase : -999)
							<< "\n";
					}
					if (enableTailCxRenormExperiment &&
						ifContract(newk1) && !ifContract(newk2) &&
						var_num == 2 &&
						r.p != nullptr && static_cast<float>(r.p->v) == newk2 &&
						e.size() == 2 &&
						e[0].w.exactlyOne() && e[1].w.approximatelyZero() &&
						r.w.exactlyOne()) {
						ComplexNumbers::mul(r.w, r.w, cn.getTemporary(std::sqrt(0.5), 0));
					}
				}

			}
			else {
				if(to_test){
				std:: cout << 1410 << int(newk2 * 2) % 2 << std::endl; 

				}
				if (ifContract(newk2)) {
					r = ResultEdge::zero;
					ResultEdge etemp;
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						auto e1 = Slicing2(xCopy, xCopy.p->v, k);
						//e2 = y.p->e[k];
						//std::cout << "1554, e1 " << e1.w << std::endl;
						//std::cout << "e1.w in: " << (e1.w.i) << " " << e1.w.r << std::endl;
						// the_maps::print_maps(e1.map);
						auto e2 = Slicing2(yCopy, yCopy.p->v, k);
						//std::cout << "1556, e1 " << e1.w << std::endl;
						//std::cout << "e1.w in: " << (e1.w.i) << " " << e1.w.r << std::endl;
						//std::cout << "1558, e2 " << e2.w << std::endl;
						//std::cout << "e2.w in: " << (e2.w.i) << " " << e2.w.r << std::endl;
						//the_maps::print_maps(e1.map);
						etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
						traceRootAggEdge("eq_contract", "input_e1", k, e1);
						traceRootAggEdge("eq_contract", "input_e2", k, e2);
						traceRootAggEdge("eq_contract", "etemp", k, etemp);
						throwOnNonFiniteResult("recurse_eq", etemp);
						if (e1.w != Complex::zero) {
							// cn.returnToCache(e1.w);
						}
						if (e2.w != Complex::zero) {
							// cn.returnToCache(e2.w);
						}
						if (etemp.w != Complex::zero) {
							traceRootAggEdge("eq_contract", "r_before", k, r);
							if (r != ResultEdge::zero) {
								traceRootAggAction("eq_contract", "T_add2", k);
								auto temp = r.w;
								r = T_add2(r, etemp);
								traceRootAggEdge("eq_contract", "r_after_add", k, r);
								//assert(temp != Complex::zero);
								//assert(etemp.w != Complex::zero);
								// cn.returnToCache(temp);
								// cn.returnToCache(etemp.w);
							}
							else {
								traceRootAggAction("eq_contract", "assign", k);
								r = etemp;
								traceRootAggEdge("eq_contract", "r_after_assign", k, r);
							}
						}
					}
				}
				else {
					std::vector<ResultEdge> e;
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						auto e1 = Slicing2(xCopy, xCopy.p->v, k);
						//e2 = y.p->e[k];
						auto e2 = Slicing2(yCopy, yCopy.p->v, k);
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						throwOnNonFiniteResult("vector_eq", e.back());
						if (e1.w != Complex::zero) {
							// cn.returnToCache(e1.w);
						}
						if (e2.w != Complex::zero) {
							// cn.returnToCache(e2.w);
						}
					}
					if(to_test){
						std::cout <<"cont2 function: " << 1456 << std::endl; 
						std::cout << "var: " << newk1 << std::endl;
						std::cout << "var(newk2): " << newk2 << std::endl;
						std::cout << "edge.node.key: " << std::endl;
						for(auto edge:e){
							std::cout << edge.p->v << " " ;
							if(edge.p->v == newk1) throw std::runtime_error("bug here");
						}
						std::cout << std::endl;

					}
					if (traceRootCall || traceFocusedContCall) {
						traceBeforeRootMakeNode = regressionDiagnostics;
						for (std::size_t edgeIndex = 0; edgeIndex < e.size(); ++edgeIndex) {
							std::cerr << (traceFocusedContCall ? "cont_focused_edges" : "cont_root_edges") << "\tstep\t" << traceStep
								<< "\tedge\t" << edgeIndex
								<< "\tw_re\t" << CTEntry::val(e[edgeIndex].w.r)
								<< "\tw_im\t" << CTEntry::val(e[edgeIndex].w.i)
								<< "\tchild_var\t" << static_cast<int>(e[edgeIndex].p->v)
								<< "\tmap_level\t" << (e[edgeIndex].map ? static_cast<int>(e[edgeIndex].map->level) : -999)
								<< "\tmap_x\t" << (e[edgeIndex].map ? static_cast<int>(e[edgeIndex].map->x) : -1)
								<< "\tmap_rot\t" << (e[edgeIndex].map ? e[edgeIndex].map->rotate : -999)
								<< "\tmap_ep\t" << (e[edgeIndex].map ? e[edgeIndex].map->extra_phase : -999)
								<< "\n";
						}
					}
					r = makeDDNode(Qubit(newk1), e, true);
					if (traceRootCall) {
						std::cerr << "cont_root_after_make\tstep\t" << traceStep
							<< "\tbranch\teq"
							<< "\tr_w_re\t" << CTEntry::val(r.w.r)
							<< "\tr_w_im\t" << CTEntry::val(r.w.i)
							<< "\tr_var\t" << static_cast<int>(r.p->v)
							<< "\tr_map_level\t" << (r.map ? static_cast<int>(r.map->level) : -999)
							<< "\tr_map_x\t" << (r.map ? static_cast<int>(r.map->x) : -1)
							<< "\tr_map_rot\t" << (r.map ? r.map->rotate : -999)
							<< "\tr_map_ep\t" << (r.map ? r.map->extra_phase : -999)
							<< "\n";
						traceSawRootMakeNode = true;
						traceAfterRootMakeNode = regressionDiagnostics;
					}
				}
			}
			// if(test_1314){
				
			// 	if(r == xCopy ){
			// 		std::cout << "true" << std::endl;
			// 	}
				
			// 	else{
			// 		if(r.p == xCopy.p) std::cout << "p right" << std::endl;
			// 		else std::cout << "p wrong" << std::endl;

			// 		std::cout << "r w:" << r.w << " x w:" << xCopy.w << std::endl;
			// 		r.map->print_maps(r.map);
			// 		r.map->print_maps(xCopy.map);
			// 		std::cout << "r size: "<< size(r) << " x size: " << size(xCopy) << std::endl;
			// 	}
			// }
			const auto traceContInsertShape = [&]() {
				if (std::getenv("LIMTDD_CONT_INSERT_TRACE") == nullptr) {
					return false;
				}
				if (const auto* startEnv = std::getenv("LIMTDD_CONT_INSERT_TRACE_STEP_START")) {
					if (contStageTraceStep < static_cast<std::size_t>(std::stoull(startEnv))) {
						return false;
					}
				}
				if (const auto* endEnv = std::getenv("LIMTDD_CONT_INSERT_TRACE_STEP_END")) {
					if (contStageTraceStep > static_cast<std::size_t>(std::stoull(endEnv))) {
						return false;
					}
				}
				const auto rRe = CTEntry::val(r.w.r);
				const auto rIm = CTEntry::val(r.w.i);
				return var_num == 1 && std::abs(rRe - 2.0) < 1e-9 && std::abs(rIm) < 1e-9 &&
					yCopy.p != nullptr && yCopy.p->v == 1 && yCopy.map != nullptr && yCopy.map->level == -1;
			}();
			if (traceContInsertShape) {
				std::cerr << "cont_cache_insert_candidate"
					<< "\tstep\t" << contStageTraceStep
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tvar_num\t" << var_num
					<< "\tx_var\t" << static_cast<int>(xCopy.p->v)
					<< "\tx_map_level\t" << (xCopy.map ? static_cast<int>(xCopy.map->level) : -999)
					<< "\tx_map_x\t" << (xCopy.map ? static_cast<int>(xCopy.map->x) : -1)
					<< "\tx_map_rot\t" << (xCopy.map ? xCopy.map->rotate : -999)
					<< "\tx_map_ep\t" << (xCopy.map ? xCopy.map->extra_phase : -999)
					<< "\ty_var\t" << static_cast<int>(yCopy.p->v)
					<< "\ty_map_level\t" << (yCopy.map ? static_cast<int>(yCopy.map->level) : -999)
					<< "\ty_map_x\t" << (yCopy.map ? static_cast<int>(yCopy.map->x) : -1)
					<< "\ty_map_rot\t" << (yCopy.map ? yCopy.map->rotate : -999)
					<< "\ty_map_ep\t" << (yCopy.map ? yCopy.map->extra_phase : -999)
					<< "\tr_w_re\t" << CTEntry::val(r.w.r)
					<< "\tr_w_im\t" << CTEntry::val(r.w.i)
					<< "\tr_var\t" << static_cast<int>(r.p->v)
					<< "\tr_map_level\t" << (r.map ? static_cast<int>(r.map->level) : -999)
					<< "\tkey1_level\t" << (temp_key_2_new_key1 ? temp_key_2_new_key1->level : -999)
					<< "\tkey1_new\t" << (temp_key_2_new_key1 ? temp_key_2_new_key1->new_key : -999)
					<< "\tkey2_level\t" << (temp_key_2_new_key2 ? temp_key_2_new_key2->level : -999)
					<< "\tkey2_new\t" << (temp_key_2_new_key2 ? temp_key_2_new_key2->new_key : -999)
					<< "\n";
			}
			if (traceChildCall) {
				traceRootAggEdge("final", "prexy_r", -1, r);
				std::cerr << "cont_child_prexy\tstep\t" << traceStep
					<< "\tparent_branch\t" << traceChildBranch
					<< "\tparent_k\t" << traceChildK
					<< "\tdepth\t" << contStageTraceDepth
					<< "\tr_w_re\t" << CTEntry::val(r.w.r)
					<< "\tr_w_im\t" << CTEntry::val(r.w.i)
					<< "\tx_w_re\t" << CTEntry::val(x.w.r)
					<< "\tx_w_im\t" << CTEntry::val(x.w.i)
					<< "\ty_w_re\t" << CTEntry::val(y.w.r)
					<< "\ty_w_im\t" << CTEntry::val(y.w.i)
					<< "\n";
			}
			if (std::getenv("LIMTDD_DISABLE_CONT_CACHE") == nullptr) {
				contTable.insert(xCopy, yCopy, { r.p, r.w,r.map }, temp_key_2_new_key1, temp_key_2_new_key2, var_num);
			}

			if (traceRootCall) {
				traceRootAggEdge("final", "prexy_r", -1, r);
				std::cerr << "cont_root_prexy\tstep\t" << traceStep
					<< "\tr_w_re\t" << CTEntry::val(r.w.r)
					<< "\tr_w_im\t" << CTEntry::val(r.w.i)
					<< "\tx_w_re\t" << CTEntry::val(x.w.r)
					<< "\tx_w_im\t" << CTEntry::val(x.w.i)
					<< "\ty_w_re\t" << CTEntry::val(y.w.r)
					<< "\ty_w_im\t" << CTEntry::val(y.w.i)
					<< "\n";
			}
			if (!r.w.exactlyZero() && (x.w.exactlyOne() || !y.w.exactlyZero())) {
				if (r.w.exactlyOne()) {
					r.w = cn.mulCached(x.w, y.w);
				}
				else {
					assert(r.w != Complex::zero);
					ComplexNumbers::mul(r.w, r.w, x.w);
					ComplexNumbers::mul(r.w, r.w, y.w);
				}
				if (r.w.approximatelyZero()) {
					//assert(r.w != Complex::zero);
					// cn.returnToCache(r.w);
					// cn.returnToCache(extra_phase);
					return ResultEdge::zero;
				}
			}
			if (r.w == Complex::zero) {
				// cn.returnToCache(extra_phase);
				return ResultEdge::zero;
			}
			else {
				if (traceChildCall) {
					std::cerr << "cont_child_build_prephase\tstep\t" << traceStep
						<< "\tparent_branch\t" << traceChildBranch
						<< "\tparent_k\t" << traceChildK
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tr_w_re\t" << CTEntry::val(r.w.r)
						<< "\tr_w_im\t" << CTEntry::val(r.w.i)
						<< "\n";
				}
				if (traceRootCall) {
					std::cerr << "cont_root_build_prephase\tstep\t" << traceStep
						<< "\tr_w_re\t" << CTEntry::val(r.w.r)
						<< "\tr_w_im\t" << CTEntry::val(r.w.i)
						<< "\n";
				}
				r.map = mapmul(r_maps->remain_map, r.map);
				assert(r.w != Complex::zero);
				// cn.mul(r.w, r.w, r.map->extra_phase);
				cn.mul(r.w, r.w, cn.getTemporary(cos(r.map->extra_phase*rotate_angle),sin(r.map->extra_phase*rotate_angle)));
				// cn.returnToCache(r.map->extra_phase);
				assert(r.w != Complex::zero);
				// cn.mul(r.w, r.w, extra_phase);
				cn.mul(r.w, r.w, cn.getTemporary(cos(extra_phase*rotate_angle),sin(extra_phase*rotate_angle)));
				// cn.returnToCache(extra_phase);
				if (traceChildCall) {
					std::cerr << "cont_child_build_postphase\tstep\t" << traceStep
						<< "\tparent_branch\t" << traceChildBranch
						<< "\tparent_k\t" << traceChildK
						<< "\tdepth\t" << contStageTraceDepth
						<< "\tr_w_re\t" << CTEntry::val(r.w.r)
						<< "\tr_w_im\t" << CTEntry::val(r.w.i)
						<< "\n";
				}
				if (traceRootCall) {
					std::cerr << "cont_root_build_postphase\tstep\t" << traceStep
						<< "\tr_w_re\t" << CTEntry::val(r.w.r)
						<< "\tr_w_im\t" << CTEntry::val(r.w.i)
						<< "\n";
				}
			}
			const auto traceAfterFinalize = traceRootCall ? regressionDiagnostics : RegressionDiagnostics{};
			if (traceRootCall) {
				if (traceSawRootMakeNode) {
					printContCounterDelta("recursive", traceAfterFindRemain, traceBeforeRootMakeNode);
					printContCounterDelta("tadd_total", RegressionDiagnostics{}, contStageTaddTotals);
					printContCounterDelta("make_node_total", RegressionDiagnostics{}, contStageMakeNodeTotals);
					printContCounterDelta("root_make_node", traceBeforeRootMakeNode, traceAfterRootMakeNode);
					printContStageDelta("finalize", traceAfterRootMakeNode, traceAfterFinalize, traceEdgeNodesBefore, size(r));
				} else {
					printContCounterDelta("tadd_total", RegressionDiagnostics{}, contStageTaddTotals);
					printContCounterDelta("make_node_total", RegressionDiagnostics{}, contStageMakeNodeTotals);
					printContStageDelta("build", traceAfterFindRemain, traceAfterFinalize, traceEdgeNodesBefore, size(r));
				}
			}

			
			//std::cout << "Case 2 " << r.w << " " << int(r.p->v) << std::endl;
			//std::cout << "<<<<<< cont2 >>>>>>>" << std::endl; 
			return r;
		}

		//==========================================我写的========================================


	public:
		///
		/// Decision diagram size
		///
		template <class Edge> unsigned int size(const Edge& e) {
			static constexpr unsigned int NODECOUNT_BUCKETS = 200000;
			static std::unordered_set<decltype(e.p)> visited{NODECOUNT_BUCKETS}; // 2e6
			visited.max_load_factor(10);
			visited.clear();
			return nodeCount(e, visited);
		}

	private:
		template <class Edge>
		unsigned int nodeCount(const Edge& e,
			std::unordered_set<decltype(e.p)>& v) const {
			v.insert(e.p);
			unsigned int sum = 1;
			if (!e.isTerminal()) {
				for (const auto& edge : e.p->e) {
					if (edge.p != nullptr && !v.count(edge.p)) {
						sum += nodeCount(edge, v);
					}
				}
			}
			return sum;
		}


		///
		/// Printing and Statistics
		///
	public:
		// print information on package and its members
		static void printInformation() {
			std::cout << "\n  compiled: " << __DATE__ << " " << __TIME__
				<< "\n  Complex size: " << sizeof(Complex) << " bytes (aligned "
				<< alignof(Complex) << " bytes)"
				<< "\n  ComplexValue size: " << sizeof(ComplexValue)
				<< " bytes (aligned " << alignof(ComplexValue) << " bytes)"
				<< "\n  ComplexNumbers size: " << sizeof(ComplexNumbers)
				<< " bytes (aligned " << alignof(ComplexNumbers) << " bytes)"
				<< "\n  mEdge size: " << sizeof(mEdge) << " bytes (aligned "
				<< alignof(mEdge) << " bytes)"
				<< "\n  mNode size: " << sizeof(mNode) << " bytes (aligned "
				<< alignof(mNode) << " bytes)"
				<< "\n  Package size: " << sizeof(Package) << " bytes (aligned "
				<< alignof(Package) << " bytes)"
				<< "\n"
				<< std::flush;
		}

		// print unique and compute table statistics

		void statistics() {
			std::cout << "DD statistics:\n";
			std::cout << "[UniqueTable] ";
			nodeUniqueTable.printStatistics();
			std::cout << "[Add] ";
			addTable.printStatistics();
			std::cout << "[Cont] ";
			contTable.printStatistics();
			std::cout << "[ComplexTable] ";
			cn.complexTable.printStatistics();
		}

	};


} // namespace dd
