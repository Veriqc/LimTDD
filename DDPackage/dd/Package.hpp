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
#include "Tensor.hpp"
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

namespace dd {

	template <class Config> class Package {
		static_assert(std::is_base_of_v<DDPackageConfig, Config>, "Config must be derived from DDPackageConfig");

		///
		/// Complex number handling
		///
	public:
		ComplexNumbers cn{};

		///
		/// Construction, destruction, information and reset
		///

		static constexpr std::size_t MAX_POSSIBLE_QUBITS = static_cast<std::make_unsigned_t<Qubit>>(std::numeric_limits<Qubit>::max()) + 1U;
		static constexpr std::size_t DEFAULT_QUBITS = 300;


		//==========================================我写的========================================
		bool to_test = false;

		std::map<std::string, int> varOrder;

		//==========================================我写的========================================

		explicit Package(std::size_t nq = DEFAULT_QUBITS) : nqubits(nq) {
			resize(nq);
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

	private:
		std::size_t nqubits;

		///
		/// Vector nodes, edges and quantum states
		///
	public:


		//==========================================我写的========================================
		//template <class Node>

		//TDD Tensor_2_TDD(const Tensor tn)
		//{

		//	TDD res;

		//	if()




		//	return res;
		//}

		TDD Matrix2TDD(const GateMatrix mat, std::vector<Index> var_out)
		{

			TDD low, high, res;
			Edge<mNode> e_temp[4];

			std::vector<Edge<mNode>> e_low(2), e_high(2), e(2);

			int Radix = 2;

			for (int i = 0; i < Radix; i++) {
				for (int j = 0; j < Radix; j++) {
					if (mat[2 * i + j] == complex_zero) {
						e_temp[i * Radix + j] = Edge<mNode>::zero;
					}
					else if (mat[2 * i + j] == complex_one) {
						e_temp[i * Radix + j] = Edge<mNode>::one;
					}
					else {
						e_temp[i * Radix + j] = Edge<mNode>::terminal(cn.lookup(mat[2 * i + j]));
					}
				}
			}


			std::vector<std::string> key_2_index;
			if (varOrder[var_out[0].key] < varOrder[var_out[1].key]) {
				e_low[0] = e_temp[0];
				e_low[1] = e_temp[1];
				e_high[0] = e_temp[2];
				e_high[1] = e_temp[3];
				key_2_index = { var_out[0].key,var_out[1].key };
			}
			else {
				e_low[0] = e_temp[0];
				e_low[1] = e_temp[2];
				e_high[0] = e_temp[1];
				e_high[1] = e_temp[3];
				key_2_index = { var_out[1].key,var_out[0].key };
			}


			if (e_low[0].p == e_low[1].p and e_low[0].w == e_low[1].w) {
				low.e = e_low[0];
			}
			else {
				low.e = makeDDNode(0, e_low, false);
			}
			if (e_high[0].p == e_high[1].p and e_high[0].w == e_high[1].w) {
				high.e = e_high[0];
			}
			else {
				high.e = makeDDNode(0, e_high, false);
			}
			if (low.e.p == high.e.p and low.e.w == high.e.w) {
				res.e = low.e;
			}
			else {
				e[0] = low.e;
				e[1] = high.e;
				res.e = makeDDNode(1, e, false); ;
			}
			res.index_set = var_out;
			res.key_2_index = key_2_index;
			return res;
		}
		//template <class Node>
		TDD diag_matrix_2_TDD(const GateMatrix mat, std::vector<Index> var_out)
		{

			TDD res;
			int Radix = 2;

			std::vector<Edge<mNode>> e_temp(2);
			for (int i = 0; i < Radix; i++) {

				if (mat[2 * i + i] == complex_zero) {
					e_temp[i] = Edge<mNode>::zero;
				}
				else {
					if (mat[2 * i + i] == complex_one) {
						e_temp[i] = Edge<mNode>::one;
					}
					else {
						e_temp[i] = Edge<mNode>::terminal(cn.lookup(mat[2 * i + i]));
					}
				}
			}

			res.e = makeDDNode(0, e_temp, false);
			res.index_set = var_out;
			res.key_2_index = { var_out[0].key };
			return res;
		}

		//template <class Node>
		TDD cnot_2_TDD(std::vector<Index> var, int ca = 1) {


			TDD low, high, res;
			std::vector<Edge<mNode>> e(2);
			if (ca == 1) {
				if (varOrder[var[0].key] > varOrder[var[3].key] && varOrder[var[0].key] > varOrder[var[4].key]) {
					low = Matrix2TDD(Imat, { var[3] ,var[4] });
					high = Matrix2TDD(Xmat, { var[3] ,var[4] });
					e[0] = low.e;
					e[1] = high.e;
					res.e = makeDDNode(2, e, false);
					res.index_set = { var[0],var[2],var[3],var[4] };
					low.key_2_index.push_back(var[0].key);
					res.key_2_index = low.key_2_index;
				}
				else if (varOrder[var[3].key] > varOrder[var[0].key] && varOrder[var[3].key] > varOrder[var[4].key]) {
					low = Matrix2TDD(Imat, { var[0] ,var[4] });
					high = Matrix2TDD(Xmat, { var[0] ,var[4] });
					e[0] = low.e;
					e[1] = high.e;
					res.e = makeDDNode(2, e, false);
					res.index_set = { var[0],var[2],var[3],var[4] };
					low.key_2_index.push_back(var[3].key);
					res.key_2_index = low.key_2_index;
				}
				else {
					low = Matrix2TDD(Imat, { var[0] ,var[3] });
					high = Matrix2TDD(Xmat, { var[0] ,var[3] });
					e[0] = low.e;
					e[1] = high.e;
					res.e = makeDDNode(2, e, false);
					res.index_set = { var[0],var[2],var[3],var[4] };
					low.key_2_index.push_back(var[4].key);
					res.key_2_index = low.key_2_index;
				}
				return res;
			}
			return res;

		}

		//==========================================我写的========================================

		///
		/// Matrix nodes, edges and quantum gates
		///
		template <class Node> Edge<Node> normalize(const Edge<Node>& e, bool cached) {

			auto argmax = -1;

			//auto zero = std::array{	e.p->e[0].w.approximatelyZero(), e.p->e[1].w.approximatelyZero()};
			int R = e.p->e.size();
			std::vector<bool> zero;
			for (int k = 0; k < R; k++) {
				zero.push_back(e.p->e[k].w.approximatelyZero());
			}

			// make sure to release cached numbers approximately zero, but not exactly
			// zero
			if (cached) {
				for (auto i = 0U; i < R; i++) {
					if (zero[i] && e.p->e[i].w != Complex::zero) {
						//assert(e.p->e[i].w != Complex::zero);
						cn.returnToCache(e.p->e[i].w);
						e.p->e[i] = Edge<Node>::zero;
					}
				}
			}

			fp max = 0;
			auto maxc = Complex::one;
			// determine max amplitude
			for (auto i = 0U; i < R; ++i) {
				if (zero[i]) {
					continue;
				}
				if (argmax == -1) {
					argmax = static_cast<decltype(argmax)>(i);
					max = ComplexNumbers::mag2(e.p->e[i].w);
					maxc = e.p->e[i].w;
				}
				else {
					auto mag = ComplexNumbers::mag2(e.p->e[i].w);
					if (mag - max > ComplexTable<>::tolerance()) {
						argmax = static_cast<decltype(argmax)>(i);
						max = mag;
						maxc = e.p->e[i].w;
					}
				}
			}

			// all equal to zero
			if (argmax == -1) {
				if (!cached && !e.isTerminal()) {
					// If it is not a cached computation, the node has to be put back into
					// the chain
					getUniqueTable<Node>().returnNode(e.p);
				}
				return Edge<Node>::zero;
			}

			short x = 0;

			auto r=e;
			if (argmax > 0) {
				r.p->e = {r.p->e[1],r.p->e[0]};
				zero = { zero[1],zero[0] };
				x = 1;
			}
			argmax = 0;
			// divide each entry by max
			for (auto i = 0U; i < R; ++i) {
				if (static_cast<decltype(argmax)>(i) == argmax) {
					if (cached) {
						if (r.w.exactlyOne()) {
							r.w = maxc;
						}
						else {
							ComplexNumbers::mul(r.w, r.w, maxc);
						}
					}
					else {
						if (r.w.exactlyOne()) {
							r.w = maxc;
						}
						else {
							auto c = cn.getTemporary();
							ComplexNumbers::mul(c, r.w, maxc);
							r.w = cn.lookup(c);
						}
					}
					r.map = r.p->e[i].map;
					r.p->e[i].w = Complex::one;
					r.p->e[i].map = the_maps::the_maps_header();
				}
				else {
					if (zero[i]) {
						if (cached && r.p->e[i].w != Complex::zero) {
							//assert(r.p->e[i].w != Complex::zero);
							cn.returnToCache(r.p->e[i].w);
						}
						//r.p->e[i] = Edge<Node>::zero;
						r.p->e[i] = { r.p->e[0].p,Complex::zero, the_maps::the_maps_header() };
						continue;
					}
					if (cached && !zero[i] && !r.p->e[i].w.exactlyOne()) {
						//assert(r.p->e[i].w != Complex::zero);
						cn.returnToCache(r.p->e[i].w);
					}
					if (r.p->e[i].w.approximatelyOne()) {
						r.p->e[i].w = Complex::one;
					}
					auto c = cn.getTemporary();
					ComplexNumbers::div(c, r.p->e[i].w, maxc);
					float angle = atan(c.i->value / c.r->value);
					c.r->value = ComplexNumbers::mag2(c);
					c.i->value = 0;
					r.p->e[i].w = cn.lookup(c);
					r.p->e[i].map = the_maps::mapdiv(r.p->e[i].map, r.p->e[0].map);
					r.p->e[i].map->extra_phase += (long int) angle/ unit_rotate_angle;
				}
			}

			r.map = the_maps::append_new_map(r.map, r.p->v, x, (r.p->e[1].map->extra_phase)% root_of_unit);

			return r;

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

			auto& uniqueTable = getUniqueTable<Node>();
			Edge<Node> e{uniqueTable.getNode(), Complex::one};
			e.p->v = var;
			e.p->e = edges;

			assert(e.p->ref == 0);


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

			e = normalize(e, cached);

			assert(e.p->v == var || e.isTerminal());

			// look it up in the unique tables
			auto l = uniqueTable.lookup(e, false);

			assert(l.p->v == var || l.isTerminal());

			return l;
		}


		///
		/// Compute table definitions
		///
	public:
		void clearComputeTables() {

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
				self->next[new_key] = new key_2_new_key_node{ self->level + 1, new_key, {}, self };

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


			if (to_test) {
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
						last_cont_idx = last_cont_idx + 1 / (3 * nqubits) * repeat_time;
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
				cn.returnToCache(res.e.w);
				res.e.w = cn.lookup(res.e.w);
			}

			[[maybe_unused]] const auto after = cn.cacheCount();
			assert(before == after);

			the_maps::print_maps(res.e.map);

			return res;
		}


	private:

		template <class Node>
		Edge<Node> Slicing(const Edge<Node>& e, int x, int c) {
			// used for add
			assert(e.w != Complex::zero);
			if (e.p->v == -1) {
				return e;
			}
			if (e.p->v < x) {
				return e;
			}
			if (e.p->v == x) {
				if (e.p->v != e.map->level) {
					auto temp = e.p->e[c];
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = the_maps::mapmul(e.map, temp.map);
						if (temp.map->extra_phase > 0) {
							double angle = temp.map->extra_phase * unit_rotate_angle;
							cn.mul(temp.w, temp.w, cn.getTemporary(cos(angle), sin(angle)));
						}
					}
					return temp;
				}
				else if (e.map->x == 0) {
					auto temp = e.p->e[c];
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = the_maps::mapmul(e.map, temp.map);
						if (c == 0) {
							if (temp.map->extra_phase > 0) {
								double angle = temp.map->extra_phase * unit_rotate_angle;
								cn.mul(temp.w, temp.w, cn.getTemporary(cos(angle), sin(angle)));
							}
						}
						else {
							if (temp.map->extra_phase + e.map->rotate > 0) {
								double angle = ((temp.map->extra_phase + e.map->rotate) % root_of_unit) * unit_rotate_angle;
								cn.mul(temp.w, temp.w, cn.getTemporary(cos(angle), sin(angle)));
							}
						}
					}
					return temp;
				}
				else {
					auto temp = e.p->e[1 - c];
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = the_maps::mapmul(e.map, temp.map);
						if (c == 1) {
							if (temp.map->extra_phase > 0) {
								double angle = temp.map->extra_phase * unit_rotate_angle;
								cn.mul(temp.w, temp.w, cn.getTemporary(cos(angle), sin(angle)));
							}
						}
						else {
							if (temp.map->extra_phase + e.map->rotate > 0) {
								double angle = ((temp.map->extra_phase + e.map->rotate) % root_of_unit) * unit_rotate_angle;
								cn.mul(temp.w, temp.w, cn.getTemporary(cos(angle), sin(angle)));
							}
						}
					}
					return temp;
				}

			}
			else {
				std::cout << "Slicing not support yet" << std::endl;
				return e;
			}

		}

		template <class Node>
		Edge<Node> Slicing2(const Edge<Node>& e, int x, int c) {

			assert(e.w != Complex::zero);
		// used for contract
			if (e.p->v == -1) {
				return e;
			}
			if (e.p->v < x) {
				return e;
			}
			if (e.p->v == x) {
				if (e.p->v != e.map->level) {
					auto temp = e.p->e[c];
					if (temp.w != Complex::zero) {
						temp.map = the_maps::mapmul(e.map, temp.map);
						double angle = temp.map->extra_phase * unit_rotate_angle;
						temp.w=cn.mulCached(temp.w, cn.getTemporary(cos(angle), sin(angle)));
					}
					return temp;
				}
				else if (e.map->x == 0) {
					auto temp = e.p->e[c];
					if (temp.w != Complex::zero) {
						temp.map = the_maps::mapmul(e.map, temp.map);
						double angle;
						if (c == 0) {
							angle = temp.map->extra_phase * unit_rotate_angle;
						}
						else {
							angle = ((temp.map->extra_phase + e.map->rotate) % root_of_unit) * unit_rotate_angle;
						}
						temp.w = cn.mulCached(temp.w, cn.getTemporary(cos(angle), sin(angle)));
					}
					return temp;
				}
				else {
					auto temp = e.p->e[1-c];
					if (temp.w != Complex::zero) {
						temp.map = the_maps::mapmul(e.map, temp.map);
						double angle;
						if (c == 1) {
							angle = temp.map->extra_phase * unit_rotate_angle;
						}
						else {
							angle = ((temp.map->extra_phase + e.map->rotate) % root_of_unit) * unit_rotate_angle;
						}
						temp.w = cn.mulCached(temp.w, cn.getTemporary(cos(angle), sin(angle)));
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

			if (x.p > y.p) {
				return T_add2(y, x);
			}

			if (x.p == nullptr) {
				return y;
			}
			if (y.p == nullptr) {
				return x;
			}

			if (x.w.exactlyZero()) {
				if (y.w.exactlyZero()) {
					return Edge<Node>::zero;
				}
				auto r = y;
				r.w = cn.getCached(CTEntry::val(y.w.r), CTEntry::val(y.w.i));
				return r;
			}
			if (y.w.exactlyZero()) {
				auto r = x;
				r.w = cn.getCached(CTEntry::val(x.w.r), CTEntry::val(x.w.i));
				return r;
			}
			if (x.p == y.p && x.map==y.map) {
				auto r = y;
				r.w = cn.addCached(x.w, y.w);
				if (r.w.approximatelyZero()) {
					//assert(r.w != Complex::zero);
					cn.returnToCache(r.w);
					return Edge<Node>::zero;
				}
				r.map = x.map;
				return r;
			}

			auto xCopy = x;
			auto yCopy = y;
			if (x.w != Complex::one) {
				xCopy.w = Complex::one;
				xCopy.map = the_maps::the_maps_header();
				yCopy.w = cn.divCached(y.w, x.w);
				yCopy.map = the_maps::mapdiv(y.map, x.map);
				if (yCopy.map->extra_phase > 0 && yCopy.w != Complex::zero) {
					double angle = yCopy.map->extra_phase * unit_rotate_angle;
					cn.mul(yCopy.w, yCopy.w, cn.getTemporary(cos(angle), sin(angle)));
				}
			}


			auto r = addTable.lookup({ xCopy.p, xCopy.w,xCopy.map }, { yCopy.p, yCopy.w,yCopy.map });

			if (r.p != nullptr) {
				if (r.w.approximatelyZero()) {
					return Edge<Node>::zero;
				}
				auto c = cn.getCached(r.w);
				if (x.w != Complex::one) {
					cn.mul(c, c, x.w);
					//assert(yCopy.w != Complex::zero);
					cn.returnToCache(yCopy.w);
				}
				auto temp_map = the_maps::mapmul(x.map, r.map);
				if (temp_map->extra_phase > 0 && c != Complex::zero) {
					double angle = temp_map->extra_phase * unit_rotate_angle;
					cn.mul(c, c, cn.getTemporary(cos(angle), sin(angle)));
				}

				return { r.p, c,temp_map };
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
					e1 = Slicing(x, x.p->v, i);
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
					e2 = Slicing(y, y.p->v, i);
				}
				else {
					e2 = yCopy;
					if (x.p->e[i].p == nullptr) {
						e2 = { nullptr, Complex::zero };
					}
				}


				edge[i] = T_add2(e1, e2);


				if (!x.isTerminal() && x.p->v == w && e1.w != Complex::zero) {
					//assert(e1.w != Complex::zero);
					cn.returnToCache(e1.w);
				}

				if (!y.isTerminal() && y.p->v == w && e2.w != Complex::zero) {
					//assert(e2.w != Complex::zero);
					cn.returnToCache(e2.w);
				}
			}

			auto e = makeDDNode(w, edge, true);

			addTable.insert({ xCopy.p,xCopy.w,xCopy.map }, { yCopy.p,yCopy.w,yCopy.map }, { e.p, e.w,e.map });
			if (x.w != Complex::one) {
				cn.mul(e.w, e.w, x.w);
				//assert(yCopy.w != Complex::zero);
				cn.returnToCache(yCopy.w);
			}

			e.map = the_maps::mapmul(x.map, e.map);
			if (e.map->extra_phase > 0 && e.w != Complex::zero) {
				double angle = e.map->extra_phase * unit_rotate_angle;
				cn.mul(e.w, e.w, cn.getTemporary(cos(angle), sin(angle)));
			}
			return e;
		}


		comm_maps* find_remain_map(the_maps* map1, the_maps* map2, key_2_new_key_node* key_2_new_key1, key_2_new_key_node* key_2_new_key2) {

			//the_maps* res[3];
			//std::cout << 868 << "   " << map1->level << " " << map2->level << std::endl;

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

			if (newk1 > newk2 && int(newk1 * 2) % 2 != 1) {
				auto res = find_remain_map(map1->father, map2, temp_key_2_new_key1, temp_key_2_new_key2);
				auto temp_pahse = res->remain_map->extra_phase;
				res->remain_map = the_maps::append_new_map(res->remain_map, newk1, map1->x, map1->rotate);
				res->remain_map->extra_phase = temp_pahse;
				return res;
			}
			if (newk1 < newk2 && int(newk2 * 2) % 2 != 1) {
				auto res = find_remain_map(map1, map2->father, temp_key_2_new_key1, temp_key_2_new_key2);
				auto temp_pahse = res->remain_map->extra_phase;
				res->remain_map = the_maps::append_new_map(res->remain_map, newk2, map2->x, map2->rotate);
				res->remain_map->extra_phase = temp_pahse;
				return res;
			}
			if (map1->level == -1 && map2->level == -1) {
				comm_maps* res=new comm_maps{ the_maps::the_maps_header(),the_maps::the_maps_header(),the_maps::the_maps_header() };
				res->remain_map->extra_phase = 0;
				return res;
			}
			if (newk1 > newk2) {
				auto res = find_remain_map(map1->father, map2, temp_key_2_new_key1, temp_key_2_new_key2);
				res->cont_map1 = the_maps::append_new_map(res->cont_map1, map1->level, map1->x, map1->rotate);
				return res;
			}
			if (newk1 < newk2) {
				auto res = find_remain_map(map1, map2->father, temp_key_2_new_key1, temp_key_2_new_key2);
				res->cont_map2 = the_maps::append_new_map(res->cont_map2, map2->level, map2->x, map2->rotate);
				return res;
			}
			auto res = find_remain_map(map1->father, map2->father, temp_key_2_new_key1, temp_key_2_new_key2);

			auto x = (map1->x + map2->x) % 2;
			res->remain_map->extra_phase += map2->rotate * x;
			auto rotate = (long int)(map1->rotate + map2->rotate * pow(-1, x));
			rotate = rotate % root_of_unit;
			res->cont_map1 = the_maps::append_new_map(res->cont_map1, map1->level, x, rotate);

			return res;
		}

		//template <class LeftOperandNode, class RightOperandNode>
		Edge<mNode> cont2(const Edge<mNode>& x, const Edge<mNode>& y, key_2_new_key_node* key_2_new_key1, key_2_new_key_node* key_2_new_key2, const int var_num) {

			//std::cout <<"838 " << x.w << " " << y.w << " " << int(x.p->v) << " " << int(y.p->v) << std::endl;

			using ResultEdge = Edge<mNode>;

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

				if (var_num > 0) {
					ComplexNumbers::mul(c, c, cn.getTemporary(pow(2, var_num), 0));
				}

				return ResultEdge::terminal(c);
			}

			key_2_new_key_node* temp_key_2_new_key2 = key_2_new_key2;
			while (temp_key_2_new_key2->level > y.p->v) {
				temp_key_2_new_key2 = temp_key_2_new_key2->father;
			}


			if (x.p->v == -1 && var_num == 0 && std::abs(temp_key_2_new_key2->new_key - y.p->v) < 1e-10) {
				return 	ResultEdge{ y.p, cn.mulCached(x.w, y.w) ,y.map};
			}

			key_2_new_key_node* temp_key_2_new_key1 = key_2_new_key1;
			while (temp_key_2_new_key1->level > x.p->v) {
				temp_key_2_new_key1 = temp_key_2_new_key1->father;
			}

			if (y.p->v == -1 && var_num == 0 && std::abs(temp_key_2_new_key1->new_key - x.p->v) < 1e-10) {
				return 	ResultEdge{ x.p, cn.mulCached(x.w, y.w) ,x.map};
			}


			auto xCopy = x;
			xCopy.w = Complex::one;
			auto yCopy = y;
			yCopy.w = Complex::one;

			
			auto r_maps = find_remain_map(x.map, y.map, key_2_new_key1, key_2_new_key2);
			xCopy.map = r_maps->cont_map1;
			yCopy.map = r_maps->cont_map2;


			auto res = contTable.lookup(xCopy, yCopy, temp_key_2_new_key1, temp_key_2_new_key2);
			if (res.e.p != nullptr) {
				if (res.e.w.approximatelyZero()) {
					return ResultEdge::zero;
				}
				auto e = ResultEdge{ res.e.p, cn.getCached(res.e.w) };
				ComplexNumbers::mul(e.w, e.w, x.w);
				ComplexNumbers::mul(e.w, e.w, y.w);
				if (e.w.approximatelyZero()) {
					//assert(e.w != Complex::zero);
					cn.returnToCache(e.w);
					return ResultEdge::zero;
				}
				if (res.cont_num != var_num) {
					ComplexNumbers::mul(e.w, e.w, cn.getTemporary(pow(2, var_num - res.cont_num), 0));//对于一般形状的tensor,以2为底数可能有问题
				}
				e.map = the_maps::mapmul(r_maps->remain_map, e.map);
				auto temp_phase = e.map->extra_phase + r_maps->remain_map->extra_phase;
				temp_phase = temp_phase % root_of_unit;
				if (temp_phase > 0 && e.w != Complex::zero) {
					double angle = temp_phase * unit_rotate_angle;
					cn.mul(e.w, e.w, cn.getTemporary(cos(angle), sin(angle)));
				}
				return e;
			}


			float newk1 = temp_key_2_new_key1->new_key;

			float newk2 = temp_key_2_new_key2->new_key;

			ResultEdge e1{}, e2{}, r{};

			if (newk1 > newk2) {
				if (int(newk1 * 2) % 2 == 1) {
					r = ResultEdge::zero;
					ResultEdge etemp{};
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						e1 = Slicing(x, x.p->v, k);
						e2 = yCopy;
						etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
						if (e1.w != Complex::zero) {
							cn.returnToCache(e1.w);
						}
						if (etemp.w != Complex::zero) {
							if (r != ResultEdge::zero) {
								auto temp = r.w;
								r = T_add2(r, etemp);
								//assert(temp != Complex::zero);
								//assert(etemp.w != Complex::zero);
								cn.returnToCache(temp);
								cn.returnToCache(etemp.w);
							}
							else {
								r = etemp;
							}
						}
					}
				}
				else {
					std::vector<ResultEdge> e;
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						e1 = Slicing(x, x.p->v, k);
						e2 = yCopy;
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						if (e1.w != Complex::zero) {
							cn.returnToCache(e1.w);
						}
					}
					r = makeDDNode(Qubit(newk1), e, true);
				}
			}
			else if (newk1 < newk2) {
				if (int(newk2 * 2) % 2 == 1) {
					r = ResultEdge::zero;
					ResultEdge etemp{};
					for (int k = 0; k < y.p->e.size(); ++k) {
						e1 = xCopy;
						//e2 = y.p->e[k];
						e2 = Slicing(y, y.p->v, k);
						etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
						if (e2.w != Complex::zero) {
							cn.returnToCache(e2.w);
						}
						if (etemp.w != Complex::zero) {
							if (r != ResultEdge::zero) {
								auto temp = r.w;
								r = T_add2(r, etemp);
								//assert(temp != Complex::zero);
								//assert(etemp.w != Complex::zero);
								cn.returnToCache(temp);
								cn.returnToCache(etemp.w);
							}
							else {
								r = etemp;
							}
						}
					}
				}
				else {
					std::vector<ResultEdge> e;
					for (int k = 0; k < y.p->e.size(); ++k) {
						e1 = xCopy;
						//e2 = y.p->e[k];
						e2 = Slicing(y, y.p->v, k);
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						if (e2.w != Complex::zero) {
							cn.returnToCache(e2.w);
						}
					}
					r = makeDDNode(Qubit(newk2), e, true);
				}

			}
			else {
				if (int(newk2 * 2) % 2 == 1) {
					r = ResultEdge::zero;
					ResultEdge etemp{};
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						e1 = Slicing(x, x.p->v, k);
						//e2 = y.p->e[k];
						e2 = Slicing(y, y.p->v, k);
						etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
						if (e1.w != Complex::zero) {
							cn.returnToCache(e1.w);
						}
						if (e2.w != Complex::zero) {
							cn.returnToCache(e2.w);
						}
						if (etemp.w != Complex::zero) {
							if (r != ResultEdge::zero) {
								auto temp = r.w;
								r = T_add2(r, etemp);
								//assert(temp != Complex::zero);
								//assert(etemp.w != Complex::zero);
								cn.returnToCache(temp);
								cn.returnToCache(etemp.w);
							}
							else {
								r = etemp;
							}
						}
					}
				}
				else {
					std::vector<ResultEdge> e;
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						e1 = Slicing(x, x.p->v, k);
						//e2 = y.p->e[k];
						e2 = Slicing(y, y.p->v, k);
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						if (e1.w != Complex::zero) {
							cn.returnToCache(e1.w);
						}
						if (e2.w != Complex::zero) {
							cn.returnToCache(e2.w);
						}
					}
					r = makeDDNode(Qubit(newk1), e, true);
				}
			}

			contTable.insert(xCopy, yCopy, { r.p, r.w,r.map }, temp_key_2_new_key1, temp_key_2_new_key2, var_num);

			if (!r.w.exactlyZero() && (x.w.exactlyOne() || !y.w.exactlyZero())) {
				if (r.w.exactlyOne()) {
					r.w = cn.mulCached(x.w, y.w);
				}
				else {
					ComplexNumbers::mul(r.w, r.w, x.w);
					ComplexNumbers::mul(r.w, r.w, y.w);
				}
				if (r.w.approximatelyZero()) {
					//assert(r.w != Complex::zero);
					cn.returnToCache(r.w);
					return ResultEdge::zero;
				}
			}
			r.map = the_maps::mapmul(r_maps->remain_map, r.map);
			auto temp_phase = r.map->extra_phase + r_maps->remain_map->extra_phase;
			temp_phase = temp_phase % root_of_unit;
			if (temp_phase > 0 && r.w != Complex::zero) {
				double angle = temp_phase * unit_rotate_angle;
				cn.mul(r.w, r.w, cn.getTemporary(cos(angle), sin(angle)));
			}
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
