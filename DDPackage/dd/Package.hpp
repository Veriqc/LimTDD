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

#include <xtensor/xarray.hpp>
#include <xtensor/xshape.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xslice.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xview.hpp>

namespace dd {

	bool mapCompare(the_maps* map1, the_maps* map2) {
		const auto maxArg = dd::PI;

		while (map1->level >= 0 || map2->level >= 0) {
			if (map1->level > map2->level) {
				if (-std::pow(-1, map1->x)*ComplexNumbers::arg(map1->rotate) > maxArg) return true;
				map1 = (map1->level > 0) ? map1->father : map1;
			} else if (map2->level > map1->level) {
				if (ComplexNumbers::arg(map2->rotate) > maxArg) return true;
				map2 = (map2->level > 0) ? map2->father : map2;
			} else { // map1->level == map2->level
				fp theta1 = ComplexNumbers::arg(map1->rotate);
				fp theta2 = ComplexNumbers::arg(map2->rotate);
				fp phaseDiff = theta2 - theta1 * std::pow(-1, map1->x ^ map2->x);
				if (phaseDiff > maxArg) return true;
				map1 = (map1->level > 0) ? map1->father : map1;
				map2 = (map2->level > 0) ? map2->father : map2;
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

		///
		/// Construction, destruction, information and reset
		///

		static constexpr std::size_t MAX_POSSIBLE_QUBITS = static_cast<std::make_unsigned_t<Qubit>>(std::numeric_limits<Qubit>::max()) + 1U;
		static constexpr std::size_t DEFAULT_QUBITS = 300;


		//==========================================我写的========================================
		bool to_test = false;

		int mode = 1;

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
						cn.returnToCache(e.p->e[i].w);
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
					max_value = e.p->e[i].w;
				}
				else {
					auto mag = ComplexNumbers::mag2(e.p->e[i].w);
					if (mag - max_mag2 > ComplexTable<>::tolerance()) {
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
			}
			else if(std::abs(ComplexNumbers::mag2(res.p->e[0].w)-ComplexNumbers::mag2(res.p->e[1].w)) < ComplexTable<>::tolerance()){
				add_x = false;
			}
			else if(res.p->e[0].p >  res.p->e[1].p){
				add_x = true ;
			}
			else if(res.p->e[0].p <  res.p->e[1].p){
				add_x = false;
			}
			else{
				add_x = mapCompare(res.p->e[0].map,res.p->e[1].map);
			}
				// if( && ComplexNumbers::arg(res.p->e[0].w) >  ComplexNumbers::arg(res.p->e[1].w)){
				// 	res.p->e = {res.p->e[1],res.p->e[0]};
				// 	isZero = { isZero[1],isZero[0] };
				// 	add_x = 1;
				// }

			if(add_x){
				res.p->e = {res.p->e[1],res.p->e[0]};
				isZero = { isZero[1],isZero[0] };
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
							cn.returnToCache(res.p->e[i].w);
						}
						//r.p->e[i] = Edge<Node>::zero;
						res.p->e[i] = { res.p->e[0].p,Complex::zero, the_maps::the_maps_header() };
						res.p->e[i].map->extra_phase = cn.getTemporary(1,0);
						continue;
					}
					if (cached && !isZero[i] && !res.p->e[i].w.exactlyOne()) {
						//assert(r.p->e[i].w != Complex::zero);
						cn.returnToCache(res.p->e[i].w);
					}
					if (res.p->e[i].w.approximatelyOne()) {
						res.p->e[i].w = Complex::one;
					}

					if (mode == 2) {
						auto c = cn.getCached();
						ComplexNumbers::div(c, res.p->e[i].w, max_value);
						res.p->e[i].map = mapdiv(res.p->e[i].map, res.map);
						cn.mul(res.p->e[i].map->extra_phase, res.p->e[i].map->extra_phase, c);
						cn.returnToCache(c);
						res.p->e[i].w = Complex::one;

					}
					else {
						auto c = cn.getTemporary();
						ComplexNumbers::div(c, res.p->e[i].w, max_value);
						auto angle = ComplexNumbers::arg(c);
						c.r->value = sqrt(ComplexNumbers::mag2(c));
						c.i->value = 0;
						res.p->e[i].w = cn.lookup(c);
						res.p->e[i].map = mapdiv(res.p->e[i].map, res.map);
						cn.mul(res.p->e[i].map->extra_phase, res.p->e[i].map->extra_phase, cn.getTemporary(cos(angle), sin(angle)));
					}
				}
			}

			res.map = append_new_map(res.map, res.p->v, add_x, cn.lookup(res.p->e[1].map->extra_phase));
			if (!isZero[1]) {
				cn.returnToCache(res.p->e[1].map->extra_phase);
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

			e = normalize(e, cached);

			assert(e.p->v == var || e.isTerminal());

			// look it up in the unique tables
			auto l = uniqueTable.lookup(e, false);

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

		the_maps* append_new_map(the_maps* self, short level, bool x, Complex rotate) {

			if (x == 0 && rotate==Complex::one) {
				return self;
			}

			std::string new_key = std::to_string(level) + "_" + std::to_string(x) + "_" + std::to_string(std::hash<dd::Complex>{}(rotate));

			auto it = self->next.find(new_key);

			if (it != self->next.end()) {
				return self->next[new_key];
			}
			else {
				self->next[new_key] = new the_maps{ level, x, rotate,Complex::one,{}, self };
				//std::cout << 570 << " " << x << " " << rotate << std::endl;
				return self->next[new_key];
			}
		}


		ComputeTable3 <the_maps*, the_maps*, the_maps*>  mapmulTable{};

		the_maps* mapmul(the_maps* self, the_maps* other) {

			if (self->level == -1) {
				other->extra_phase = cn.getCached(1,0);
				return other;
			}

			if (other->level == -1) {
				self->extra_phase = cn.getCached(1, 0);
				return self;
			}

			auto r = mapmulTable.lookup(self, other);
			if (r != nullptr) {
				r->extra_phase = cn.getCached(CTEntry::val(r->extra_phase.r), CTEntry::val(r->extra_phase.i));
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
				auto rotate = cn.getTemporary();
				if (other->x == 0) {
					cn.mul(rotate, other->rotate, self->rotate);
				}
				else {
					cn.div(rotate, other->rotate, self->rotate);
				}

				res = append_new_map(r, self->level, (self->x + other->x) % 2, cn.lookup(rotate));
				res->extra_phase = r->extra_phase;
				if (other->x) {
					cn.mul(res->extra_phase, res->extra_phase, self->rotate);
				}
			}

			//auto temp = res->extra_phase;
			//res->extra_phase = cn.lookup(res->extra_phase);
			mapmulTable.insert(self, other, res, cn.lookup(res->extra_phase));
			//res->extra_phase = temp;

			return res;
		}

		ComputeTable3 <the_maps*, the_maps*, the_maps*>  mapdivTable{};

		the_maps* mapdiv(the_maps* self, the_maps* other) {

			if (other->level == -1) {
				self->extra_phase = cn.getCached(1, 0);
				return self;
			}
			if (self == other) {
				auto the_maps_header = the_maps::the_maps_header();
				the_maps_header->extra_phase = cn.getCached(1, 0);
				return the_maps_header;
			}
			
			auto r = mapdivTable.lookup(self, other);
			if (r != nullptr) {
				r->extra_phase = cn.getCached(CTEntry::val(r->extra_phase.r), CTEntry::val(r->extra_phase.i));
				return r;
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
						auto temp = cn.getTemporary();
						cn.div(temp, Complex::one, other->rotate);
						res = append_new_map(r, other->level, other->x, cn.lookup(temp));
					}
					else {
						res = append_new_map(r, other->level, other->x, cn.conj(other->rotate));
					}
					
					res->extra_phase = r->extra_phase;
				}
				else {
					res = append_new_map(r, other->level, other->x, other->rotate);
					res->extra_phase= r->extra_phase;
					cn.div(res->extra_phase, res->extra_phase, other->rotate);
				}
			}
			else {
				auto r = mapdiv(self->father, other->father);

				bool x = (self->x + other->x) % 2;

				auto rotate = cn.getTemporary();

				if (x==1) {
					cn.mul(rotate, self->rotate, other->rotate);
				}
				else {
					cn.div(rotate, self->rotate, other->rotate);
				}
				res = append_new_map(r, self->level, x, cn.lookup(rotate));
				res->extra_phase = r->extra_phase;
				if (x == 1) {
					cn.div(res->extra_phase, res->extra_phase, other->rotate);
				}
			}
			//auto temp = res->extra_phase;
			//res->extra_phase = cn.lookup(res->extra_phase);
			mapdivTable.insert(self, other, res, cn.lookup(res->extra_phase));
			//res->extra_phase = temp;
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
				cn.returnToCache(res.e.w);
				res.e.w = cn.lookup(res.e.w);
			}

			[[maybe_unused]] const auto after = cn.cacheCount();

			assert(before == after);

			
			//the_maps::print_maps(tdd1.e.map);
			//the_maps::print_maps(tdd2.e.map);
			//the_maps::print_maps(res.e.map);
			//std::cout << tdd1.e.w << " " << tdd2.e.w <<" "<<res.e.w << std::endl;
			// std::cout << size(res.e)<<" "<< res.e.p->v << std::endl;
			// std::cout << "-----------" << std::endl;
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
						temp.map = mapmul(e.map, temp.map);

						assert(temp.w != Complex::zero);
						cn.mul(temp.w, temp.w, temp.map->extra_phase);
						cn.returnToCache(temp.map->extra_phase);
					}
					return temp;
				}
				else if (e.map->x == 0) {
					auto temp = e.p->e[c];
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = mapmul(e.map->father, temp.map);

						assert(temp.w != Complex::zero);
						cn.mul(temp.w, temp.w, temp.map->extra_phase);
						cn.returnToCache(temp.map->extra_phase);

						if (c == 1) {
								assert(temp.w != Complex::zero);
								cn.mul(temp.w, temp.w, e.map->rotate);
						}

					}
					return temp;
				}
				else {
					auto temp = e.p->e[1 - c];
					if (temp.w != Complex::zero) {
						temp.w = cn.mulCached(temp.w, e.w);
						temp.map = mapmul(e.map->father, temp.map);

						assert(temp.w != Complex::zero);
						cn.mul(temp.w, temp.w, temp.map->extra_phase);
						cn.returnToCache(temp.map->extra_phase);
						if (c == 0) {

							assert(temp.w != Complex::zero);
							cn.mul(temp.w, temp.w, e.map->rotate);

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
						temp.map = mapmul(e.map, temp.map);
						temp.w=cn.mulCached(temp.w, temp.map->extra_phase);
						cn.returnToCache(temp.map->extra_phase);
					}
					return temp;
				}
				else if (e.map->x == 0) {
					auto temp = e.p->e[c];
					if (temp.w != Complex::zero) {
						temp.map = mapmul(e.map->father, temp.map);
						temp.w = cn.mulCached(temp.w, temp.map->extra_phase);
						cn.returnToCache(temp.map->extra_phase);
						if (c == 1) {
							assert(temp.w != Complex::zero);
							cn.mul(temp.w, temp.w, e.map->rotate);
						}
					}
					return temp;
				}
				else {
					auto temp = e.p->e[1-c];
					if (temp.w != Complex::zero) {
						temp.map = mapmul(e.map->father, temp.map);
						temp.w = cn.mulCached(temp.w, temp.map->extra_phase);
						cn.returnToCache(temp.map->extra_phase);
						if (c == 0) {
							assert(temp.w != Complex::zero);
							cn.mul(temp.w, temp.w, e.map->rotate);
						}
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

			//std::cout <<"879 " << x.w << " " << y.w << " " << int(x.p->v) << " " << int(y.p->v)<<" " << x.map << " " << y.map << std::endl;
			//the_maps::print_maps(x.map);
			//the_maps::print_maps(y.map);

			if (x.p > y.p) {
				return T_add2(y, x);
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
				//std::cout << "Case 0" << std::endl;
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


			xCopy.w = Complex::one;
			xCopy.map = the_maps::the_maps_header();
			yCopy.w = cn.divCached(y.w, x.w);
			yCopy.map = mapdiv(y.map, x.map);
			if (yCopy.w != Complex::zero) {
				cn.mul(yCopy.w, yCopy.w, yCopy.map->extra_phase);
			}
			cn.returnToCache(yCopy.map->extra_phase);


			auto r = addTable.lookup({ xCopy.p, xCopy.w,xCopy.map }, { yCopy.p, yCopy.w,yCopy.map });

			if (r.p != nullptr) {
				//std::cout << "Case 1" << std::endl;
				//assert(yCopy.w != Complex::zero);
				cn.returnToCache(yCopy.w);
				if (r.w.approximatelyZero()) {
					return Edge<Node>::zero;
				}
				auto c = cn.getCached(r.w);

				if (c != Complex::zero) {
					cn.mul(c, c, x.w);
				}

				auto temp_map = mapmul(x.map, r.map);
				if (c != Complex::zero) {
					cn.mul(c, c, temp_map->extra_phase);
				}
				cn.returnToCache(temp_map->extra_phase);
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
				cn.mul(e.w, e.w, e.map->extra_phase);
				cn.returnToCache(e.map->extra_phase);
			}

			cn.returnToCache(yCopy.w);
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
				res->remain_map->extra_phase = temp_pahse;
				return res;
			}
			if (newk1 < newk2 && !ifContract(newk2)) {
				auto res = find_remain_map(map1, map2->father, temp_key_2_new_key1, temp_key_2_new_key2);
				auto temp_pahse = res->remain_map->extra_phase;
				res->remain_map = append_new_map(res->remain_map, newk2, map2->x, map2->rotate);
				res->remain_map->extra_phase = temp_pahse;
				return res;
			}
			if (map1->level == -1 && map2->level == -1) {
				comm_maps* res=new comm_maps{ the_maps::the_maps_header(),the_maps::the_maps_header(),the_maps::the_maps_header() };
				res->remain_map->extra_phase = cn.getCached(1,0);
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
				assert(res->remain_map->extra_phase != Complex::zero);
				cn.mul(res->remain_map->extra_phase, res->remain_map->extra_phase, map2->rotate);
			}

			auto rotate = cn.getTemporary();
			if (x == 0) {
				assert(rotate != Complex::zero);
				cn.mul(rotate, map1->rotate, map2->rotate);
			}
			else {
				cn.div(rotate, map1->rotate, map2->rotate);
			}

			res->cont_map1 = append_new_map(res->cont_map1, map1->level, x, cn.lookup(rotate));

			return res;
		}

		//template <class LeftOperandNode, class RightOperandNode>
		Edge<mNode> cont2(const Edge<mNode>& x, const Edge<mNode>& y, key_2_new_key_node* key_2_new_key1, key_2_new_key_node* key_2_new_key2, const int var_num) {
			// auto& id = this->identity;

			//std::cout <<"838 " << x.w << " " << y.w << " " << int(x.p->v) << " " << int(y.p->v) << std::endl;
			//the_maps::print_maps(x.map);
			//the_maps::print_maps(y.map);

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
					assert(c != Complex::zero);
					ComplexNumbers::mul(c, c, cn.getTemporary(pow(2, var_num), 0));
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

			xCopy.map = r_maps->cont_map1;
			yCopy.map = r_maps->cont_map2;
			yCopy.map->print_maps(yCopy.map);
			//auto extra_phase = cn.getCached(r_maps->remain_map->extra_phase.r->value, r_maps->remain_map->extra_phase.i->value);
			auto extra_phase = r_maps->remain_map->extra_phase;

			auto res = contTable.lookup(xCopy, yCopy, temp_key_2_new_key1, temp_key_2_new_key2);
			if (res.e.p != nullptr) {
				if (res.e.w.approximatelyZero()) {
					cn.returnToCache(extra_phase);
					return ResultEdge::zero;
				}
				auto e = ResultEdge{ res.e.p, cn.getCached(res.e.w),res.e.map };
				assert(e.w != Complex::zero);
				ComplexNumbers::mul(e.w, e.w, x.w);
				ComplexNumbers::mul(e.w, e.w, y.w);
				if (e.w.approximatelyZero()) {
					//assert(e.w != Complex::zero);
					cn.returnToCache(e.w);
					cn.returnToCache(extra_phase);
					return ResultEdge::zero;
				}
				//std::cout << "1160 " << var_num << " " << res.cont_num << std::endl;
				if (res.cont_num != var_num) {
					assert(e.w != Complex::zero);
					ComplexNumbers::mul(e.w, e.w, cn.getTemporary(pow(2, var_num - res.cont_num), 0));//对于一般形状的tensor,以2为底数可能有问题
					// TODO: pow(2,n) can be optimized by 1<<n
				}
				e.map = mapmul(r_maps->remain_map, e.map);
				assert(e.w != Complex::zero);
				cn.mul(e.w, e.w, e.map->extra_phase);
				cn.returnToCache(e.map->extra_phase);
				assert(e.w != Complex::zero);
				cn.mul(e.w, e.w, extra_phase);
				cn.returnToCache(extra_phase);
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
			if(to_test){
				std::cout << 1298 << std::endl;
				std::cout << "newk1: " << newk1 << " newk2: " << newk2 << std::endl;
			}
			ResultEdge e1{}, e2{}, r{};

			if (newk1 > newk2) {
				// TODO: half integer?
				// bool isHalfInteger = std::abs(newk1 - std::round(newk1)) > 0.4999999 && std::abs(newk1 - std::round(newk1)) < 0.5000001;
				if (ifContract(newk1)) {
					r = ResultEdge::zero;
					ResultEdge etemp{};
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						e1 = Slicing2(xCopy, xCopy.p->v, k);
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
						e1 = Slicing2(xCopy, xCopy.p->v, k);
						e2 = yCopy;
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						if (e1.w != Complex::zero) {
							cn.returnToCache(e1.w);
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
					r = makeDDNode(Qubit(newk1), e, true);
				}
			}
			else if (newk1 < newk2) {
				if (ifContract(newk2)) {
					r = ResultEdge::zero;
					ResultEdge etemp{};
					for (int k = 0; k < y.p->e.size(); ++k) {
						e1 = xCopy;
						//e2 = y.p->e[k];
						e2 = Slicing2(yCopy, yCopy.p->v, k);
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
						e2 = Slicing2(yCopy, yCopy.p->v, k);
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						if (e2.w != Complex::zero) {
							cn.returnToCache(e2.w);
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
					r = makeDDNode(Qubit(newk2), e, true);
				}

			}
			else {
				if(to_test){
				std:: cout << 1410 << int(newk2 * 2) % 2 << std::endl; 

				}
				if (ifContract(newk2)) {
					r = ResultEdge::zero;
					ResultEdge etemp{};
					for (int k = 0; k < x.p->e.size(); ++k) {
						//e1 = x.p->e[k];
						e1 = Slicing2(xCopy, xCopy.p->v, k);
						//e2 = y.p->e[k];
						e2 = Slicing2(yCopy, yCopy.p->v, k);
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
						e1 = Slicing2(xCopy, xCopy.p->v, k);
						//e2 = y.p->e[k];
						e2 = Slicing2(yCopy, yCopy.p->v, k);
						e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
						if (e1.w != Complex::zero) {
							cn.returnToCache(e1.w);
						}
						if (e2.w != Complex::zero) {
							cn.returnToCache(e2.w);
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
					r = makeDDNode(Qubit(newk1), e, true);
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
			contTable.insert(xCopy, yCopy, { r.p, r.w,r.map }, temp_key_2_new_key1, temp_key_2_new_key2, var_num);
			
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
					cn.returnToCache(r.w);
					cn.returnToCache(extra_phase);
					return ResultEdge::zero;
				}
			}
			if (r.w == Complex::zero) {
				cn.returnToCache(extra_phase);
				return ResultEdge::zero;
			}
			else {
				r.map = mapmul(r_maps->remain_map, r.map);
				assert(r.w != Complex::zero);
				cn.mul(r.w, r.w, r.map->extra_phase);
				cn.returnToCache(r.map->extra_phase);
				assert(r.w != Complex::zero);
				cn.mul(r.w, r.w, extra_phase);
				cn.returnToCache(extra_phase);
			}

			
			//std::cout << "Case 2 " << r.w << " " << int(r.p->v) << std::endl;
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
