#include "dd/Edge.hpp"

#include "dd/Node.hpp"

#include "dd/ComputeTable.hpp"

#include <math.h>

namespace dd {

	the_maps the_maps_header_element={ -1,0,0,0,{},nullptr };
	the_maps* the_maps_header = &the_maps_header_element;

	the_maps* get_the_maps_header() {
		return the_maps_header;
	}

	the_maps* append_new_map(the_maps* self, short level, short x, long int rotate) {

		if(x==0 && rotate==0){
			return self;
		}

		std::string new_key = std::to_string(level) + "_" + std::to_string(x) + "_" + std::to_string(rotate);

		auto it = self->next.find(new_key);

		if (it != self->next.end()) {
			return self->next[new_key];
		}
		else {
			self->next[new_key] = new the_maps{ level, x, rotate,0,{}, self };

			return self->next[new_key];
		}
	}

	ComputeTable3 <the_maps*, the_maps*, the_maps*>  mapmulTable{};

	the_maps* mapmul(the_maps* self, the_maps* other) {
		if (self->level == -1) {
			other->extra_phase = 0;
			return other;
		}

		if (other->level == -1) {
			self->extra_phase = 0;
			return self;
		}

		auto r = mapmulTable.lookup(self, other);
		if (r != nullptr) {
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
			long int rotate = other->rotate + self->rotate * pow(-1, other->x);
			rotate = rotate % root_of_unit;
			res = append_new_map(r, self->level, (self->x + other->x) % 2, rotate);
			res->extra_phase = r->extra_phase + self->rotate * other->x;
		}
		mapmulTable.insert(self, other, res);
		return res;
	}
	ComputeTable3 <the_maps*, the_maps*, the_maps*>  mapdivTable{};

	the_maps* mapdiv(the_maps* self, the_maps* other) {

		if (other->level == -1) {
			self->extra_phase = 0;
			return self;
		}
		if (self == other) {
			the_maps_header->extra_phase = 0;
			return the_maps_header;
		}

		auto r = mapdivTable.lookup(self, other);
		if (r != nullptr) {
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
			long int rotate = other->rotate * pow(-1, 1 - other->x);
			rotate = rotate % root_of_unit;
			res = append_new_map(r, other->level, other->x, rotate);
			res->extra_phase = r->extra_phase - other->rotate * other->x;
		}
		else {
			auto r = mapdiv(self->father, other->father);
			r->extra_phase += self->rotate * other->x;
			short x = (self->x + other->x) % 2;
			long int rotate = self->rotate + other->rotate * pow(-1, 1 - x);
			rotate = rotate % root_of_unit;
			res = append_new_map(r, self->level, x, rotate);
			res->extra_phase = r->extra_phase - other->rotate * x;
		}
		mapdivTable.insert(self, other, res);
		return res;
	}

	void print_maps(the_maps* map) {

		return;
	}


} // namespace dd

namespace dd {

	template <class Node>
	const Edge<Node> Edge<Node>::zero{Node::getTerminal(), Complex::zero, the_maps_header};
	template <class Node>
	const Edge<Node> Edge<Node>::one{Node::getTerminal(), Complex::one, the_maps_header};

	template <class Node> Edge<Node> Edge<Node>::terminal(const Complex& w) {
		return { Node::getTerminal(), w ,the_maps_header };
	}

	template <class Node> bool Edge<Node>::isTerminal() const {
		return Node::isTerminal(p);
	}

	template <class Node> bool Edge<Node>::isZeroTerminal() const {
		return isTerminal() && w == Complex::zero;
	}

	template <class Node> bool Edge<Node>::isOneTerminal() const {
		return isTerminal() && w == Complex::one;
	}

	// Explicit instantiations

	template struct Edge<mNode>;


} // namespace dd


