#include "dd/Maps.hpp"

#include <math.h>

#include "dd/ComputeTable.hpp"

namespace dd {

	the_maps the_maps::the_maps_header_element{ -1, 0, 0, 0, {}, nullptr };
	//the_maps* the_maps_header = &the_maps_header_element;

	//the_maps* get_the_maps_header() {
	//	return the_maps_header;
	//}

	the_maps* the_maps::append_new_map(the_maps* self, short level, short x, int rotate) {

		rotate = rotate % root_of_unit;

		//std::cout << rotate << " " << rotate % root_of_unit << std::endl;
		if (rotate < 0) {
			rotate += root_of_unit;
		}

		if (x == 0 && rotate == 0) {
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

	the_maps* the_maps::mapmul(the_maps* self, the_maps* other) {

		//std::cout << self->level << " " << other->level << std::endl;

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
			//long int rotate = other->rotate + self->rotate * pow(-1, other->x);
			auto rotate = other->rotate;
			if (other->x == 0) {
				rotate += self->rotate;
			}
			else {
				rotate -= self->rotate;
			}

			rotate = rotate % root_of_unit;
			res = append_new_map(r, self->level, (self->x + other->x) % 2, rotate);
			res->extra_phase = r->extra_phase + self->rotate * other->x;
		}
		mapmulTable.insert(self, other, res);
		return res;
	}
	ComputeTable3 <the_maps*, the_maps*, the_maps*>  mapdivTable{};

	the_maps* the_maps::mapdiv(the_maps* self, the_maps* other) {

		//std::cout << self->level << "  " << other->level << std::endl;

		if (other->level == -1) {
			self->extra_phase = 0;
			return self;
		}
		if (self == other) {
			auto the_maps_header = the_maps::the_maps_header();
			the_maps_header->extra_phase = 0;
			//assert(the_maps::the_maps_header()->extra_phase == 0);
			return the_maps_header;
		}

		auto r = mapdivTable.lookup(self, other);
		if (r != nullptr) {
			return r;
		}

		the_maps* res;
		if (self->level > other->level) {
			auto r = the_maps::mapdiv(self->father, other);
			res = append_new_map(r, self->level, self->x, self->rotate);
			res->extra_phase = r->extra_phase;
		}
		else if (self->level < other->level) {
			auto r = the_maps::mapdiv(self, other->father);
			//long int rotate = other->rotate * pow(-1, 1 - other->x);
			auto rotate = other->rotate;
			if (1 - other->x == 1) {
				rotate *= -1;
			}
			rotate = rotate % root_of_unit;
			res = append_new_map(r, other->level, other->x, rotate);
			res->extra_phase = r->extra_phase - other->rotate * other->x;
		}
		else {
			auto r = the_maps::mapdiv(self->father, other->father);
			r->extra_phase += self->rotate * other->x;
			short x = (self->x + other->x) % 2;
			//long int rotate = self->rotate + other->rotate * pow(-1, 1 - x);
			auto rotate = self->rotate;
			if (1 - x == 0) {
				rotate += other->rotate;
			}
			else {
				rotate -= other->rotate;
			}
			rotate = rotate % root_of_unit;
			res = append_new_map(r, self->level, x, rotate);
			res->extra_phase = r->extra_phase - other->rotate * x;
		}
		mapdivTable.insert(self, other, res);
		return res;
	}

	void the_maps::print_maps(the_maps* map) {
		if (map->level == -1) {
			std::cout<<"  ." << std::endl;
		}
		else {
			std::cout << map->level << ":";
			if (map->x == 1) {
				std::cout << "x ";
			}
			if (map->rotate != 0) {
				std::cout << map->rotate;
			}
			std::cout << ";";
			print_maps(map->father);
		}
	}

	std::string the_maps::to_string(the_maps* map) {
		if (map->level == -1) {
			return "  .";
		}
		else {
			std::string s = "";
			s += std::to_string(map->level);
			s+=":";
			if (map->x == 1) {
				s+= "x ";
			}
			if (map->rotate != 0) {
				s+= std::to_string(map->rotate);
			}
			s+= ";";
			s+=to_string(map->father);
			return s;
		}
	}


} // namespace dd