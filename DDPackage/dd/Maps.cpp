#include "dd/Maps.hpp"

#include <math.h>

#include "dd/ComputeTable.hpp"


namespace dd {

	the_maps the_maps::the_maps_header_element{ -1, 0, 0, 0, {}, nullptr };

	int the_maps::normalize_phase(int phase) {
		return (phase % root_of_unit + root_of_unit) % root_of_unit;
	}

	void the_maps::print_maps(the_maps* map) {
		if (map->level == -1) {
			if (map->extra_phase != 0) {
				std::cout << "  .@" << map->extra_phase << std::endl;
			} else {
				std::cout<<"  ." << std::endl;
			}
		}
		else {
			std::cout << map->level << ":";
			if (map->x) {
				std::cout << "x ";
			}
			if (map->rotate != 0) {
				std::cout << map->rotate;
			}
			if (map->extra_phase != 0) {
				std::cout << "@" << map->extra_phase;
			}
			std::cout << ";";
			print_maps(map->father);
		}
	}

	std::string the_maps::to_string(the_maps* map) {
		if (map->level == -1) {
			if (map->extra_phase != 0) {
				return "  .@" + std::to_string(map->extra_phase);
			}
			return "  .";
		}
		else {
			std::string s = "";
			s += std::to_string(map->level);
			s+=":";
			if (map->x) {
				s+= "x ";
			}
			if (map->rotate != 0) {
				s+= std::to_string(map->rotate);
			}
			if (map->extra_phase != 0) {
				s += "@";
				s += std::to_string(map->extra_phase);
			}
			s+= ";";
			s+=to_string(map->father);
			return s;
		}
	}


} // namespace dd