#pragma once

#include "Complex.hpp"
#include "ComplexValue.hpp"
#include "ComplexNumbers.hpp"
#include "Definitions.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace dd {

	struct the_maps {
		short level;
		bool x;
		Complex rotate;// rotate始终是一个complexTable里的元素，在中间计算过程，可以在temporary里面
		Complex extra_phase;// rotate始终是一个temporary里的元素

		std::map<std::string, the_maps*> next;
		the_maps* father;

		static the_maps the_maps_header_element;

		static constexpr the_maps* the_maps_header() { return &the_maps_header_element; }

		//static the_maps* mapdiv(the_maps* self, the_maps* other);

		//static the_maps* mapmul(the_maps* self, the_maps* other);

		//static the_maps* append_new_map(the_maps* self, short level, bool x, Complex rotate);

		static void print_maps(the_maps* map);

		static std::string to_string(the_maps* map);

		//static the_maps** find_remain_map(the_maps* map1, the_maps* map2, key_2_new_key_node* key_2_new_key);

	};

	struct comm_maps {
		the_maps* remain_map;
		the_maps* cont_map1;
		the_maps* cont_map2;

	};

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


}