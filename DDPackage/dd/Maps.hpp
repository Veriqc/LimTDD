#pragma once

#include "Complex.hpp"
#include "ComplexValue.hpp"
#include "Definitions.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace dd {

	struct the_maps {
		short level;
		short x;
		int rotate;
		int extra_phase;

		std::map<std::string, the_maps*> next;
		the_maps* father;

		static the_maps the_maps_header_element;

		static constexpr the_maps* the_maps_header() { return &the_maps_header_element; }

		static the_maps* mapdiv(the_maps* self, the_maps* other);

		static the_maps* mapmul(the_maps* self, the_maps* other);

		static the_maps* append_new_map(the_maps* self, short level, short x, int rotate);

		static void print_maps(the_maps* map);

		static std::string to_string(the_maps* map);

		//static the_maps** find_remain_map(the_maps* map1, the_maps* map2, key_2_new_key_node* key_2_new_key);

	};

	struct comm_maps {
		the_maps* remain_map;
		the_maps* cont_map1;
		the_maps* cont_map2;

	};

}