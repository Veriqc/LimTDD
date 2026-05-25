#pragma once

#include "Complex.hpp"
#include "ComplexValue.hpp"
#include "Definitions.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <unordered_map>
#include <functional>

namespace dd {

	struct the_maps {
		const short level;
		const bool x;
		// Complex rotate;// rotate始终是一个complexTable里的元素，在中间计算过程，可以在temporary里面
		// Complex extra_phase;// rotate始终是一个temporary里的元素
		const int rotate;// rotate始终是一个complexTable里的元素，在中间计算过程，可以在temporary里面
		const int extra_phase;// rotate始终是一个temporary里的元素


		struct MapKey {
			short level;
			bool x;
			int rotate;
			bool operator==(MapKey const& o) const noexcept {
				return level == o.level && x == o.x && rotate == o.rotate;
			}
		};

		struct MapKeyHash {
			size_t operator()(MapKey const& k) const noexcept {
				std::size_t h = std::hash<short>()(k.level);
				h ^= std::hash<int>()(k.rotate) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
				h ^= std::hash<bool>()(k.x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
				return h;
			}
		};

		mutable std::unordered_map<MapKey, the_maps*, MapKeyHash> next;
		the_maps* const father;

		static the_maps the_maps_header_element;

		static constexpr the_maps* the_maps_header() { return &the_maps_header_element; }

		[[nodiscard]] static int normalize_phase(int phase);

		//static the_maps* mapdiv(the_maps* self, the_maps* other);

		//static the_maps* mapmul(the_maps* self, the_maps* other);

		//static the_maps* append_new_map(the_maps* self, short level, bool x, Complex rotate);

		static void print_maps(the_maps* map);

		static std::string to_string(the_maps* map);

		//static the_maps** find_remain_map(the_maps* map1, the_maps* map2, key_2_new_key_node* key_2_new_key);

	};

	struct comm_maps {
		the_maps* remain_map = the_maps::the_maps_header();
		the_maps* cont_map1 = the_maps::the_maps_header();
		the_maps* cont_map2 = the_maps::the_maps_header();
		int phase = 0;

	};


}