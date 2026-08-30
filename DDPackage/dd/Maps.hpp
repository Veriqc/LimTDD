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
		short level;
		bool x;
		// Complex rotate;// rotate始终是一个complexTable里的元素，在中间计算过程，可以在temporary里面
		// Complex extra_phase;// rotate始终是一个temporary里的元素
		int rotate;// rotate始终是一个complexTable里的元素，在中间计算过程，可以在temporary里面

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
		int remain_phase;

	};

	// Return type of mapmul / mapdiv: the structural map plus the "pending global
	// phase" (overflow) produced by the combination, in units of rotate_angle (π/4).
	struct map_res {
		the_maps* map;
		int phase;
	};


}

namespace std {
	// Content-based hash for the_maps*: hash the (level, x, rotate) chain so that
	// structurally-identical maps (deduplicated by append_new_map's string key)
	// hash deterministically regardless of the ASLR-dependent node addresses. This
	// is exact (not tolerance-rounded), so it stays consistent with the interning.
	template <> struct hash<dd::the_maps*> {
		std::size_t operator()(dd::the_maps* m) const noexcept {
			std::size_t key = 0;
			for (const dd::the_maps* cur = m; cur != nullptr; cur = cur->father) {
				key = dd::combineHash(key, dd::murmur64(static_cast<std::size_t>(cur->level)));
				key = dd::combineHash(key, dd::murmur64(static_cast<std::size_t>(cur->x ? 1 : 0)));
				key = dd::combineHash(key, dd::murmur64(static_cast<std::size_t>(cur->rotate)));
				if (cur->level == -1) {
					break;
				}
			}
			return key;
		}
	};
} // namespace std
