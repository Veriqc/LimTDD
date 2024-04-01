#include <map>
#include <string>
#include <iostream>
#include <functional>
#include <ctime>


#include "dd/Export.hpp"
#include "dd/Tensor.hpp"
#include "dd/Package.hpp"
#include <xtensor/xio.hpp>
#include <xtensor/xarray.hpp>

int main(){
    dd::ComplexValue one = { 1,0 };
    dd::ComplexValue zero = { 0,0 };
    auto dd1 = std::make_unique<dd::Package<>>(5);
    dd1->varOrder = { {"x0",0},{"y0",1} };

    xt::xarray<dd::ComplexValue> U = {{one, zero}, {zero, one}};
    std::vector<dd::Index> indexs = {{"x0",0},{"y0",1}};

    dd::Tensor tn = { U,indexs, "I" };
    dd::TDD tdd = tn.to_tdd(dd1.get());


  //   std::cout << tdd.e.p->v << std::endl;
  //   std::cout << tdd.e.w << std::endl;
  //   std::cout << tdd.e.map->level << std::endl;
  //   std::cout << tdd.e.map->x << std::endl;
  //   tdd.e.map->print_maps(tdd.e.map);
  //   for (auto & i: tdd.index_set) {
	// 	std::cout << "(" << i.key << ", " << i.idx << ") ,";
	// }
	// std::cout << std::endl;
  //   for (auto & i: tdd.key_2_index) {
	// 	std::cout << i << " ";
	// }
	std::cout << std::endl;
  dd::export2Dot(dd1->identity, "identity",true,true);
}