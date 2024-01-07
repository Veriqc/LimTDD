#pragma once

#include "ComplexValue.hpp"
#include "Definitions.hpp"
#include "xtensor/xarray.hpp"

#include <cmath>

namespace dd {
// Complex constants
// NOLINTBEGIN(readability-identifier-naming) As these constants are used by
// other projects, we keep the naming
ComplexValue complex_one = {1., 0.};
ComplexValue complex_mone = {-1., 0.};
ComplexValue complex_zero = {0., 0.};
ComplexValue complex_i = {0., 1.};
ComplexValue complex_mi = {0., -1.};
ComplexValue complex_SQRT2_2 = {SQRT2_2, 0.};
ComplexValue complex_mSQRT2_2 = {-SQRT2_2, 0.};
ComplexValue complex_iSQRT2_2 = {0., SQRT2_2};
ComplexValue complex_miSQRT2_2 = {0., -SQRT2_2};
ComplexValue complex_1plusi = {SQRT2_2, SQRT2_2};
ComplexValue complex_1minusi = {SQRT2_2, -SQRT2_2};
ComplexValue complex_1plusi_2 = {0.5, 0.5};
ComplexValue complex_1minusi_2 = {0.5, -0.5};

// Gate matrices
using GateMatrix = xt::xarray<dd::ComplexValue>;

GateMatrix Imat{
		{complex_one, complex_zero},
		{complex_zero, complex_one}
		};
GateMatrix Hmat{
  {complex_SQRT2_2, complex_SQRT2_2},
  {complex_SQRT2_2, complex_mSQRT2_2}
  };
GateMatrix Xmat{
  {complex_zero, complex_one},
  {complex_one, complex_zero}
  };
GateMatrix Ymat{
  {complex_zero, complex_mi},
  {complex_i, complex_zero}
  };
GateMatrix Zmat{
  {complex_one, complex_zero},
  {complex_zero,complex_mone}
  };
GateMatrix Smat{
  {complex_one, complex_zero},
  {complex_zero, complex_i}
  };
GateMatrix Sdagmat{
  {complex_one, complex_zero},
  {complex_zero,complex_mi}
  };
GateMatrix Tmat{
  {complex_one, complex_zero},
  {complex_zero,complex_1plusi}
  };
GateMatrix Tdagmat{
  {complex_one, complex_zero},
  {complex_zero,complex_1minusi}
  };
GateMatrix SXmat{
      {complex_1plusi_2, complex_1minusi_2},
      {complex_1minusi_2, complex_1plusi_2}
      };
GateMatrix SXdagmat{
      {complex_1minusi_2, complex_1plusi_2},
      {complex_1plusi_2, complex_1minusi_2}
      };
GateMatrix Vmat{
      {complex_SQRT2_2, complex_miSQRT2_2},
      {complex_miSQRT2_2,complex_SQRT2_2}
      };
GateMatrix Vdagmat{
      {complex_SQRT2_2, complex_iSQRT2_2},
      {complex_iSQRT2_2, complex_SQRT2_2}
      };
GateMatrix CXmat{
      {
          {
              {complex_one, complex_zero},
              {complex_zero, complex_one}
          },
          {
              {complex_zero, complex_zero},
              {complex_zero, complex_zero}
          }
      },
      {
          {
              {complex_zero, complex_zero},
              {complex_zero, complex_zero}
          },
          {
              {complex_zero, complex_one},
              {complex_one, complex_zero}
          }
      }
  };

inline GateMatrix U3mat(double lambda, double phi, double theta) {
return GateMatrix{
          {
              {std::cos(theta / 2.),
              0.},
              {-std::cos(lambda) * std::sin(theta / 2.),
              -std::sin(lambda) * std::sin(theta / 2.)}
          },
          {
              {std::cos(phi) * std::sin(theta / 2.),
              std::sin(phi) * std::sin(theta / 2.)},
              {std::cos(lambda + phi) * std::cos(theta / 2.),
              std::sin(lambda + phi) * std::cos(theta / 2.)}
          }
      };
}

inline GateMatrix Phasemat(double lambda) {
return GateMatrix{
          {complex_one,complex_zero},
          {
              complex_zero,
              {std::cos(lambda),
              std::sin(lambda)
              }
          }
      };
}
}// namespace dd
