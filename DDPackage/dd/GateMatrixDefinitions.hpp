#pragma once

#include "ComplexValue.hpp"
#include "Definitions.hpp"
#include "xtensor/xarray.hpp"
#include <xtensor/xview.hpp>

#include <cmath>

namespace dd {
// Complex constants
// NOLINTBEGIN(readability-identifier-naming) As these constants are used by
// other projects, we keep the naming
std::set<std::string> singleGate = {"X","Y","Z","S","Sdag"};
std::set<std::string> twoGate = {"cx","cu1","swap"};
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

inline GateMatrix U2mat(double phi, double lambda){
        return GateMatrix{
            {complex_SQRT2_2,{-SQRT2_2 *std::cos(lambda),-SQRT2_2 * std::sin(lambda)}},
            {{SQRT2_2*std::cos(phi), SQRT2_2*std::sin(phi)},{SQRT2_2* std::cos(lambda+phi), SQRT2_2* std::sin(lambda+phi)}}
        };
    }

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

GateMatrix Rzmat(double lambda) {
return GateMatrix{
        {{std::cos(lambda/2.),
            -1*std::sin(lambda/2.)}
            ,complex_zero},
        {
            complex_zero,
            {std::cos(lambda/2.),
            std::sin(lambda/2.)
            }
        }
    };
}

GateMatrix Rxmat(double lambda) {
return GateMatrix{
        {{std::cos(lambda/2.),0.},{0.,-1*std::sin(lambda/2)}},
        {{0.,-1*std::sin(lambda/2.)},{cos(lambda/2.),0.}}
    };
}

GateMatrix Rymat(double lambda) {
return GateMatrix{
        {{std::cos(lambda/2.),0.},{-1*std::sin(lambda/2),0.}},
        {{std::sin(lambda/2.),0.},{cos(lambda/2.),0.}}
    };
}
GateMatrix SWAPmat{
      {
          {
              {complex_one, complex_zero},
              {complex_zero, complex_zero}
          },
          {
              {complex_zero, complex_zero},
              {complex_one, complex_zero}
          }
      },
      {
          {
              {complex_zero, complex_one},
              {complex_zero, complex_zero}
          },
          {
              {complex_zero, complex_zero},
              {complex_zero, complex_one}
          }
      }
  };

GateMatrix CU1mat(double lambda){
return GateMatrix{
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
              {complex_one, complex_zero},
              {complex_zero, {std::cos(lambda), std::sin(lambda)}}
          }
      }
  };
}
GateMatrix identity(int n) {
    if (n == 1) {
        return Imat;
    }

    GateMatrix I = identity(n - 1);

    // Directly calculate the new shape without intermediate steps
    std::vector<size_t> newShape(2+ I.dimension(), 2); // Fills the new shape with 2s

    GateMatrix combined = GateMatrix::from_shape(newShape); 
    combined.fill(complex_zero);
    // Initialized to 0s
    for (size_t i = 0; i < 2; ++i) {
        xt::view(combined, i, i, xt::all(), xt::all()) = I; // Place I in the diagonal blocks
    }
    return combined;
}
GateMatrix controlMat(GateMatrix mat, int c) {
    if (c == 0) return mat;

    mat = controlMat(mat, c - 1);

    std::vector<size_t> newShape(2 + mat.dimension(), 2); 

    GateMatrix combined = GateMatrix::from_shape(newShape); 
    combined.fill(complex_zero);

    xt::view(combined, 0, 0, xt::all(), xt::all()) = identity(mat.dimension() / 2); // Place identity matrix
    xt::view(combined, 1, 1, xt::all(), xt::all()) = mat; // Place mat in the bottom-right block

    return combined;
}
}// namespace dd
