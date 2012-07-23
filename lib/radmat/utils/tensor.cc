// tensor.cc -
//
// Thursday, May 10 2012
//

#include "tensor.h"
#include "hadron/irrep_util.h"

namespace radmat
{

  Tensor<double,2> g_dd(void)
  {
    Tensor<double,2> gmunu((TensorShape<2>())[4][4],0.);
    gmunu[0][0] = 1;
    gmunu[1][1] = gmunu[2][2] = gmunu[3][3] = -1;
    gmunu.lower_index(0);
    gmunu.lower_index(1);
    return gmunu;
  }

  Tensor<double,2> g_uu(void)
  {
    Tensor<double,2> gmunu((TensorShape<2>())[4][4],0.);
    gmunu[0][0] = 1;
    gmunu[1][1] = gmunu[2][2] = gmunu[3][3] = -1;
    return gmunu;
  }

 Tensor<double,2> genRotationMatrix(const XMLArray::Array<int> &mom)
  {
    Hadron::CubicCanonicalRotation_t eulerangles = Hadron::cubicCanonicalRotation(mom);

    Tensor<double,2> A((TensorShape<2>())[4][4],0.),B((TensorShape<2>())[4][4],0.),C((TensorShape<2>())[4][4],0.);

    for(idx_t i = 0; i < 4; ++i)
      A[i][i] = B[i][i] = C[i][i] = 1.;

    A[1][1] = A[2][2] = cos(eulerangles.gamma);
    A[2][1] = sin(eulerangles.gamma);
    A[1][2] = -A[2][1];

    B[1][1] = B[3][3] = cos(eulerangles.beta);
    B[1][3] = sin(eulerangles.beta);
    B[3][1] = -B[1][3];

    C[1][1] = C[2][2] = cos(eulerangles.alpha);
    C[2][1] = sin(eulerangles.alpha);
    C[1][2] = -C[2][1];

    A.lower_index(1);
    B.lower_index(1);
    C.lower_index(1);

    return C*B*A;
  }

}
