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

  Tensor<double,2> I4(void)
  {
    Tensor<double,2> I((TensorShape<2>())[4][4],0.);
    for(int i =0; i < 4; ++i)
      I[i][i] = 1.;
    I.lower_index(1);
    return I;
  }

  bool isRest(const XMLArray::Array<int> &mom)
  {
    return ((mom[0] == 0) && (mom[1] == 0) && (mom[2] == 0));
  }


  // active z-y-z rotations
  Tensor<double,2> genRotationMatrix(const XMLArray::Array<int> &mom)
  {
    if(isRest(mom))
      return I4();

    Hadron::CubicCanonicalRotation_t eulerangles = Hadron::cubicCanonicalRotation(mom);

    Tensor<double,2> A((TensorShape<2>())[4][4],0.),B((TensorShape<2>())[4][4],0.),C((TensorShape<2>())[4][4],0.);


    double a,b,g;
    a = eulerangles.alpha;
    b = eulerangles.beta;
    g = eulerangles.gamma;


    A[0][0] = 1.;
    B[0][0] = 1.;
    C[0][0] = 1.;


    A[1][1] = cos(a)   ;    A[1][2] = -sin(a)    ;    A[1][3] =  0.      ; 
    A[2][1] = sin(a)  ;    A[2][2] = cos(a)    ;    A[2][3] =  0.      ; 
    A[3][1] = 0.       ;    A[3][2] = 0.        ;    A[3][3] =  1.      ; 


#if 0
    // B_x
    B[1][1] = 1.       ;    B[1][2] = 0.        ;    B[1][3] =  0.      ; 
    B[2][1] = 0.       ;    B[2][2] = cos(b)    ;    B[2][3] = sin(b)   ; 
    B[3][1] = 0.       ;    B[3][2] = -sin(b)   ;    B[3][3] = cos(b)   ; 
#endif 

    // B_y 
    B[1][1] = cos(b)   ;    B[1][2] = 0.        ;    B[1][3] = sin(b)   ; 
    B[2][1] = 0.       ;    B[2][2] = 1.        ;    B[2][3] =  0.      ; 
    B[3][1] = -sin(b)  ;    B[3][2] = 0.        ;    B[3][3] = cos(b)   ; 



    C[1][1] = cos(g)   ;    C[1][2] = -sin(g)    ;    C[1][3] =  0.      ; 
    C[2][1] = sin(g)  ;    C[2][2] = cos(g)    ;    C[2][3] =  0.      ; 
    C[3][1] = 0.       ;    C[3][2] = 0.        ;    C[3][3] =  1.      ; 



#if 0

    for(idx_t i = 0; i < 4; ++i)
      A[i][i] = 1.;

    B = A;
    C = A;


    A[1][1] = A[2][2] = cos(eulerangles.gamma);
    A[2][1] = sin(eulerangles.gamma);
    A[1][2] = -A[2][1];

    B[1][1] = B[3][3] = cos(eulerangles.beta);
    B[3][1] = sin(eulerangles.beta);
    B[1][3] = -B[3][1];

    C[1][1] = C[2][2] = cos(eulerangles.alpha);
    C[2][1] = sin(eulerangles.alpha);
    C[1][2] = -C[2][1];
#endif 

    A.lower_index(1);
    B.lower_index(1);
    C.lower_index(1);

    // return C * B * A;
    return A*B*C;
  }


  Tensor<double,2> genRotationMatrix3D(const XMLArray::Array<int> &mom)
  {
    Tensor<double,2> Four = genRotationMatrix(mom);
    Tensor<double,2> Three((TensorShape<2>())[3][3],0.);
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        Three[i][j] = Four[i+1][j+1];

    Three.lower_index(1); 

    return Three;
  }

}
