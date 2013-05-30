/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : phases2.cc

 * Purpose :

 * Creation Date : 05-05-2013

 * Last Modified : Mon May  6 12:29:47 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "hadron/clebsch.h"
#include "ensem/ensem.h"
#include "hadron/irrep_util.h"
#include "itpp/itbase.h"
#include "semble/semble_meta.h"
#include <complex>
#include "radmat/utils/tensor.h"
#include <sstream>
#include <iostream>
#include <string>

using namespace radmat;


std::string print(const itpp::Mat<std::complex<double> > &m)
{
  std::stringstream ss; 
  ss << itpp::round_to_zero(m,0.001); 
  return ss.str();
}

itpp::Vec<std::complex<double> > epsz(int m_row)
{
  itpp::Vec<std::complex<double> > ret(3); 
  ret.zeros(); 

  switch(m_row)
  {

    case 3:
      ret[0] = std::complex<double>(1./sqrt(2.),0.); 
      ret[1] = std::complex<double>(0.,-1./sqrt(2.));
      break;
    case 2:
      ret[2] = std::complex<double>(1.,0.); 
      break; 
    case 1:
      ret[0] = std::complex<double>(-1./sqrt(2.),0.); 
      ret[1] = std::complex<double>(0.,-1./sqrt(2.));
      break;
    default:
      std::cerr << __func__ << ": unexpected row\n";
      exit(1); 

  }
  // skip phase hacking
  // return std::complex<double>(0.,1.)*ret; 
  return ret; 
}



itpp::Mat<std::complex<double> > epsz3(void)
{
  itpp::Mat<std::complex<double> > ret; 
  for(int i = 1; i <= 3; ++i)
    ret.append_row(epsz(i)); 
  return ret; 
}





//  ++ +0 +-  
//  0+ 00 0-
//  -+ -0 --
itpp::Mat<std::complex<double> > wig_D(const ADATXML::Array<int> &mom)
{

  itpp::Mat<std::complex<double> > ret(3,3);
  ret.zeros(); 

  // angles to rotate
  Hadron::CubicCanonicalRotation_t rot;

  if((mom[0] == 0) && (mom[1] == 0) && (mom[2] == 0))
  {
    rot.alpha = 0.;
    rot.beta = 0.;
    rot.gamma = 0.;
  }
  else
  {
    rot = Hadron::cubicCanonicalRotation(mom); 
  }


  for(int m = 0; m < 3; ++m)
    for(int mm = 0; mm < 3; ++mm)
    {
      int m2 = -2*m + 2; 
      int mm2 = -2*mm + 2; 

      ret(m,mm) = SEMBLE::toScalar(Hadron::Wigner_D(2,m2,mm2,rot.alpha,rot.beta,rot.gamma)); 

    }

  return ret; 
}


// active z-y-z rotations
  Tensor<double,2> rotMat(const XMLArray::Array<int> &mom)
  {


    Hadron::CubicCanonicalRotation_t eulerangles = Hadron::cubicCanonicalRotation(mom);

      if((mom[0] == 0) && (mom[1] == 0) && (mom[2] == 0))
  {
    eulerangles.alpha = 0.;
    eulerangles.beta = 0.;
    eulerangles.gamma = 0.;
  }

    Tensor<double,2> A((TensorShape<2>())[4][4],0.),B((TensorShape<2>())[4][4],0.),C((TensorShape<2>())[4][4],0.);


    double a,b,g;
    a = eulerangles.alpha;
    b = eulerangles.beta;
    g = eulerangles.gamma;


    A[0][0] = 1.;
    B[0][0] = 1.;
    C[0][0] = 1.;


    A[1][1] = cos(a)   ;    A[1][2] = sin(a)    ;    A[1][3] =  0.      ; 
    A[2][1] = -sin(a)  ;    A[2][2] = cos(a)    ;    A[2][3] =  0.      ; 
    A[3][1] = 0.       ;    A[3][2] = 0.        ;    A[3][3] =  1.      ; 


    
    B[1][1] = cos(b)   ;    B[1][2] = 0.        ;    B[1][3] = sin(b)   ; 
    B[2][1] = 0.       ;    B[2][2] = 1.        ;    B[2][3] =  0.      ; 
    B[3][1] = -sin(b)  ;    B[3][2] = 0.        ;    B[3][3] = cos(b)   ; 



    C[1][1] = cos(g)   ;    C[1][2] = sin(g)    ;    C[1][3] =  0.      ; 
    C[2][1] = -sin(g)  ;    C[2][2] = cos(g)    ;    C[2][3] =  0.      ; 
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

    return A*B*C;
  }


  Tensor<double,2> rotMat3D(const XMLArray::Array<int> &mom)
  {
    Tensor<double,2> Four = genRotationMatrix(mom);
    Tensor<double,2> Three((TensorShape<2>())[3][3],0.);
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j)
        Three[i][j] = Four[i+1][j+1];

    Three.lower_index(1); 
    
    return Three;
  }


  //! copy a rank 2 tensor into an itpp mat 
  //  -- don't know why I bothered writing this but it may prove to be useful
  template<typename T>
    itpp::Mat<T> toItpp(const Tensor<T,2> &m)
    {
      idx_t row = m.getDim(0);
      idx_t col = m.getDim(1);
      itpp::Mat<T> foo(row,col);

      for(idx_t i = 0; i < row; i++)
        for(idx_t j = 0; j < col; j++)
          foo(i,j) = m[i][j];

      return foo;
    }


itpp::Mat<std::complex<double> >  R3D(const XMLArray::Array<int> &mom)
{
  Tensor<std::complex<double> , 2> foo;
 foo = convertTensorUnderlyingType<std::complex<double>,double,2>(rotMat3D(mom));
 
  return toItpp(foo); 
}




  int 
main(int argc , char * argv[])
{

  if(argc != 4)
  {
    std::cerr << "usage: <px> <py> <pz>" << std::endl;
    exit(1); 
  }


  int px,py,pz; 

  {std::stringstream val(argv[1]); val >> px;}
  {std::stringstream val(argv[2]); val >> py;}
  {std::stringstream val(argv[3]); val >> pz;}


  ADATXML::Array<int> tmom;
  tmom.resize(3);
  tmom[0] = px;
  tmom[1] = py;
  tmom[2] = pz; 


  itpp::Mat<std::complex<double> > D, R, destD, destR, teps = epsz3(); 
  D = wig_D(tmom); 
  R = R3D(tmom); 
  destD = D * teps; 
  destR = itpp::transpose(R * itpp::transpose(teps)); 



  std::cout << "epsz = \n" << print(teps) << std::endl;

  std::cout << "D = \n" << print(D) << std::endl;
  
  std::cout << "R = \n" << print(R) << std::endl;

  std::cout << "destD = \n" << print(destD) << std::endl;  

  std::cout << "destR = \n" << print(destR) << std::endl;  

  std::cout << "destR - destD = \n" << print(destD - destR) << std::endl;


  return 0; 
}
