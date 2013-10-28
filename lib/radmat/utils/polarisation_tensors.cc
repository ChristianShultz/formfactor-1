/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 04-05-2013

 * Last Modified : Mon 21 Oct 2013 10:51:50 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "polarisation_tensors.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"

#include "semble/semble_meta.h"

#include "tensor.h"

namespace radmat
{


  // so the case where the box rotates but the photon is at rest kind of sucks..



  namespace
  {


    //
    //  3D for solving photon helicities -> cartesian, this is simply a definition from adat
    //


    // all phases and thus all definitions ultimately result from this function and the 
    // one that multiplies by the d-matrix

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
      // ret *= std::complex<double>(0.,1.);  
       return ret; 
    }


    // this is 1 based
    itpp::Vec<std::complex<double> > eps_lambda(const ADATXML::Array<int> &mom, 
        const int lambda_row)
    {
      itpp::Vec<std::complex<double> > eps(3); 
      eps.zeros();

      int J2 = 2; 
      int lambda2 = -2*(lambda_row - 1) + J2; 

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
      Tensor<std::complex<double>, 2> RotMat;
      RotMat = convertTensorUnderlyingType<std::complex<double> , double , 2 > (genRotationMatrix(mom));
      itpp::Mat<std::complex<double> > R(3,3); 
      R.zeros();

      for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
          R(i,j) = RotMat[i+1][j+1];


#if 0
#warning " USING R *epsz(lambda) for 3D polarization "

      eps = R*epsz(lambda_row); 
 
#else

#warning " USING  D_{m,lambda} * epsz(m) for 3D polarization "

      for(int i = 1; i <= 3; ++i)
      {

        int M2 = -2*(i - 1) +J2;
        double phase = 1.;

        // attempted phase cooking
        // int M2 = 2*(i -1) - J2; // m -> -m
        // double phase = (i %2 == 0) ? 1. : -1.; 

        ENSEM::Complex D = Hadron::Wigner_D(J2,M2,lambda2,rot.alpha,rot.beta,rot.gamma); 
        if(std::norm(SEMBLE::toScalar(D)) > 0.00001) 
          eps = eps + phase*SEMBLE::toScalar(D) * epsz(i); 
      }
#endif 

      return eps;

    }




    //
    //  4D - this encodes the transformation properties of a polarization vector
    //





    itpp::Vec<std::complex<double> > epsz(int m_row, const double p, const double E)
    {
      // itpp::Vec<std::complex<double> > eps3 = epsz(m_row); 
      itpp::Vec<std::complex<double> > ret(4);
      ret.zeros(); 
      double m = sqrt(E*E - p*p);

      POW2_ASSERT(m > 1e-14); 

      double t = p/m;
      double z = E/m; 

      switch(m_row)
      {

        case 3:
          ret[1] = std::complex<double>(1./sqrt(2.),0.); 
          ret[2] = std::complex<double>(0.,-1./sqrt(2.));
          // ret[1] = eps3[0];
          // ret[2] = eps3[1];
          break;
        case 2:
          ret[0] = std::complex<double>(1.,0.)*t; 
          ret[3] = std::complex<double>(1.,0.)*z; 
          // ret[0] = t * eps3[2];
          // ret[3] = z * eps3[2];
          break; 
        case 1:
          ret[1] = std::complex<double>(-1./sqrt(2.),0.); 
          ret[2] = std::complex<double>(0.,-1./sqrt(2.));          
          // ret[1] = eps3[0];
          // ret[2] = eps3[1];
          break;
        default:
          std::cerr << __func__ << ": unexpected row\n";
          exit(1); 

      }

      // this phase is the guy for spatial insertions, for temporal guys it introduces a factor of i 
      // return std::complex<double>(0.,1.)*ret; 

       return ret; 
    }




    itpp::Vec<std::complex<double> > eps_lambda(const ADATXML::Array<int> &mom, 
        const int lambda_row, 
        const double p,
        const double E)
    {
      itpp::Vec<std::complex<double> > eps(4); 
      eps.zeros();

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

#if 1

#warning "using R^u_v eps^v(p,lam) = eps^u(Rp,lam)"

      Tensor<std::complex<double>, 2> RotMat;
      RotMat = convertTensorUnderlyingType<std::complex<double> , double , 2 > (genRotationMatrix(mom));
      itpp::Mat<std::complex<double> > R(4,4); 
      R.zeros();



      for(int mu = 0; mu < 4; ++mu)
        for(int nu = 0; nu < 4; ++nu)
          R(mu,nu) = RotMat[mu][nu];


      eps = R*epsz(lambda_row,p,E);  


    /*
      std::cout << __func__ << " sent in mom = " << mom[0] << mom[1] << mom[2] << std::endl;
      std::cout << __func__ << " euler angles = " << rot.alpha << " " << rot.beta << " " << rot.gamma << std::endl;
      std::cout << __func__ << " rot mat \n" << R << std::endl; 
      std::cout << __func__ << " eps z " << epsz(lambda_row,p,E) << std::endl;
      std::cout << __func__ << " R*epsz " << eps << std::endl;
      */

#else 
#warning "using eps^u(Rp,lam) = D_m,lam * R * eps(pz,m)"

      for(int i = 1; i <= 3; ++i)
      {
        int M2 = -2*(i - 1) +J2;
        double phase = 1.; 

        ENSEM::Complex D = Hadron::Wigner_D(J2,M2,lambda2,rot.alpha,rot.beta,rot.gamma);
        if(std::norm(SEMBLE::toScalar(D)) > 0.00001) 
          eps = eps + phase*SEMBLE::toScalar(D) * R * epsz(i,p,E); 
      }
      
       return itpp::conj(eps);
#endif 

      return eps;
    }


#if 0 

    // form the matrix eps
    //
    // |j_+|    | eps_x eps_y eps_z |    |j_x|              
    // |j_0|  = | eps_x eps_y eps_z |  x |j_y|               
    // |j_-|    | eps_x eps_y eps_z |    |j_z|     
    //

    itpp::Mat<std::complex<double> > eps3d(const ADATXML::Array<int> &mom, const bool create)
    {
      itpp::Mat<std::complex<double> > ret(3,3);
      ret.zeros();   

      for(int i = 0; i < 3; ++i)
      {
        ret.set_row(i,eps_lambda(mom,i+1)); 
      }
      return ret; 
    }

#endif 

  } // anonomyous 



  int remapHelicity_1based(const int H, const int J)
  {
    POW2_ASSERT( abs(H) <= J ); 
    return J - H + 1;
  }



  Tensor<std::complex<double> , 1> creation_op_J1_3(const ADATXML::Array<int> &mom, 
      const int lambda)
  {
    int lambda_row = remapHelicity_1based(lambda,1); 
    Tensor<std::complex<double> , 1> dest((TensorShape<1>())[3],0.); 
    itpp::Vec<std::complex<double> > eps(eps_lambda(mom,lambda_row));

    dest[0] = eps[0];
    dest[1] = eps[1];
    dest[2] = eps[2];

    return dest; 
  }

  Tensor<std::complex<double> , 1> creation_op_J1_4(const ADATXML::Array<int> &mom, 
      const int lambda,
      const double E, 
      const double mom_fac) 
  {
    int lambda_row = remapHelicity_1based(lambda,1); 
    double p(0.);
    Tensor<std::complex<double> , 1> dest((TensorShape<1>())[4],0.);

    for(int i = 0; i < 3; ++i)
      p += double(mom[i]*mom[i]); 

    p = sqrt(mom_fac*mom_fac*p);  
    itpp::Vec<std::complex<double> > eps(eps_lambda(mom,lambda_row,p,E)); 

    dest[0] = eps[0];
    dest[1] = eps[1];
    dest[2] = eps[2];
    dest[3] = eps[3];

    return dest; 
  }


}
