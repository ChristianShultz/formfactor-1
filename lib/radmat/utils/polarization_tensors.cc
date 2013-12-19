/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 04-05-2013

 * Last Modified : Mon 09 Dec 2013 03:08:45 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "polarization_tensors.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"

#include "semble/semble_meta.h"

#include "tensor.h"

namespace radmat
{
  namespace
  {
    // match to the Rho_J1_Op in adat, phase is -complex_one * complex_i
    std::complex<double> phase_hack(std::complex<double>(0.,-1.)); 


    // 1 based
    itpp::Vec<std::complex<double> > epsz3(int m_row)
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
          std::cerr << __func__ << ": unexpected row" << m_row << "\n";
          exit(1); 
      }

      // match to the Rho_J1_Op in adat, phase is -complex_one * complex_i
      return phase_hack * ret;  
    }


    // this is 1 based
    itpp::Vec<std::complex<double> > rotated_eps3(const ADATXML::Array<int> &mom, 
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
        return epsz3(lambda_row); 
      }
      else
      {
        rot = Hadron::cubicCanonicalRotation(mom); 
      }

      for(int i = 1; i <= 3; ++i)
      {
        int M2 = -2*(i - 1) +J2;
        double phase = 1.;

        ENSEM::Complex D = Hadron::Wigner_D(J2,M2,lambda2,rot.alpha,rot.beta,rot.gamma); 
        if(std::norm(SEMBLE::toScalar(D)) > 0.00001) 
          eps = eps + phase*SEMBLE::toScalar(D) * epsz3(i); 
      }

      return eps;

    }


    // 1 based
    itpp::Vec<std::complex<double> > epsz4(int m_row, const double p, const double E)
    {
      itpp::Vec<std::complex<double> > eps3 = epsz3(m_row); 
      itpp::Vec<std::complex<double> > ret(4);
      ret.zeros(); 
      double m = sqrt(E*E - p*p);

      POW2_ASSERT(m > 1e-14); 

      double t = p/m;
      double z = E/m; 

      switch(m_row)
      {

        case 3:
          ret[1] = eps3[0];
          ret[2] = eps3[1];
          break;
        case 2:
          ret[0] = t * eps3[2];
          ret[3] = z * eps3[2];
          break; 
        case 1:
          ret[1] = eps3[0];
          ret[2] = eps3[1];
          break;
        default:
          std::cerr << __func__ << ": unexpected row\n";
          exit(1); 

      }

      // undo the phase hack from the redstar op eps3
      return std::conj(phase_hack) * ret; 
    }



    itpp::Vec<std::complex<double> > rotated_eps4(const ADATXML::Array<int> &mom, 
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
        return epsz4(lambda_row,p,E); 
      }
      else
      {
        rot = Hadron::cubicCanonicalRotation(mom); 
      }

      Tensor<std::complex<double>, 2> RotMat;
      RotMat = convertTensorUnderlyingType<std::complex<double> , double , 2 > (genRotationMatrix(mom));
      itpp::Mat<std::complex<double> > R(4,4); 
      R.zeros();

      for(int mu = 0; mu < 4; ++mu)
        for(int nu = 0; nu < 4; ++nu)
          R(mu,nu) = RotMat[mu][nu];

      eps = R*epsz4(lambda_row,p,E);  

      return eps;
    }

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
    itpp::Vec<std::complex<double> > eps(rotated_eps3(mom,lambda_row));

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
    itpp::Vec<std::complex<double> > eps(rotated_eps4(mom,lambda_row,p,E)); 

    dest[0] = eps[0];
    dest[1] = eps[1];
    dest[2] = eps[2];
    dest[3] = eps[3];

    return dest; 
  }

  Tensor<std::complex<double> , 1> creation_op_J1_3z(const int lambda)
  {
    int lambda_row = remapHelicity_1based(lambda,1); 
    Tensor<std::complex<double> , 1> dest((TensorShape<1>())[3],0.); 
    itpp::Vec<std::complex<double> > eps(epsz3(lambda_row));

    dest[0] = eps[0];
    dest[1] = eps[1];
    dest[2] = eps[2];

    return dest; 
  }

  Tensor<std::complex<double> , 1> creation_op_J1_4z(const double mod_p, 
      const int lambda,
      const double E) 
  {
    int lambda_row = remapHelicity_1based(lambda,1); 
    Tensor<std::complex<double> , 1> dest((TensorShape<1>())[4],0.);

    itpp::Vec<std::complex<double> > eps(epsz4(lambda_row,mod_p,E)); 

    dest[0] = eps[0];
    dest[1] = eps[1];
    dest[2] = eps[2];
    dest[3] = eps[3];

    return dest; 
  }


} // radmat