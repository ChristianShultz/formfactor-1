#ifndef LORENTZFF_POLARIZATION_EMBEDDING_H
#define LORENTZFF_POLARIZATION_EMBEDDING_H

#include <complex>

#include "lorentzff_polarization_tensors.h"
#include "lorentzff_canonical_rotations.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "hadron/clebsch.h"
#include "hadron/irrep_util.h"


namespace radmat
{


  // provide an interface to churn out helicity polarization 
  // tensors from the momentum tensors 
  template<idx_t J>
    struct HelicityPolarizationTensor
    {

      HelicityPolarizationTensor(void) {}
      virtual ~HelicityPolarizationTensor(void) {}

      inline int 
        round(const double &d)
        {
          return d < 0 ? -1*int(fabs(d) + 0.5) : int(fabs(d) + 0.5);
        }

      mom_t 
        int_based_mom(const Tensor<double,1> &p, const double mom_factor)
      {
        mom_t ret;
        ret.resize(3);
        ret[0] = round(p[1]/mom_factor);
        ret[1] = round(p[2]/mom_factor);
        ret[2] = round(p[3]/mom_factor);

        return ret; 
      };

      double 
        round_to_zero(const double &d, const double thresh) const
        {
          if ( fabs(d) <  thresh ) 
            return 0.;
          else
            return d; 
        }

      std::complex<double> 
        round_to_zero(const std::complex<double> &d) const
        {
          double thresh = 1e-6; 
          return std::complex<double>(round_to_zero(d.real(),thresh) , round_to_zero(d.imag(),thresh)); 
        }


      // do work 
      Tensor<std::complex<double>, J> 
        z_axis_helicity_tensor(const Tensor<double,1> &t,   // target 
            const int hel, 
            const double mom_factor)
        {
          double mod_p = sqrt( t[1]*t[1] + t[2]*t[2] + t[3]*t[3] );
          ZAxisHelicityTensor<J> foo;

          Tensor<std::complex<double> , J> tens = foo(mod_p,hel,t[0]); 

          typename Tensor<std::complex<double> , J>::iterator it;
          for(it = tens.begin(); it != tens.end(); ++it)
            *it = round_to_zero(*it); 

          return tens; 
        }
    };


  ////////////////////////////////////////////////////////////////////////
  //
  //  embed target helicity
  //
  ////////////////////////////////////////////////////////////////////////
  template<idx_t J, int hel>
    struct embedHelicityPolarizationTensor
    : public HelicityPolarizationTensor<J> 
    {

      virtual ~embedHelicityPolarizationTensor() {}

      virtual Tensor<std::complex<double>, J>
        rotate(const Tensor<std::complex<double>,J> &z_axis,
            const rHandle<Tensor<double,2> > &R)
        {
          Tensor<std::complex<double>,J> rot(z_axis); 
          Tensor<std::complex<double>,2> Rloc;
          Rloc = convertTensorUnderlyingType<std::complex<double>,double,2>(*R); 

          // rotate all indicies 
          for(int i = 0; i < J; ++i)
            rot = contract(Rloc,rot,1,i); 
         
          // clean up 
          typename Tensor<std::complex<double> , J>::iterator it;
          for(it = rot.begin(); it != rot.end(); ++it)
            *it = round_to_zero(*it); 

          return rot; 
        }


      virtual Tensor<std::complex<double> , J> 
        conjugate(const Tensor<std::complex<double> , J> &inp) const
        {
          Tensor<std::complex<double> , J> foo = inp;
          typename Tensor<std::complex<double> , J>::iterator it;

          for(it = foo.begin(); it != foo.end(); ++it)
            *it = std::conj(*it); 

          return foo;
        }
    };



  //  Classes that need polarization tensors should derive from the 
  //  leftPTensor and rightPTensor class to enforce consistency 
  //
  //  take care of the left/right who is complex convention
  //
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  //  left polarization tensors (final state, annih ops)
  template<idx_t J_left, int hel_left>
    struct leftPTensor
    {

      virtual ~leftPTensor() {}

      virtual Tensor<std::complex<double> , J_left > 
        left_p_tensor(const Tensor<double,1> &lefty, 
            const Tensor<double,1> &righty, 
            const double mom_factor) const
        {
          embedHelicityPolarizationTensor<J_left,hel_left> foo; 
          Tensor<std::complex<double> , J_left> eps_z = foo.z_axis_helicity_tensor(lefty,hel_left,mom_factor); 
          mom_t l = foo.int_based_mom(lefty,mom_factor);
          mom_t r = foo.int_based_mom(righty,mom_factor); 
        
          rHandle<Tensor<double,2> > R = radmat::LatticeRotationEnv::get_left_rotation(l,r);      
          Tensor<std::complex<double>, J_left> eps = foo.rotate(eps_z,R); 

          return foo.conjugate( eps ); 
        }
    };


  ////////////////////////////////////////////////////////////////////////
  //  right polarization tensors (initial state, creation ops)
  template<idx_t J_right, int hel_right>
    struct rightPTensor
    {

      virtual ~rightPTensor() {}

      virtual Tensor<std::complex<double> , J_right > 
        right_p_tensor(const Tensor<double,1> &lefty, 
            const Tensor<double,1> &righty, 
            const double mom_factor) const
        {
          embedHelicityPolarizationTensor<J_right,hel_right> foo; 
          Tensor<std::complex<double> , J_right> eps_z = foo.z_axis_helicity_tensor(righty,hel_right,mom_factor); 
          mom_t l = foo.int_based_mom(lefty,mom_factor);
          mom_t r = foo.int_based_mom(righty,mom_factor); 
        
          rHandle<Tensor<double,2> > R = radmat::LatticeRotationEnv::get_right_rotation(l,r);      

          return foo.rotate(eps_z,R); 
        }
    };

}



#endif /* LORENTZFF_POLARIZATION_EMBEDDING_H */
