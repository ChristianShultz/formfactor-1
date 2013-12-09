#ifndef LORENTZFF_POLARIZATION_EMBEDDING_H
#define LORENTZFF_POLARIZATION_EMBEDDING_H

#include <complex>

#include "radmat/utils/polarisation_tensors.h"
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
      typedef typename genPolTens<J>::mom_t mom_t; // Array<int> 

      HelicityPolarizationTensor(void) {}
      virtual ~HelicityPolarizationTensor(void) {}

      // this took me awhile to find, one needs to be careful about how to round
      // this take the fabs to get magnitude and then attach phase 
      inline int round(const double &d) {return d < 0 ? -1*int(fabs(d) + 0.5) : int(fabs(d) + 0.5);}

      mom_t int_based_mom(const Tensor<double,1> &p, const double mom_factor)
      {
        mom_t ret;
        ret.resize(3);
        ret[0] = round(p[1]/mom_factor);
        ret[1] = round(p[2]/mom_factor);
        ret[2] = round(p[3]/mom_factor);

        //  std::cout << __func__ << ": passed in " << p 
        //    << "\n passing out " << ret[0] << ret[1] << ret[2] << std::endl;

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
        operator()(const Tensor<double,1> &p, 
            const int hel, 
            const double mom_factor)
        {
          mom_t mom = int_based_mom(p,mom_factor);
          genPolTens<J> foo(mom);
          Tensor<std::complex<double> , J> tens = foo(p[0],hel,mom_factor); 
          typename Tensor<std::complex<double> , J>::iterator it;

          for(it = tens.begin(); it != tens.end(); ++it)
            *it = round_to_zero(*it); 

          return tens; 
        }
    };


  //  embed the target helicity as a template parameter
  //      and set apply the outside world rest rotation 
  //      convention here 
  //
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  template<idx_t J, int hel>
    struct embedHelicityPolarizationTensor
    {
      typedef typename genPolTens<J>::mom_t mom_t; // Array<int> 
      virtual ~embedHelicityPolarizationTensor() {}

      // return the polarization tensor associated with p,
      //    if the particle with momentum p is at rest then 
      //    apply the rotation convention for particle p_prime 
      //    to the rest polarization vector to force a rotation
      //
      //    this corresponds to taking the z frame of the 
      //    particle in flight as the reference frame
      virtual Tensor<std::complex<double> , J>
        ptensor(const Tensor<double,1> &p,
            const double mom_factor, 
            const Tensor<double,1> p_prime) const
        {
          HelicityPolarizationTensor<J> foo;
          mom_t p_mom = foo.int_based_mom(p,mom_factor);  
          mom_t p_prime_mom = foo.int_based_mom(p_prime,mom_factor);  
          Tensor<std::complex<double>,J> epsilon(foo(p,hel,mom_factor)); 

          if ( is_rest(p_mom) ) 
            return apply_rest_rotation_convention( epsilon , p_prime_mom ); 

          return epsilon;
        }


      virtual Tensor<std::complex<double> , J> 
        apply_rest_rotation_convention(
            const Tensor<std::complex<double> ,J> &eps,
            const mom_t &p_prime) const
        {
          // redstar allows us to twist one end w/o twisting the other 
          //   so the idea of a rotation is a bit more complicated
          //
          //  when we spin about the momentum direction of a particle
          //  redstar does not spin the definition of the corresponding 
          //  matrix element, thus we are effectively measuring 
          //  the matrix elems in different frames.  I think of this 
          //  as twisting each end of the matrix element independently 
          //
          //  the result then is that the phases that occurr between 
          //  unallowed rotations (spinning about a momentum direction)
          //  need to be accounted for in the decompositions
          //
          //  if in fact the form factor is just a number then all possible
          //  complex phase information must be stored somehow in the 
          //  polarization vectors and getting these correct will allow
          //  us to use the unrelated measurements in the same llsq
          //

          return eps;  
        }


      bool is_rest(const mom_t &mom) const
      {
        return ((mom[0] == 0) && (mom[1] == 0) && (mom[2] == 0)); 
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
        left_p_tensor(const Tensor<double,1> &p, 
            const double mom_factor ,
            const Tensor<double,1> &pp) const
        {
          embedHelicityPolarizationTensor<J_left,hel_left> foo; 
          return foo.conjugate( foo.ptensor(p,mom_factor,pp) ); 
        }
    };


  ////////////////////////////////////////////////////////////////////////
  //  right polarization tensors (initial state, creation ops)
  template<idx_t J_right, int hel_right>
    struct rightPTensor
    {
      virtual ~rightPTensor() {}
      virtual Tensor<std::complex<double> , J_right > 
        right_p_tensor(const Tensor<double,1> &p, 
            const double mom_factor,
            const Tensor<double,1> &pp ) const
        {
          embedHelicityPolarizationTensor<J_right,hel_right> foo; 
          return foo.ptensor(p,mom_factor,pp); 
        }
    };

}





#endif /* LORENTZFF_POLARIZATION_EMBEDDING_H */
