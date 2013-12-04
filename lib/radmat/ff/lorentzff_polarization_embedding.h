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



  // 
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

        //        std::cout << __func__ << " sent in " << p << std::endl; 
        mom_t ret;
        ret.resize(3);
        ret[0] = round(p[1]/mom_factor);
        ret[1] = round(p[2]/mom_factor);
        ret[2] = round(p[3]/mom_factor);

        //        std::cout << __func__ << " sending out " << ret[0] << ret[1] << ret[2] << std::endl;

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

  template<idx_t J, int hel>
    struct embedHelicityPolarizationTensor
    {
      virtual ~embedHelicityPolarizationTensor() {}

      virtual Tensor<std::complex<double> , J>
        ptensor(const Tensor<double,1> &p,
            const double mom_factor) const
        {
          HelicityPolarizationTensor<J> foo;
          return foo(p,hel,mom_factor); 
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


  template<idx_t J_left, int hel_left>
    struct leftPTensor
    {
      virtual ~leftPTensor() {}

      virtual Tensor<std::complex<double> , J_left > 
        left_p_tensor(const Tensor<double,1> &p, const double mom_factor ) const
        {
          embedHelicityPolarizationTensor<J_left,hel_left> foo; 
          return foo.conjugate( foo.ptensor(p,mom_factor) ); 
        }
    };


  template<idx_t J_right, int hel_right>
    struct rightPTensor
    {
      virtual ~rightPTensor() {}
      virtual Tensor<std::complex<double> , J_right > 
        right_p_tensor(const Tensor<double,1> &p, const double mom_factor ) const
        {
          embedHelicityPolarizationTensor<J_right,hel_right> foo; 
          return foo.ptensor(p,mom_factor); 
        }
    };


}





#endif /* LORENTZFF_POLARIZATION_EMBEDDING_H */
