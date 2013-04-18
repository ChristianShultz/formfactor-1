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

      inline int round(const double &d) {return int(d + 0.5);}

      mom_t int_based_mom(const Tensor<double,1> &p, const double mom_factor)
      {
        mom_t ret;
        ret.resize(3);
        ret[0] = round(p[1]/mom_factor);
        ret[1] = round(p[2]/mom_factor);
        ret[2] = round(p[3]/mom_factor);
        return ret; 
      };

      Tensor<std::complex<double>, J> operator()(const Tensor<double,1> &p, const int hel, const double mom_factor)
      {
        mom_t mom = int_based_mom(p,mom_factor);
        genPolTens<J> foo(mom);
        return foo.get(p[0],hel,mom_factor); 
      }
    };

  template<idx_t J, int hel>
    struct embedHelicityPolarizationTensor
    {
      virtual Tensor<std::complex<double> , J> ptensor(const Tensor<double,1> &p, const double mom_factor) const
      {
        HelicityPolarizationTensor<J> foo;
        return foo(p,hel,mom_factor); 
      }

      virtual Tensor<std::complex<double> , J> conjugate(const Tensor<std::complex<double> , J> &inp) const
      {
        Tensor<std::complex<double> , J> foo = inp;
        typename Tensor<std::complex<double> , J>::iterator it;

        for(it = foo.begin(); it != foo.end(); ++it)
          *it = std::conj(*it); 

        return foo;
      }

    };

}





#endif /* LORENTZFF_POLARIZATION_EMBEDDING_H */
