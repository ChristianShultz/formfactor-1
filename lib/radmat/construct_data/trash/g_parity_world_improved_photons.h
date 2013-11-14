#ifndef G_PARITY_WORLD_IMPROVED_PHOTONS_H
#define G_PARITY_WORLD_IMPROVED_PHOTONS_H 

#include <complex>
#include "io/adat_xmlio.h"

namespace radmat
{

  namespace Rho2ImprovementTerm
  {
    struct hWeight_t
    {
      hWeight_t(const int hh, const std::complex<double> &ww)
        : h(hh) , w(ww) 
      { }

      int h;
      std::complex<double> w; 
    };



    std::vector<hWeight_t> 
      spatial_rho2_insertion(const int h, 
          const ADATXML::Array<int> mom, 
          const double m_f, 
          const double m_i,
          const double mom_unit);

    temporal_rho2_insertion(const int h, 
        const ADATXML::Array<int> mom, 
        const double m_f, 
        const double m_i,
        const double mom_unit);
  }
}





#endif /* G_PARITY_WORLD_IMPROVED_PHOTONS_H */
