#ifndef LORENTZFF_WIGNER_D_MATRIX_FACTORY_H
#define LORENTZFF_WIGNER_D_MATRIX_FACTORY_H 

#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "io/adat_xmlio.h"
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include <complex>
#include <map>
#include <string>

namespace radmat
{
  typedef Tensor<std::complex<double>,2> WignerMatrix_t; 


  namespace WignerDMatrixEnv
  {
    typedef Util::SingletonHolder<
      std::map<std::string,WignerMatrix_t> >
      TheWignerDMatrixFactory;

    bool registerAll(const int Jmax); 

    // the most positive index is mapped to zero, next most is 1..
    WignerMatrix_t* call_factory(const mom_t &, const int); 

    template<int J> 
      WignerMatrix_t* call_factory(const mom_t &p)
      {
        return call_factory(p,J); 
      } 

  } // WignerDMatrixEnv 


} // radmat



#endif /* LORENTZFF_WIGNER_D_MATRIX_FACTORY_H */
