#ifndef LORENTZFF_ROTATED_LEVI_CIVITA_TENSOR_H
#define LORENTZFF_ROTATED_LEVI_CIVITA_TENSOR_H 

#include "lorentzff_formfac_utils.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/tensor.h"
#include "adat/singleton.h"

namespace radmat
{

  namespace RotatedLeviCivitaTensorEnv
  {
    rHandle<Tensor<double,4> >
      call_factory(const Tensor<double,1> &left,
          const Tensor<double,1> &right,
          const double mom_kick);

    rHandle<Tensor<double,4> >
      call_factory(const mom_t &l, const mom_t &r); 

    typedef Util::SingletonHolder<std::map<std::string,Tensor<double,4> > >
      TheRotatedLeviCivitaTensorGenerator; 


    bool registerAll(void); 

  } // RotatedLeviCivitaTensorEnv


} // radmat



#endif /* LORENTZFF_ROTATED_LEVI_CIVITA_TENSOR_H */
