#ifndef THREE_POINT_H_H_GUARD
#define THREE_POINT_H_H_GUARD

#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/fake_data/fake_3pt_function.h"



namespace radmat
{


  template<typename T>
    struct ThreePointCorrelator
    {
      int lorentz;
      int hel_source;
      int hel_sink;
      int t_source;
      int t_sink;
      std::string elemIDBase;

      ENSEM::EnsemReal E_source;
      ENSEM::EnsemReal E_sink;
      typename SEMBLE::PromoteEnsem<T>::Type Z_source;
      typename SEMBLE::PromoteEnsem<T>::Type Z_sink;
      ENSEM::EnsemReal Q2;
      typename SEMBLE::PromoteEnsemVec<T>::Type C3pt;
      double mom_factor;
      pProp_t mom;
    };


  template<typename T>  
    ThreePointCorrelator<T> makeThreePointFromFake(const Fake3ptCorr<T> &c)
    {
      ThreePointCorrelator<T> C3;
      C3.lorentz = c.getKey().lorentz;
      C3.hel_source = c.getKey().hel_source;
      C3.hel_sink = c.getKey().hel_sink;
      C3.elemIDBase = c.getKey().elemIDBase;
      C3.E_source = c.getKey().E_source;
      C3.E_sink = c.getKey().E_sink;
      C3.t_source = c.getKey().t_source;
      C3.t_sink = c.getKey().t_sink;

      // NB !!
      // normalization factor included here, the 1/2E from the lorentz invariant metric
      C3.Z_source = c.getKey().Z_source/c.getKey().E_source/SEMBLE::toScalar(double(2));
      C3.Z_sink = c.getKey().Z_sink/c.getKey().E_sink/SEMBLE::toScalar(double(2));

      /*
         std::cout << __func__ << std::endl;
         std::cout << SEMBLE::toScalar(ENSEM::mean(c.getKey().E_source)) << std::endl;
         std::cout << SEMBLE::toScalar(ENSEM::mean(c.getKey().E_sink)) << std::endl;
         std::cout << SEMBLE::toScalar(ENSEM::mean(c.getKey().Z_source)) << std::endl;
       */

      C3.Q2 = c.getKey().Q2;
      C3.C3pt = c.get3pt();
      C3.mom_factor = c.getKey().mom_factor;
      C3.mom = c.getKey().mom;

      return C3;
    }


} // namespace radmat

#endif
