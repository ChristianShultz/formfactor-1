#ifndef FAKE_3PT_FUNCTION_H_H_GUARD
#define FAKE_3PT_FUNCTION_H_H_GUARD

#include "ensem/ensem.h"
#include "fake_3pt_function_aux.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include <string>

namespace radmat
{

  struct Fake3ptKey
  {
    int lorentz;
    int hel_source;
    int hel_sink;
    std::string elemIDBase;
    ENSEM::EnsemReal E_source;
    ENSEM::EnsemReal E_sink;
    pProp_t mom;
  };



  template<typename T>
    struct Fake3ptCorr
    {

      Fake3ptCorr(const typename ADAT::Handle<FakeDataInputs<T> > &inputs, 
          const int lorentz, 
          const int hel_source, 
          const int hel_sink, 
          const pProp_t &mom)
        : m_inputs(inputs)
      { 
        m_key.lorentz = lorentz;
        m_key.hel_source = hel_source;
        m_key.hel_sink = hel_sink;
        m_key.mom = mom;

        if(inputs->working->ini.matElemProps.left_target > inputs->working->ini.matElemProps.right_target)  
          m_key.elemIDBase = inputs->working->ini.matElemProps.upper;
        else if(inputs->working->ini.matElemProps.left_target < inputs->working->ini.matElemProps.right_target)
          m_key.elemIDBase = inputs->working->ini.matElemProps.lower;
        else
          m_key.elemIDBase = inputs->working->ini.matElemProps.diag;


        m_key.E_source = getE(inputs->working->specsource,inputs->working->ini.matElemProps.right_target);
        m_key.E_sink = getE(inputs->working->specsink,inputs->working->ini.matElemProps.left_target);


        m_threePtFcn.resize(m_key.E_source.size());
        m_threePtFcn.resizeObs(inputs->working->specsink.size());
        m_threePtFcn = SEMBLE::toScalar(T(0.));

        for(unsigned int t = 0; t < inputs->working->specsink.size(); t++)
          ENSEM::pokeObs(m_threePtFcn, makeFakeDataPoint(inputs,mom,hel_sink,hel_source,t,lorentz) ,t);

      }

      public:
      Fake3ptKey& getKey(void) const {return m_key;}
      Fake3ptKey& getKey(void) {return m_key;}
      typename SEMBLE::PromoteEnsemVec<T>::Type& get3pt(void) const {return m_threePtFcn;}
      typename SEMBLE::PromoteEnsemVec<T>::Type& get3pt(void) {return m_threePtFcn;}
      typename ADAT::Handle<FakeDataInputs<T> >& getInputs(void) const {return m_inputs;}
      typename ADAT::Handle<FakeDataInputs<T> >& getInputs(void) {return m_inputs;}



      private:  
      ENSEM::EnsemReal getE(const std::vector<SEMBLE::SembleVector<double> > &spec, const int target)
      {
        ENSEM::EnsemReal ret; 
        ret.resize(spec[0].getB());
        std::vector<SEMBLE::SembleVector<double> >::const_iterator it;
        for(it = spec.begin(); it != spec.end(); it++)
          ret += (*it)(target);

        ret /= SEMBLE::toScalar(double(spec.size()));

        return ret;      
      }

      private:
      Fake3ptCorr(void);
      Fake3ptCorr(const Fake3ptCorr<T> &o);

      Fake3ptKey m_key;
      typename SEMBLE::PromoteEnsemVec<T>::Type m_threePtFcn;
      typename ADAT::Handle<FakeDataInputs<T> > m_inputs;
    };



}

#endif
