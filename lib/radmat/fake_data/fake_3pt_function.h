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
    int t_source;
    int t_sink;

    std::string elemIDBase;
    ENSEM::EnsemReal E_source;
    ENSEM::EnsemReal E_sink;
    ENSEM::EnsemComplex Z_source;
    ENSEM::EnsemComplex Z_sink;
    ENSEM::EnsemReal Q2;
    pProp_t mom;

    double mom_factor;
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

        int source = inputs->original->ini.matElemProps.right_target;
        int sink = inputs->original->ini.matElemProps.left_target;

        double Esink,Esource,meanvar;
        Esink = SEMBLE::toScalar(ENSEM::mean(inputs->original->specsink[0](sink)));
        Esource = SEMBLE::toScalar(ENSEM::mean(inputs->original->specsource[0](source)));
        meanvar = inputs->original->ini.stateProps.mProps.sourceVarO + inputs->original->ini.stateProps.mProps.sinkVarO;

        // try to figure out if we are looking at a diagonal guy or not -- duct tape and dreams
        if( fabs(Esink - Esource)/std::max(Esink,Esource) <= meanvar)
          m_key.elemIDBase = inputs->working->ini.matElemProps.diag;
        else  
          m_key.elemIDBase = inputs->working->ini.matElemProps.off;

        // fill in the key
        m_key.lorentz = lorentz;
        m_key.hel_source = hel_source;
        m_key.hel_sink = hel_sink;
        m_key.mom = mom;
        m_key.E_source = getE(inputs->working->specsource,inputs->working->ini.matElemProps.right_target);
        m_key.E_sink = getE(inputs->working->specsink,inputs->working->ini.matElemProps.left_target);
        m_key.Z_source = getZ(inputs->working->zsource[inputs->working->ini.timeProps.tsource],
            inputs->working->ini.matElemProps.right_target);
        m_key.Z_sink = getZ(inputs->working->zsink[inputs->working->ini.timeProps.tsink],
            inputs->working->ini.matElemProps.left_target);

        SemblePInv_t mom_inv =   makeMomInvariants( m_key.E_sink,
            m_key.E_source,
            mom.momSink,
            mom.momSource,
            inputs->working->mom_factor);

        SEMBLE::SembleVector<double> q = (mom_inv.first - mom_inv.second);
        ENSEM::EnsemReal Q2 = (-q(0)*q(0) + q(1)*q(1) + q(2)*q(2) + q(3)*q(3));

        m_key.Q2 = Q2;
        m_key.t_source = inputs->working->ini.timeProps.tsource;
        m_key.t_sink = inputs->working->ini.timeProps.tsink;
        m_key.mom_factor = inputs->working->mom_factor;

        // fake up the three pt function
        m_threePtFcn.resize(m_key.E_source.size());
        m_threePtFcn.resizeObs(inputs->working->specsink.size());
        m_threePtFcn = SEMBLE::toScalar(T(0.));

        for(unsigned int t = 0; t < inputs->working->specsink.size(); t++)
          ENSEM::pokeObs(m_threePtFcn, makeFakeDataPoint(inputs,mom,hel_sink,hel_source,t,lorentz) ,t);

      }


      Fake3ptCorr(const Fake3ptCorr<T> &o)
        :m_key(o.m_key) , m_threePtFcn(o.m_threePtFcn) , m_inputs(o.m_inputs)
      {

      }


      public:
      const Fake3ptKey& getKey(void) const {return m_key;}
      Fake3ptKey& getKey(void) {return m_key;}
      const typename SEMBLE::PromoteEnsemVec<T>::Type& get3pt(void) const {return m_threePtFcn;}
      typename SEMBLE::PromoteEnsemVec<T>::Type& get3pt(void) {return m_threePtFcn;}
      const typename ADAT::Handle<FakeDataInputs<T> >& getInputs(void) const {return m_inputs;}
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

      ENSEM::EnsemComplex getZ(const SEMBLE::SembleVector<std::complex<double> > &zz, const int state)
      {
        return zz(state);
      } 


      private:
      Fake3ptCorr(void);

      Fake3ptKey m_key;
      typename SEMBLE::PromoteEnsemVec<T>::Type m_threePtFcn;
      typename ADAT::Handle<FakeDataInputs<T> > m_inputs;
    };



}

#endif
