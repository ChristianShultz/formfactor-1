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

      Fake3ptCorr(const rHandle<FakeDataInputs<T> > &inputs, 
          const int lorentz, 
          const int hel_source, 
          const int hel_sink, 
          const pProp_t &mom);

      Fake3ptCorr(const Fake3ptCorr<T> &o);

      Fake3ptCorr<T>& operator=(const Fake3ptCorr<T> &o);

      public:
      // peek and poke the key
      const Fake3ptKey& getKey(void) const {return m_key;}
      Fake3ptKey& getKey(void) {return m_key;}

      // peek and poke the 3pt corr
      const typename SEMBLE::PromoteEnsemVec<T>::Type& get3pt(void) const {return m_threePtFcn;}
      typename SEMBLE::PromoteEnsemVec<T>::Type& get3pt(void) {return m_threePtFcn;}

      // peek and poke the input parameters
      const rHandle<FakeDataInputs<T> >& getInputs(void) const {return m_inputs;}
       rHandle<FakeDataInputs<T> >& getInputs(void) {return m_inputs;}


      // get one of the guys in the sum 
      typename SEMBLE::PromoteEnsemVec<T>::Type get3ptComponent(const int source_idx, const int sink_idx) const; 
      int getNSource(void) const;
      int getNSink(void) const;
      int getLt(void) const;
      int getNCfg(void) const;

      private:  
      ENSEM::EnsemReal getE(const std::vector<SEMBLE::SembleVector<double> > &spec, const int target);
      ENSEM::EnsemComplex getZ(const std::vector< SEMBLE::SembleVector<std::complex<double> > > &zz, const int target);

      private:
      Fake3ptCorr(void); // hide ctor

      Fake3ptKey m_key;
      typename SEMBLE::PromoteEnsemVec<T>::Type m_threePtFcn; 
      typename std::vector<SEMBLE::SembleMatrix<T> > m_three; // vector index is time, matrix is source, sink
      rHandle<FakeDataInputs<T> > m_inputs;
    };


  // impl

  template<typename T>
    Fake3ptCorr<T>::Fake3ptCorr(const rHandle<FakeDataInputs<T> > &inputs, 
        const int lorentz, 
        const int hel_source, 
        const int hel_sink, 
        const pProp_t &mom)

    {

#pragma omp critical
      {
        m_inputs = inputs; 
      }

      double meanvar = (inputs->original->ini.stateProps.mProps.sourceVarO
          + inputs->original->ini.stateProps.mProps.sinkVarO)/2.;


      // fill in the key
      m_key.lorentz = lorentz;
      m_key.hel_source = hel_source;
      m_key.hel_sink = hel_sink;
      m_key.mom = mom;
      m_key.E_source = getE(inputs->working->specsource,inputs->working->ini.matElemProps.right_target);
      m_key.E_sink = getE(inputs->working->specsink,inputs->working->ini.matElemProps.left_target);
      m_key.Z_source = getZ(inputs->working->zsource,inputs->working->ini.matElemProps.right_target);
      m_key.Z_sink = getZ(inputs->working->zsink,inputs->working->ini.matElemProps.left_target);

      SemblePInv mom_inv =   makeMomInvariants( m_key.E_sink,
          m_key.E_source,
          mom.momSink,
          mom.momSource,
          inputs->working->mom_factor);

      // try to figure out if we are looking at a diagonal guy or not -- duct tape and dreams
      SEMBLE::SembleVector<double> p1(mom_inv.pf()), p2(mom_inv.pi());
      double m1,m2;
      m1 = sqrt(SEMBLE::toScalar(ENSEM::mean(p1(0)*p1(0) - p1(1)*p1(1) - p1(2)*p1(2) - p1(3)*p1(3))));
      m2 = sqrt(SEMBLE::toScalar(ENSEM::mean(p2(0)*p2(0) - p2(1)*p2(1) - p2(2)*p2(2) - p2(3)*p2(3))));
      if( fabs(m1 - m2)/std::max(m1,m2) <= meanvar)
        m_key.elemIDBase = inputs->working->ini.matElemProps.diag;
      else  
        m_key.elemIDBase = inputs->working->ini.matElemProps.off;


      ENSEM::EnsemReal Q2 = mom_inv.Q2(); 

      m_key.Q2 = Q2;
      m_key.t_source = inputs->working->ini.timeProps.tsource;
      m_key.t_sink = inputs->working->ini.timeProps.tsink;
      m_key.mom_factor = inputs->working->mom_factor;

      // fake up the three pt function
      m_three.clear();
      for(unsigned int t = 0; t < inputs->working->specsource.size(); t++)
        m_three.push_back( makeFakeDataPoint(inputs,mom,hel_sink,hel_source,t,lorentz) );

      m_threePtFcn.resize(m_three[0].getB());
      m_threePtFcn.resizeObs(m_three.size());

      typename SEMBLE::PromoteEnsem<T>::Type  sum; 
      sum.resize(m_three[0].getB());


      for(unsigned int t=0; t < m_three.size(); ++t)
      {
        sum = SEMBLE::toScalar(T(0)); 
        for(int row = 0; row < m_three[t].getN(); ++row)
          for(int col = 0; col < m_three[t].getM(); ++col)
            sum += m_three[t].getEnsemElement(row,col);

        ENSEM::pokeObs(m_threePtFcn,sum,t);
      }


    }

  template<typename T>
    Fake3ptCorr<T>::Fake3ptCorr(const Fake3ptCorr<T> &o)
    :m_key(o.m_key) , m_threePtFcn(o.m_threePtFcn)  , m_three(o.m_three) 
    {

#pragma omp critical
      {
        m_inputs = o.m_inputs; 
      }
    }


  template<typename T>
    Fake3ptCorr<T>& Fake3ptCorr<T>::operator=(const Fake3ptCorr<T> &o)
    {
      if(this != &o)
      {
        m_key = o.m_key;
        m_threePtFcn = o.m_threePtFcn;
        m_three = o.m_three;
#pragma omp critical
        {
          m_inputs = o.m_inputs;
        }
      }

      return *this;
    }


  template<typename T>
    int Fake3ptCorr<T>::getNSource(void) const
    {
      POW2_ASSERT(!!!m_three.empty()); 
      return m_three[0].getN();
    }

  template<typename T>
    int Fake3ptCorr<T>::getNSink(void) const 
    {
      POW2_ASSERT(!!!m_three.empty()); 
      return m_three[0].getM();
    }

  template<typename T>
    int  Fake3ptCorr<T>::getLt(void) const 
    {
      POW2_ASSERT(!!!m_three.empty()); 
      return m_three.size();
    }

  template<typename T>
    int  Fake3ptCorr<T>::getNCfg(void) const 
    {
      POW2_ASSERT(!!!m_three.empty()); 
      return m_three[0].getB();
    }

  // get one of the guys in the sum 
  template<typename T> 
    typename SEMBLE::PromoteEnsemVec<T>::Type Fake3ptCorr<T>::get3ptComponent(const int source_idx, const int sink_idx) const
    {
      typename SEMBLE::PromoteEnsemVec<T>::Type ret; 
      ret = m_threePtFcn;
      ret = SEMBLE::toScalar(T(0)); 

      for(unsigned int t = 0; t < m_three.size(); ++t)
        ENSEM::pokeObs(ret,m_three[t].getEnsemElement(source_idx,sink_idx),t);

      return ret; 
    } 



  // private
  template<typename T>
    ENSEM::EnsemReal Fake3ptCorr<T>::getE(const std::vector<SEMBLE::SembleVector<double> > &spec, const int target)
    {
      ENSEM::EnsemReal ret; 
      ret.resize(spec[0].getB());
      ret = SEMBLE::toScalar(double(0));
      std::vector<SEMBLE::SembleVector<double> >::const_iterator it;
      for(it = spec.begin(); it != spec.end(); it++)
        ret += (*it)(target);

      ret /= SEMBLE::toScalar(double(spec.size()));

      return ret;      
    }


  template<typename T>
    ENSEM::EnsemComplex Fake3ptCorr<T>::getZ(const std::vector< SEMBLE::SembleVector<std::complex<double> > > &zz, const int target)
    {
      ENSEM::EnsemComplex ret; 
      ret.resize(zz[0].getB());
      ret = SEMBLE::toScalar(std::complex<double>(0,0));
      std::vector<SEMBLE::SembleVector<std::complex<double> > >::const_iterator it;
      for(it = zz.begin(); it != zz.end(); ++it)
        ret += (*it)(target);

      ret /= SEMBLE::toScalar(double(zz.size()));

      return ret; 
    } 



}

#endif
