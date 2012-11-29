#ifndef FAKE_3PT_FUNCTION_AUX_H_H_GUARD
#define FAKE_3PT_FUNCTION_AUX_H_H_GUARD

#include "ensem/ensem.h"
#include "semble/semble_semble.h"
#include "radmat/utils/aux.h"
#include "radmat/utils/harmonic.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include "radmat/ff/formfactor_factory.h"
#include "radmat/ff/ff_gen_llsq_row.h"
#include "adat/handle.h"
#include "fake_data_ini.h"
#include "fake_spectrum.h"
#include "fake_overlaps.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace radmat
{


  // FWD
  ///////

  struct FakeMatrixElement;
  struct Fake3ptFactor;
  struct Fake3ptSumElement;
  template<typename T>
    struct FakeDataInputs_p;
  template<typename T>
    struct FakeDataInputs; 


  // make a copy of the originals and apply dispersion etc to them
  template<typename T>
    typename ADAT::Handle< FakeDataInputs<T> >
    copyFakeInput(const typename ADAT::Handle<FakeDataInputs_p<T> > &orig);

  // get the first round
  template<typename T>
    typename ADAT::Handle<FakeDataInputs_p<T> > generateOriginalInputs(const FakeDataIni_t &ini);

  // apply Z suppression.. basically try to make something that looks like optimized operators were used
  // its a hack and i don't know at the time of typing this if its even worth the effort
  template<typename T>
    void applyZSuppression(typename ADAT::Handle<FakeDataInputs_p<T> > &input);

  // assume all states are stable against decay and use a dispersion relation to 
  // determine the energies at non-zero momentum 
  // aaEE = aamm + 1/(xi*xi) *(2pi/L_s)^2 * p*p   || a is the temporal spacing
  template<typename T>
    void applyDispersion(typename ADAT::Handle<FakeDataInputs_p<T> > &input, const pProp_t &mom);

  template<typename T>
    SEMBLE::SembleMatrix<T>  
    makeFakeDataPoint(const ADAT::Handle<FakeDataInputs<T> > &inputs,
        const pProp_t &mom,
        const int hel_sink,
        const int hel_source,
        const int t_ins,
        const int lorentz);


  // SOME FUNCTOR CLASSES TO MAKE LIFE EASIER
  /////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////




  // #define DEBUG_FAKE_MAT_ELEM_CJS

  // NB: moved matrix extra time dependence into propagation factor
  struct FakeMatrixElement
  {
    typedef bind1st_2ParFunction_cc<double,int,double,&One> ffFunction;

    ENSEM::EnsemComplex operator()(const std::string &elemID, 
        const int lorentz,
        const SemblePInv &mom,
        const std::vector<ffFunction> &ffgen,
        const double Q2)
    {
      ENSEM::EnsemComplex matelem;
      ffKinematicFactors_t<std::complex<double> >
        K(FormFactorDecompositionFactoryEnv::callFactory(elemID));

      SEMBLE::SembleVector<std::complex<double> > K_k = (K.genFactors(mom)).getRow(lorentz);
      matelem.resize(K_k.getB());
      matelem = SEMBLE::toScalar(std::complex<double>(0.,0.));

      POW2_ASSERT(K_k.getN() == (int)ffgen.size());

      for(unsigned int k = 0; k < ffgen.size(); k++)
        matelem += K_k(k) * SEMBLE::toScalar(ffgen[k](Q2));

#ifdef DEBUG_FAKE_MAT_ELEM_CJS
      // -- DEBUG
      if(SEMBLE::toScalar(ENSEM::mean(matelem) ) != 0.)
      {
        std::cout << __func__ << std::endl;
        std::cout << "matelem " << SEMBLE::toScalar(ENSEM::mean(matelem)) << std::endl;
        std::cout << "exp " << SEMBLE::toScalar(ENSEM::mean(exp)) << std::endl;
        std::cout << "mom.first(0) " << SEMBLE::toScalar(ENSEM::mean(mom.first(0))) << std::endl;
        std::cout << "mom.second(0) " << SEMBLE::toScalar(ENSEM::mean(mom.second(0))) << std::endl;
        std::cout << "K_k(0) " << SEMBLE::toScalar(ENSEM::mean(K_k(0))) << std::endl;
        std::cout << "ffgen[0](Q2) " << ffgen[0](Q2) << std::endl;
        std::cout << "Q2 " << Q2 << std::endl;
        std::cout << std::endl;
      }
#endif

      return matelem ;
    }
  };

#undef DEBUG_FAKE_MAT_ELEM_CJS



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<typename T> 
    struct ThreePtPropagationFactor
    {
      typename SEMBLE::PromoteEnsem<T>::Type operator()(const ENSEM::EnsemReal &E_sink,
          const typename SEMBLE::PromoteEnsem<T>::Type &Z_sink,
          const int t_sink,
          const int t_ins,
          const ENSEM::EnsemReal &E_source, 
          const typename SEMBLE::PromoteEnsem<T>::Type &Z_source, 
          const int t_source)
      {
        return ((ENSEM::conj(Z_source)*Z_sink/ (E_source * E_sink * SEMBLE::toScalar(4.)))
            * ENSEM::exp(-E_sink*(SEMBLE::toScalar(double(t_sink - t_ins))))
            * ENSEM::exp(-E_source*(SEMBLE::toScalar(double(t_ins - t_source))))
            );

      }

    };



#if 0
  // returns 1/4EnEm exp(-E_n t_n) exp(E_m t_m) Z_n^* Z_m ~ all the factors not 
  // associated with the matrix elem
  // NB : we do a complex conjugation in here !!!

  // DONT USE ME!!!
  struct Fake3ptFactor
  {

    ENSEM::EnsemComplex operator()(const ENSEM::EnsemReal &E_n,
        const ENSEM::EnsemComplex &Z_n,
        const int t_n,
        const ENSEM::EnsemReal &E_m,
        const ENSEM::EnsemComplex &Z_m,
        const int t_m)
    {
      __builtin_trap();
      return  ((ENSEM::exp(-E_n * SEMBLE::toScalar(double(t_n))) 
            * ENSEM::exp(E_m * SEMBLE::toScalar(double(t_m)))
            * ENSEM::conj(Z_n)
            * Z_m) 
          / ( E_n * E_m * SEMBLE::toScalar(double(4))));

    }
  };



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  // NB: dont use me!!!!!
  struct Fake3ptSumElement
  {

    ENSEM::EnsemComplex operator()(const std::string &elemID,
        const int t_sink, 
        const int t_ins,
        const int t_source,
        const int lorentz, 
        const std::vector<FakeMatrixElement::ffFunction> &ffgen,
        const SemblePInv &mom,
        const double Q2, 
        const ENSEM::EnsemComplex &Z_sink,
        const ENSEM::EnsemComplex &Z_source)
    {
      __builtin_trap();
      FakeMatrixElement matelem;
      Fake3ptFactor factor;

      return (matelem(elemID,t_ins,lorentz,mom,ffgen,Q2) 
          * factor(mom.first(0), Z_sink,t_sink, mom.second(0),Z_source,t_source));
    }

  };

#endif 

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<typename T>
    struct FakeDataInputs_p
    {
      typename std::vector<SEMBLE::SembleVector<T> > zsource; // vector index is time, semble index is state
      typename std::vector<SEMBLE::SembleVector<T> > zsink;   // same as above
      std::vector<SEMBLE::SembleVector<double> > specsource;  // ''
      std::vector<SEMBLE::SembleVector<double> > specsink;    // ''
      itpp::Mat<std::vector<FakeMatrixElement::ffFunction> > ffgenerator;
      FakeDataIni_t ini;
      double mom_factor;
    };

  template<typename T>
    struct FakeDataInputs
    {
      typename ADAT::Handle<FakeDataInputs_p<T> > working;   // after dispersion and suppression
      typename ADAT::Handle<FakeDataInputs_p<T> > original;  // the original guy
    };



  // FUNCTION IMPLEMENTATION 
  ////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  template<typename T>
    ADAT::Handle<FakeDataInputs_p<T> > generateOriginalInputs(const FakeDataIni_t &ini)
    {
      ADAT::Handle<FakeDataInputs_p<T> > handle(new FakeDataInputs_p<T> );
      POW2_ASSERT(&*handle);

      // Z
      FakeOverlaps FakeLaps(ini);
      int szSource = ini.stateProps.mProps.sourceMasses.size();
      int szSink = ini.stateProps.mProps.sinkMasses.size();

      handle->zsource = pullRow(FakeLaps.generate<T>(szSource,std::string("source")),0);
      if(!!!ini.stateProps.sameOp)
        handle->zsink = pullRow(FakeLaps.generate<T>(szSink,std::string("sink")),0);
      else
        if( (szSource == szSink) &&
            (ini.stateProps.mProps.sourceMasses == ini.stateProps.mProps.sinkMasses))
          handle->zsink = handle->zsource;
        else
        {
          SPLASH("incorrect specification of source and sink masses");
          SPLASH("stateProps.source/sinkMasses conflicts with matElemProps.sameOp, exiting");
          exit(1);
        }

      // spectrum
      FakeSpectrum FakeSpec(ini);

      handle->specsource = FakeSpec.generate(std::string("source"));

      /*   // DEBUG
           std::cout << __func__ << std::endl;
           std::vector<SEMBLE::SembleVector<double> >::const_iterator it;

           std::cout << "size: " << handle->specsource.size() << std::endl;
           for(it = handle->specsource.begin(); it != handle->specsource.end(); it++)
           std::cout << "B: " << it->getB() << "\n N: " << it->getN() << "\n mean: " << it->mean() << std::endl;

           __builtin_trap();

       */

      if(!!!ini.stateProps.sameOp)
        handle->specsink = FakeSpec.generate(std::string("sink"));
      else
        handle->specsink = handle->specsource;

      // ini
      handle->ini = ini;
      handle->mom_factor = 1./ini.dispersionProps.xi * 2. * acos(-1.)/ini.dispersionProps.L_s;

      // ffgenerator
      itpp::Mat<std::vector<FakeMatrixElement::ffFunction> > ffmv(szSource,szSink);

      typename ADAT::Handle<ffBase_t<T> > foobar =
        FormFactorDecompositionFactoryEnv::callFactory(ini.matElemProps.diag);
      int nff = foobar->nFacs();

      for(int n = 0; n < szSource; n++)
        for(int m = 0; m < szSink; m++)
        {
          std::vector<FakeMatrixElement::ffFunction> ffv;
          for(int k =0; k < nff; k++)
          {
            FakeMatrixElement::ffFunction bar;
            bar.bind1st(itpp::randi(0,12));
            ffv.push_back(bar);
          }
          ffmv(n,m) = ffv;
        }

      // non-square symmetric..
      for(int n = 0 ; n < szSource; ++n)
        if(n < szSink)
          for(int m = n ; m < szSink; ++m)
            if(m < szSource)
              ffmv(m,n) = ffmv(n,m);


      handle->ffgenerator = ffmv;
      return handle;
    }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<typename T>
    typename ADAT::Handle<FakeDataInputs<T> >
    copyFakeInput(const typename ADAT::Handle<FakeDataInputs_p <T> > &orig)
    {

      FakeDataInputs<T> *ret = new FakeDataInputs<T>();
      ret->original = orig;
      ret->working = ADAT::Handle<FakeDataInputs_p<T> >(new FakeDataInputs_p<T>);

      // need to copy 'by hand' here since we want a "new" set of inputs
      ret->working->zsource = orig->zsource;
      ret->working->zsink = orig->zsink;
      ret->working->specsource = orig->specsource;
      ret->working->specsink = orig->specsink;
      ret->working->ffgenerator = orig->ffgenerator;
      ret->working->ini = orig->ini;
      ret->working->mom_factor = orig->mom_factor;

      return ADAT::Handle<FakeDataInputs<T> >(ret);
    }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<typename T>
    void applyZSuppression(ADAT::Handle<FakeDataInputs_p<T> > &input)
    {
      // NB don't do anything if we're reading them from an xml file
      if(!!!input->ini.stateProps.readZ)
      {
        int left = input->ini.matElemProps.left_target;
        int right = input->ini.matElemProps.right_target;
        SEMBLE::PromoteScalar<double>::Type suppression;
        suppression = SEMBLE::toScalar(double(input->ini.stateProps.zProps.suppressionOrder));
        typename std::vector<SEMBLE::SembleVector<T> >::iterator it;


        if(input->ini.stateProps.zProps.suppress)
        {
          for(it = input->zsource.begin(); it != input->zsource.end(); it++)
            for(int elem = 0; elem < it->getN(); elem++)
              if(elem == right)
              {
                if(input->ini.stateProps.zProps.targetZ_at_order1)
                  it->loadEnsemElement(elem,(it->getEnsemElement(elem))/(it->getEnsemElement(elem)));
              }
              else
                it->loadEnsemElement(elem,suppression*it->getEnsemElement(elem));

          for(it = input->zsink.begin(); it != input->zsink.end(); it++)
            for(int elem = 0; elem < it->getN(); elem++)
              if(elem == left)
              {
                if(input->ini.stateProps.zProps.targetZ_at_order1)
                  it->loadEnsemElement(elem,(it->getEnsemElement(elem))/(it->getEnsemElement(elem)));
              }
              else
                it->loadEnsemElement(elem,suppression*it->getEnsemElement(elem));
        }
        else if(input->ini.stateProps.zProps.targetZ_at_order1)
        {	

          for(it = input->zsource.begin(); it != input->zsource.end(); it++)
            it->loadEnsemElement(right,(it->getEnsemElement(right))/(it->getEnsemElement(right)));

          for(it = input->zsink.begin(); it != input->zsink.end(); it++)
            it->loadEnsemElement(left,(it->getEnsemElement(left))/(it->getEnsemElement(left)));
        }
      }
    }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  // assume all states are stable against decay and use a dispersion relation to 
  // determine the energies at non-zero momentum 
  // aaEE = aamm + 1/(xi*xi) *(2pi/L_s)^2 * p*p   || a is the temporal spacing
  template<typename T>
    void applyDispersion(typename ADAT::Handle<FakeDataInputs_p<T> > &input, const pProp_t &mom)
    {
      SEMBLE::PromoteScalar<double>::Type factor,sink,source;
      factor = SEMBLE::toScalar(( 2. * acos(-1.))
          /double(input->ini.dispersionProps.L_s)
          /(input->ini.dispersionProps.xi));         // 1/xi * 2pi/L_s

      double mom2source = ( mom.momSource[0] * mom.momSource[0]
          +mom.momSource[1] * mom.momSource[1]
          +mom.momSource[2] * mom.momSource[2]);

      double mom2sink = ( mom.momSink[0] * mom.momSink[0]
          + mom.momSink[1] * mom.momSink[1]
          + mom.momSink[2] * mom.momSink[2]);


      sink = factor*factor*SEMBLE::toScalar(mom2sink);
      source = factor*factor*SEMBLE::toScalar(mom2source);

      std::vector<SEMBLE::SembleVector<double> >::iterator it;
      for(it = input->specsource.begin();  it != input->specsource.end(); it++)
        for(int elem = 0; elem < it->getN(); elem++)
          it->loadEnsemElement(elem,sqrt(source + (*it)(elem) * (*it)(elem)));

      for(it = input->specsink.begin(); it != input->specsink.end(); it++)
        for(int elem = 0; elem < it->getN(); elem++)
          it->loadEnsemElement(elem,sqrt(sink + (*it)(elem)*(*it)(elem)));
    }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<typename T>
    SEMBLE::SembleMatrix<T> 
    makeFakeDataPoint(const typename ADAT::Handle<FakeDataInputs<T> > &inputs,
        const pProp_t &mom,
        const int hel_sink,
        const int hel_source,
        const int t_ins,
        const int lorentz)
    {
      const SEMBLE::SembleVector<T> *zsource, *zsink;
      const SEMBLE::SembleVector<double> *specsink, *specsource, *specsink_ins, *specsource_ins;
      typename ADAT::Handle<FakeDataInputs_p<T> > handle(inputs->working);

      zsource = &(handle->zsource[handle->ini.timeProps.tsource]);
      zsink = &(handle->zsink[handle->ini.timeProps.tsink]);
      specsource_ins = &(handle->specsource[t_ins]);
      specsink_ins = &(handle->specsink[t_ins]);
      specsource = &(handle->specsource[handle->ini.timeProps.tsource]);
      specsink = &(handle->specsink[handle->ini.timeProps.tsink]);


      const int nsource = specsource->getN();
      const int nsink = specsink->getN();

      // Fake3ptFactor facGen;
      ThreePtPropagationFactor<T> facGen;
      FakeMatrixElement matGen;
     SEMBLE::SembleMatrix<T> ret(specsink->getB(),nsource,nsink);
      ret.zeros();


      for(int source = 0; source < nsource; source++)
        for(int sink = 0; sink < nsink; sink++)
        {
          std::stringstream ss;

          double meanvar;
          meanvar = (handle->ini.stateProps.mProps.sourceVarO + handle->ini.stateProps.mProps.sinkVarO)/2.;


          // there is some implicit covariance structure built in here by using covarrying energy draws
          SemblePInv mom_inv =   makeMomInvariants( (*specsink_ins)(sink),
              (*specsource_ins)(source),
              mom.momSink,
              mom.momSource,
              handle->mom_factor);

          ENSEM::EnsemReal Q2 = mom_inv.Q2();



          // try to figure out if we are looking at a diagonal guy or not -- duct tape and dreams
          SEMBLE::SembleVector<double> p1(mom_inv.pf()), p2(mom_inv.pi());
          double m1,m2;
          m1 = sqrt(SEMBLE::toScalar(ENSEM::mean(p1(0)*p1(0) - p1(1)*p1(1) - p1(2)*p1(2) - p1(3)*p1(3))));
          m2 = sqrt(SEMBLE::toScalar(ENSEM::mean(p2(0)*p2(0) - p2(1)*p2(1) - p2(2)*p2(2) - p2(3)*p2(3))));

          /* -- DEBUG          
             std::cout << __func__ << std::endl;
             std::cout << m1 - m2 << std::endl;
             std::cout << fabs(m1 - m2)/std::max(m1,m2) << " " << meanvar << std::endl;
           */ 

          if( fabs(m1 - m2)/std::max(m1,m2) <= meanvar)
            ss << handle->ini.matElemProps.diag;
          else 
            ss << handle->ini.matElemProps.off;

          ss << "_" << hel_source << "_" << hel_sink;


          ENSEM::EnsemComplex mat;
          mat = matGen(ss.str(),
              lorentz,
              mom_inv,
              handle->ffgenerator(source,sink),
              SEMBLE::toScalar(ENSEM::mean(Q2)));

          ENSEM::EnsemComplex factor;
          factor = facGen( (*specsink)(sink), (*zsink)(sink), handle->ini.timeProps.tsink,
              t_ins,(*specsource)(source), (*zsource)(source), handle->ini.timeProps.tsource);

          ret.loadEnsemElement(source,sink, (factor*mat));
        }  

      return ret;

    }

} // namespace

#endif

