#ifndef GEN_FAKE_DATASET_H_H_GUARD
#define GEN_FAKE_DATASET_H_H_GUARD


#include "fake_3pt_function.h"
#include "fake_3pt_function_aux.h"
#include "ensem/ensem.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include "adat/handle.h"
#include <vector>



namespace radmat
{

  template<typename T> 
    struct GenFakeDataSet
    {
      typedef Fake3ptCorr<T> Corr;
      typedef typename ADAT::Handle<FakeDataInputs_p<T> > PrimInput_h;
      typedef typename ADAT::Handle<FakeDataInputs<T> > Input_h;
      typedef FakeMatrixElement::ffFunction ffFunction;


      GenFakeDataSet(void);
      GenFakeDataSet(const FakeDataIni_t &);


      void generate(void);
      void generate(const FakeDataIni_t &);
      std::vector<ffFunction> get_FF_inputs(void) const {return m_ff_input_functions;}

      typename ADAT::Handle<std::vector<Corr> > get(void);

      private:
      bool m_haveIni;
      bool m_have3Pt;
      FakeDataIni_t m_ini;
      typename ADAT::Handle<std::vector<Corr> > m_3Pt;
      std::vector<ffFunction> m_ff_input_functions;

    };



  //
  // Impl
  /////////////

  template<typename T>
    GenFakeDataSet<T>::GenFakeDataSet(void)
    : m_haveIni(false) , m_have3Pt(false)
    { }


  template<typename T>
    GenFakeDataSet<T>::GenFakeDataSet(const FakeDataIni_t &ini)
    : m_haveIni(true) , m_have3Pt(false) , m_ini(ini)
    { }

  template<typename T>
    void GenFakeDataSet<T>::generate(void)
    {
      POW2_ASSERT(m_haveIni);

      // basically safely dispose of the old guy/ allocate a new guy
      typename ADAT::Handle<std::vector<Corr> > foo(new std::vector<Corr>() );
      m_3Pt = foo;
      POW2_ASSERT(&*m_3Pt);

      // cook up some three point functions -- since handle isn't thread safe we can't
      // use omp here which sucks a lot since it should be possible
      PrimInput_h orig = generateOriginalInputs<T>(m_ini);
      Array<int> h_f = m_ini.dataProps.hel_sink;
      Array<int> h_i = m_ini.dataProps.hel_source;
      Array<pProp_t> moms = m_ini.dataProps.momenta;

      for(int p = 0; p < moms.size(); p++)
      {
        Input_h work = copyFakeInput(orig);
        applyZSuppression(work->working);
        applyDispersion(work->working,moms[p]);

        for(int hf = 0; hf < h_f.size(); hf++)
          for(int hi = 0; hi < h_i.size(); hi++)
            for(int lorentz = 0; lorentz < 4; lorentz++)
              m_3Pt->push_back(  Corr(work,lorentz,h_i[hi],h_f[hf],moms[p])  );
      }

      const int left_target = orig->ini.matElemProps.left_target;
      const int right_target = orig->ini.matElemProps.right_target;
  

      m_ff_input_functions = orig->ffgenerator(left_target,right_target);

      m_have3Pt = true;
    }


  template<typename T>
    void GenFakeDataSet<T>::generate(const FakeDataIni_t &ini)
    {
      m_ini = ini;
      m_haveIni = true;
      generate();
    }


  template<typename T>
    ADAT::Handle<std::vector< typename GenFakeDataSet<T>::Corr> > GenFakeDataSet<T>::get(void)
    {
      if(m_have3Pt)
        return m_3Pt;

      POW2_ASSERT(m_haveIni);
      generate();  
      return m_3Pt;    
    }

} // close namespace radmat





#endif
