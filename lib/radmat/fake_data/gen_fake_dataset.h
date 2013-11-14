#ifndef GEN_FAKE_DATASET_H_H_GUARD
#define GEN_FAKE_DATASET_H_H_GUARD


#include "fake_3pt_function.h"
#include "fake_3pt_function_aux.h"
#include "ensem/ensem.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include "radmat/utils/handle.h"
#include "semble/semble_file_management.h"
#include <string>
#include <sstream>
#include <vector>



namespace radmat
{

  template<typename T> 
    struct GenFakeDataSet
    {
      typedef Fake3ptCorr<T> Corr;
      typedef  rHandle<FakeDataInputs_p<T> > PrimInput_h;
      typedef rHandle<FakeDataInputs<T> > Input_h;
      typedef FakeMatrixElement::ffFunction ffFunction;


      GenFakeDataSet(void);
      GenFakeDataSet(const FakeDataIni_t &);


      void generate(void);
      void generate(const FakeDataIni_t &);
      void write_spectrum_log(const Input_h &iptm, const pProp_t &mom);  
      std::vector<ffFunction> get_FF_inputs(void) const {return m_ff_input_functions;}

     rHandle<std::vector<Corr> > get(void);

      private:
      bool m_haveIni;
      bool m_have3Pt;
      FakeDataIni_t m_ini;
      rHandle<std::vector<Corr> > m_3Pt;
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
     rHandle<std::vector<Corr> > foo(new std::vector<Corr>() );
      m_3Pt = foo;
      POW2_ASSERT(&*m_3Pt);

      // cook up some three point functions -- since handle isn't thread safe we can't
      // use omp here which sucks a lot since it should be possible
      PrimInput_h orig = generateOriginalInputs<T>(m_ini);
      Array<int> h_f = m_ini.dataProps.hel_sink;
      Array<int> h_i = m_ini.dataProps.hel_source;
      Array<pProp_t> moms = m_ini.dataProps.momenta;

      int p, sz;
      sz = moms.size(); 
#if 0
#pragma omp parallel for shared(p,sz) 
#endif 

      for( p = 0; p < sz; p++)
      {
        Input_h work; 
#pragma omp critical
        {
          work = copyFakeInput(orig);
        }
        applyZSuppression(work->working);
        applyDispersion(work->working,moms[p]);
        write_spectrum_log(work,moms[p]); 


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
    void GenFakeDataSet<T>::write_spectrum_log(const Input_h &ipt, const pProp_t &mom)
    {
      std::string pth = SEMBLE::SEMBLEIO::getPath();
      pth += std::string("fake_data_logs");
      SEMBLE::SEMBLEIO::makeDirectoryPath(pth);
      std::stringstream ss; 
      ss << "/pi_" << mom.momSource[0] << "_" << mom.momSource[1] << "_" << mom.momSource[2] 
        << "__pf_" << mom.momSink[0] << "_" << mom.momSink[1] << "_" << mom.momSink[2];
      pth += ss.str();

      std::ofstream out(pth.c_str());

      SEMBLE::SembleVector<double> source_rest(ipt->working->specsink[0]), source_flight, sink_rest(ipt->working->specsource[0]), sink_flight; 
      source_rest = 0;
      source_flight = source_rest;
      sink_rest = 0;
      sink_flight = sink_rest;

      std::vector<SEMBLE::SembleVector<double> >::const_iterator it; 

      for(it = ipt->original->specsource.begin(); it != ipt->original->specsource.end(); ++it)
        source_rest += *it; 
      for(it = ipt->original->specsink.begin(); it != ipt->original->specsink.end(); ++it)
        sink_rest += *it; 
      for(it = ipt->working->specsource.begin(); it != ipt->working->specsource.end(); ++it)
        source_flight += *it;
      for(it = ipt->working->specsink.begin(); it != ipt->working->specsink.end(); ++it)
        sink_flight += *it; 

      source_rest /= double(ipt->original->specsource.size());
      source_flight /= double(ipt->working->specsource.size());
      sink_rest /= double(ipt->original->specsink.size());
      sink_flight /= double(ipt->working->specsink.size());

      out << "mean source rest \n" << source_rest.mean() << "\n";
      out << "mean source flight \n" << source_flight.mean() << "\n";
      out << "mean sink rest \n" << sink_rest.mean() << "\n";
      out << "mean sink flight \n" << sink_flight.mean() << "\n";


      out.close();
    }  



  template<typename T>
    void GenFakeDataSet<T>::generate(const FakeDataIni_t &ini)
    {
      m_ini = ini;
      m_haveIni = true;
      generate();
    }


  template<typename T>
    rHandle<std::vector< typename GenFakeDataSet<T>::Corr> > GenFakeDataSet<T>::get(void)
    {
      if(m_have3Pt)
        return m_3Pt;

      POW2_ASSERT(m_haveIni);
      generate();  
      return m_3Pt;    
    }

} // close namespace radmat





#endif
