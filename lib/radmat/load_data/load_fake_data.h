#ifndef LOAD_FAKE_DATA_H_H_GUARD
#define LOAD_FAKE_DATA_H_H_GUARD

#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/fake_data/gen_fake_dataset.h"
#include "three_point.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include <vector>

namespace radmat
{

  template<typename T>
    struct LoadFake3pt
    { 

      typedef typename GenFakeDataSet<T>::ffFunction ffFunction;

      LoadFake3pt(const FakeDataIni_t &ini)
        : m_ini(ini) 
      {  }

      std::vector<ThreePointCorrelator<T> > genData(void)
      {
        typename ADAT::Handle<std::vector< typename GenFakeDataSet<T>::Corr> > fake_data;
        GenFakeDataSet<T> gen_data(m_ini);
        gen_data.generate();
        fake_data = gen_data.get();

        m_ff_inputs = gen_data.get_FF_inputs();

        typename std::vector< typename GenFakeDataSet<T>::Corr >::const_iterator it;
        typename std::vector<ThreePointCorrelator<T> > ret;

        for(it = fake_data->begin(); it != fake_data->end(); it++)
          ret.push_back(makeThreePointFromFake(*it));

        return ret;
      }


      typename std::vector<ffFunction> get_FF_inputs(void) const
      {
        return m_ff_inputs;
      }


      private:
      LoadFake3pt(void);

      FakeDataIni_t m_ini;
      typename std::vector<ffFunction> m_ff_inputs;

    };


} // radmat namespace


#endif
