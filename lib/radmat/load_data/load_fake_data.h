#ifndef LOAD_FAKE_DATA_H_H_GUARD
#define LOAD_FAKE_DATA_H_H_GUARD

#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/fake_data/gen_fake_dataset.h"
#include "three_point.h"


namespace radmat
{

  template<typename T>
    struct LoadFake3pt
    { 
      LoadFake3pt(const FakeDataIni_t &ini)
        : m_ini(ini) 
      {  }

      std::vector<ThreePointCorrelator<T> > genData(void) const
      {
        typename ADAT::Handle<std::vector< typename GenFakeDataSet<T>::Corr> > fake_data;
        GenFakeDataSet<T> gen_data(m_ini);
        gen_data.generate();
        fake_data = gen_data.get();

        typename std::vector< typename GenFakeDataSet<T>::Corr >::const_iterator it;
        typename std::vector<ThreePointCorrelator<T> > ret;

        for(it = fake_data->begin(); it != fake_data->end(); it++)
          ret.push_back(makeThreePointFromFake(*it));

        return ret;
      }


      private:
      LoadFake3pt(void);

      FakeDataIni_t m_ini;

    };


} // radmat namespace


#endif
