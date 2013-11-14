#ifndef LOAD_FAKE_DATA_H_H_GUARD
#define LOAD_FAKE_DATA_H_H_GUARD

#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/fake_data/gen_fake_dataset.h"
#include "three_point.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/fitting/axis_plotter.h"
#include <vector>
#include <omp.h>
namespace radmat
{

  template<typename T>
    struct LoadFake3pt
    { 

      typedef typename GenFakeDataSet<T>::ffFunction ffFunction;

      LoadFake3pt(const FakeDataIni_t &ini)
        : m_ini(ini) 
      {  }

      // parallelize me please ~ 95% of the work oes throuugh me
      std::vector<ThreePointCorrelator<T> > genData(void)
      {
        rHandle<std::vector< typename GenFakeDataSet<T>::Corr> > fake_data;
        GenFakeDataSet<T> gen_data(m_ini);
        gen_data.generate();
        fake_data = gen_data.get();

        m_ff_inputs = gen_data.get_FF_inputs();

        typename std::vector< typename GenFakeDataSet<T>::Corr >::const_iterator it;
        typename std::vector<ThreePointCorrelator<T> > ret;
        ThreePointCorrelator<T> fill = makeThreePointFromFake((*fake_data)[0]); 

        // need to provide a fill else checkResize on ensem will cause failure .. why is it built like this?
        ret.resize(fake_data->size(), fill); 
        unsigned int index, stop; 
        stop = ret.size(); 

#if 0
        omp_set_num_threads(omp_get_max_threads()); 

#pragma omp parallel for shared(index,stop)

#endif 
        for(index = 0; index < stop; ++index)
        {
          MakeAxisPlots plotter; 
          plotter.plot( (*fake_data)[index] );
          ret[index] = makeThreePointFromFake( (*fake_data)[index] ); 
        }
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
