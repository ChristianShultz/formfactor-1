// test_make_fake_spectrum.cc -
//
// Monday, July  2 2012
//

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include "radmat/fake_data/fake_spectrum.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_matrix.h"
#include <vector>
#include <string>


namespace radmat
{
  
  tester test_make_fake_spectrum(void)
  {
    tester m_test(__func__);
    TESTER_TEST(m_test,true,"foobar");
    
    Array<double> spectrum; 
    spectrum.resize(3);
    spectrum[0] = spectrum[1] = spectrum[2] = acos(-1.)/2.;


    FakeDataIni_t fakeIni;
    fakeIni.dataProps.ncfg = 50;
    fakeIni.timeProps.tsink = 5;
    fakeIni.timeProps.tsource = 0;
    fakeIni.stateProps.sinkUpdateCovariance = true;
    fakeIni.stateProps.sinkUpdateVariance = true;
    fakeIni.stateProps.sinkVarO = 0.1;
    fakeIni.stateProps.sinkMasses = spectrum;
    
    FakeSpectrum fakeSpec(fakeIni);

    std::vector<SEMBLE::SembleVector<double> > foo = fakeSpec.generate(std::string("sink"));

    return m_test;
  }



}
