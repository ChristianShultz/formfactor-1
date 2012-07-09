// test_make_fake_overlaps.cc -
//
// Monday, July  2 2012
//

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include "radmat/fake_data/fake_overlaps.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_matrix.h"
#include <complex>
#include <vector>
#include <string>

namespace radmat
{

  tester test_make_fake_overlaps(void) 
  {
    tester m_test(__func__);
    TESTER_TEST(m_test,true,"foobar"); // just need it to go into this function so a fake test is present
    
    FakeDataIni_t fakeIni;
    fakeIni.stateProps.overlapGenerator = std::string("unitary");
    fakeIni.dataProps.ncfg = 50;
    fakeIni.matElemProps.updateCovariance = true;
    fakeIni.matElemProps.updateVariance = true;
    fakeIni.matElemProps.varianceOrder = 0.1;
    fakeIni.timeProps.tsink = 5;
    fakeIni.timeProps.tsource = 0;

    FakeOverlaps fakeLaps(fakeIni);

    std::vector<SEMBLE::SembleMatrix<double> > doublefoo = fakeLaps.generate<double>(5);
    std::vector<SEMBLE::SembleMatrix<std::complex<double> > > complexfoo = fakeLaps.generate<std::complex<double> >(5);

    return m_test;
  }





}
