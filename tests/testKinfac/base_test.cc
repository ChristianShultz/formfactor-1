// base_test.cc -
//
// Wednesday, April 11 2012
//


#include "kinematic_factors/kinfacbase.h"
#include "itpp/itbase.h"
#include "ensem/ensem.h"

using namespace kinfac;
using namespace ENSEM;


int 
main(void)
{

  const int ncfg = 600;
  
  const double m_ = 0.21, _m = 0.21, err = 1e-2;

  kinIni ini;
  ini.massCut = 2*err + err/100.;

  KinFacABC<double> *foo = new kf0p0p<double>(ini);
  Array<int> _mom(3), mom_(3);
  _mom[0] = 1;
  _mom[1] = 0;
  _mom[2] = 0;
  mom_[0] = 0;
  mom_[1] = 1;
  mom_[2] = 1;


  EnsemReal E_, _E;
  const double e_ = sqrt(m_*m_
			 + mom_[0]*mom_[0]
			 + mom_[1]*mom_[1]
			 + mom_[2]*mom_[2]
			 );
  const double _e = sqrt(_m*_m
			 + _mom[0]*_mom[0]
			 + _mom[1]*_mom[1]
			 + _mom[2]*_mom[2]
			 );

  E_.resize(ncfg);
  _E.resize(ncfg);

  for(int cfg = 0; cfg < ncfg; ++cfg)
    {
      E_.elem(cfg) = e_ + err*(2.*itpp::randu() -1.);
      _E.elem(cfg) = _e + err*(2.*itpp::randu() -1.);
    }


  std::cout << "toDouble(mean(E_)) " << toDouble(mean(E_)) << std::endl;
  std::cout << "toDouble(mean(_E)) " << toDouble(mean(_E)) << std::endl;

  KinKey_p src(E_,mom_,0,0,"+");
  KinKey_p sink(_E,_mom,0,0,"+");

  
  KinKey key(sink,src);

  // doesn't actually do anything in this example
  registerPolarisationTensors(key);
 
  
  LinSysRetType<double>::type sys;

  sys = (*foo)(key);

  delete foo;


  std::cout << sys[0].mean() << std::endl;
  std::cout << sys[1].mean() << std::endl;
  std::cout << sys[2].mean() << std::endl;
  std::cout << sys[3].mean() << std::endl;
  
  return 0;
}
