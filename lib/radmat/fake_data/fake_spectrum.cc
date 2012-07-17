// fake_spectrum.cc -
//
// Friday, June 29 2012
//

#include "io/adat_xmlio.h"
#include "fake_spectrum.h"
#include "fake_data_ini.h"
#include "covarrying_vectors.h"
#include "semble/semble_vector.h"
#include "itpp/itbase.h"
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"
#include "io/adat_xmlio.h"
#include <vector>
#include <string>

using namespace ADATXML;
using namespace ADATIO;

typedef std::vector<SEMBLE::SembleVector<double> > mt_t;
typedef std::vector<itpp::Vec<double> > mitpp_t;

namespace 
{

  struct pack_t
  {
    int ncfg;
    bool uCov;
    bool uVar;
    double varO;
    int Lt;
    Array<double> spectrum;
  };


  pack_t pull(const radmat::FakeDataIni_t &ini, const std::string &s)
  {
    pack_t pack;
    pack.ncfg = ini.dataProps.ncfg;
    pack.Lt = abs(ini.timeProps.tsink - ini.timeProps.tsource + 1);

    if(s == std::string("sink"))
    {
      pack.uCov = ini.stateProps.mProps.sinkUpdateCovariance;
      pack.uVar = ini.stateProps.mProps.sinkUpdateVariance;
      pack.varO = ini.stateProps.mProps.sinkVarO;
      pack.spectrum = ini.stateProps.mProps.sinkMasses;

      return pack;
    }
    else if(s == std::string("source"))
    {
      pack.uCov = ini.stateProps.mProps.sourceUpdateCovariance;
      pack.uVar = ini.stateProps.mProps.sourceUpdateVariance;
      pack.varO = ini.stateProps.mProps.sourceVarO;
      pack.spectrum = ini.stateProps.mProps.sourceMasses;

      return pack;
    }
    else
      SPLASH("failed to specify source or sink");
    exit(1);
  }

}



namespace radmat
{

  std::vector<SEMBLE::SembleVector<double> > FakeSpectrum::generate(const std::string &s) const
  {
    pack_t pack = pull(m_ini,s);
    int dim = pack.spectrum.size();
    SEMBLE::SembleVector<double> zero(pack.ncfg,dim);
    zero.zeros();
    std::vector<SEMBLE::SembleVector<double> > Lambdat(pack.Lt,zero);
    itpp::Vec<double> mean(pack.Lt);
    itpp::Vec<double> var = itpp::abs(pack.varO*itpp::randn(pack.Lt));
    corMat genCor;
    itpp::Mat<double> cor = genCor.genMat<double>(pack.Lt);    
    std::vector<itpp::Vec<double> > work;

    for(int elem = 0; elem < dim; elem++)
    {
      mean = pack.spectrum[elem];

      if(pack.uCov)
        cor = genCor.genMat<double>(pack.Lt);
      if(pack.uVar)
        var = itpp::abs(pack.varO * itpp::randn(pack.Lt));

      work = genCovarryingDist(mean,var,cor,pack.ncfg);

      for(int t = 0; t < pack.Lt; t++)
        for(int bin = 0; bin < pack.ncfg; bin++)
          Lambdat[t].setElement(bin,elem,work[bin][t]);
    }

    //  std::cout << __func__ << std::endl;
    //  std::vector<SEMBLE::SembleVector<double> >::const_iterator it;
    //  for(it = Lambdat.begin(); it != Lambdat.end(); it++)
    //      std::cout << it->mean() << std::endl;

    return Lambdat;
  }




}
