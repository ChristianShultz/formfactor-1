#ifndef FAKE_DATA_H_H_GUARD
#define FAKE_DATA_H_H_GUARD

#include "io/adat_xmlio.h"
#include "fake_data_ini.h"
#include "fake_overlaps.h"
#include "fake_spectrum.h"
#include "semble/semble_meta.h"
#include "semble/semble_matrix.h"
#include "semble/semble_vector.h"
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"
#include <vector>
#include <list>

using namespace ADATXML;
using namespace ADATIO;

namespace radmat
{

  template<typename T>
  struct FakeDataPoint
  {
    typedef typename SEMBLE::PromoteEnsemVec<T> EnsemVecType;

    // vector index is time
    EnsemVecType measurement;  
    EnsemVectorReal initial_mass;
    EnsemVectorReal final_mass;
    Array<int> momSource;
    Array<int> momSink;
    std::vector<double> FF_input;
  };


  template<typename T>
  struct FakeInputs
  {
    typedef typename SEMBLE::PromoteEnsemVec<T> EnsemVecType;

    std::vector<SEMBLE::SembleMatrix<T> > zsource;
    std::vector<SEMBLE::SembleMatrix<T> > zsink;
    std::Vector<SEMBLE::SembleVector<double> > specsource;
    std::Vector<SEMBLE::SembleVector<double> > specsink;
    int Lt;
    int ncfg;
    
    void fill(void)
    {
      POW2_ASSERT( (zsource.size() == zsink.size()) &&
		   (zsink.size() == specsource.size()) && 
		   (specsource.size() == specsink.size()) );
      Lt = specsink.size();
      ncfg = ini.dataProps.ncfg;
    }

    EnsemVecType<double> pullSpec(const std::string &source_or_sink, int target)
    {
      EnsemVectorReal spec
    }
  };

  
  template<typename T>
  FakeInputs<T> getInputs(const FakeDataIni_t &ini)
  {
    FakeInputs ret;
    FakeOverlaps FakeLaps;
    int szSource = ini.stateProps.sourceMasses.size();
    int szSink = ini.stateProps.sinkMasses.size();
    ret.zsource = FakeLaps.generate<T>(szSource);
    if(!!!ini.matElemProps.sameOp)
      ret.zsink = FakeLaps.generate<T>(szSink);
    else
      if( (szSource == szSink) &&
	  (ini.stateProps.sourceMasses == ini.stateProps.sinkMasses))
	ret.zsink = ret.zsource;
      else
	{
	  SPLASH("incorrect specification of source and sink masses");
	  SPLASH("stateProps.source/sinkMasses conflicts with matElemProps.sameOp, exiting");
	  exit(1);
	}
    FakeSpectrum FakeSpec(ini);
    specsource = FakeSpec.generate(std::string("source"));
    specsink = FakeSpec.generate(std::string("sink"));
    
    ret.fill();

    return ret;
  }


  template<typename T>
  struct Fake3ptFunction
  {
    Fake3ptFunction(const FakeDataIni_t &ini)
    : m_ini(ini)
    {   }

    std::list<FakeDataPoint<T> > generateData(void) const;

  private:
    Fake3ptFunction(void); // hide ctor;
    FakeDataIni_t m_ini;
  };

  
  template<typename T>
  std::list<FakeDataPoint<T> > Fake3ptFuction<T>::generateData(void) const
  {
    FakeInputs inputs = getInputs(m_ini);
    
  }


}

#endif
