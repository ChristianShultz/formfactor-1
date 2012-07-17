#ifndef FAKE_SPECTRUM_H_H_GUARD
#define FAKE_SPECTRUM_H_H_GUARD

#include "fake_data_ini.h"
#include "covarrying_vectors.h"
#include "semble/semble_vector.h"
#include "radmat/utils/type_computations.h"
#include <vector>


namespace radmat
{

  struct FakeSpectrum
  {
    FakeSpectrum(const FakeDataIni_t &ini)
    : m_ini(ini)
    {  }

    // vector index is time, semble index is state elem
    std::vector<SEMBLE::SembleVector<double> > generate(const std::string &source_or_sink) const;

  private:
    FakeDataIni_t m_ini;
    FakeSpectrum(void); // hide ctor
  };  

}

#endif
