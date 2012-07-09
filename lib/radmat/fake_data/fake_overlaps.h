#ifndef FAKE_OVERLAPS_H_H_GUARD
#define FAKE_OVERLAPS_H_H_GUARD

#include "fake_data_ini.h"
#include "covarrying_vectors.h"
#include "semble/semble_matrix.h"
#include <complex>
#include <vector>

namespace radmat
{

  struct FakeOverlaps
  {
    FakeOverlaps(const FakeDataIni_t &ini)
    : m_ini(ini) 
    {  }

    template<typename T>
    typename std::vector<SEMBLE::SembleMatrix<T> > generate(const int) const;

  private:
    FakeDataIni_t m_ini;
    FakeOverlaps(void); // hide ctor
  };

  template<> 
  std::vector<SEMBLE::SembleMatrix<std::complex<double> > > 
  FakeOverlaps::generate<std::complex<double> >(const int dim) const;

  template<> 
  std::vector<SEMBLE::SembleMatrix<double> > 
  FakeOverlaps::generate<double>(const int dim) const;

}

#endif
