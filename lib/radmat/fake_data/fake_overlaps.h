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

    template<typename T>  // vector index is time, matrix index is state,op
    typename std::vector<SEMBLE::SembleMatrix<T> > generate(const int, const std::string &source_or_sink) const; 

  private:
    FakeDataIni_t m_ini;
    FakeOverlaps(void); // hide ctor
  };

  template<> 
  std::vector<SEMBLE::SembleMatrix<std::complex<double> > > 
  FakeOverlaps::generate<std::complex<double> >(const int dim, const std::string &source_or_sink) const;

  template<> 
  std::vector<SEMBLE::SembleMatrix<double> > 
  FakeOverlaps::generate<double>(const int dim, const std::string &source_or_sink) const;

  template<typename T>
  std::vector<SEMBLE::SembleVector<T> > pullRow(const std::vector<SEMBLE::SembleMatrix<T> > &stuff, const int row)
  {
    std::vector<SEMBLE::SembleVector<T> > foo;
    typename std::vector<SEMBLE::SembleMatrix<T> >::const_iterator it;
    for(it = stuff.begin(); it != stuff.end(); it++)
      foo.push_back(it->getRow(row));
    return foo;
  }

}

#endif
