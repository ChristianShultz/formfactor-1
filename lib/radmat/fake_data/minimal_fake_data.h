#ifndef MINIMAL_FAKE_DATA_H_H_GUARD
#define MINIMAL_FAKE_DATA_H_H_GUARD


#include "adat/handle.h"
#include "radmat/ff/ff_gen_llsq_row.h"
#include "semble/semble_matrix.h"
#include <string>
#include "itpp/itbase.h"
#include "radmat/utils/tensor.h"


// comment me someday

namespace radmat
{

  struct genMinimalFakeData_c
  {
    genMinimalFakeData_c(const std::string &matElemID)
      : m_matGen_h(FormFactorDecompositionFactoryEnv::callFactory(matElemID) )
    {    }

    itpp::Mat<std::complex<double> > operator()(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i) const
    {
      return (*m_matGen_h)(p_f,p_i);
    }

    SEMBLE::SembleMatrix<std::complex<double> > operator()(const Tensor<double,1> &p_f, 
							   const Tensor<double,1> &p_i,
							   const int ncfg,
							   const double &noiseOrder,
							   bool normal=true) const
    {
      SEMBLE::SembleMatrix<std::complex<double> > bar(ncfg,4,m_matGen_h->nFacs());
      bar = (*m_matGen_h)(p_f,p_i);
      int N = bar.getN();

      if(noiseOrder != 0.)
	for(int i = 0; i < ncfg; i++)
	  if(normal)
	    bar[i] += noiseOrder*itpp::randn_c(4,N);
	  else
	    bar[i] += noiseOrder*itpp::randn_c(4,N);

      return bar;
    }

    void swapMatElem(const std::string &newMatElemID)
    {
      m_matGen_h = FormFactorDecompositionFactoryEnv::callFactory(newMatElemID);
    }

    // hide ctor
  private:
    genMinimalFakeData_c(void);

    ADAT::Handle<ffBase_t<std::complex<double> > > m_matGen_h;
  };



} // close namespace


#endif
