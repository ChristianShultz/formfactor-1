#ifndef COVARRYING_VECTORS_H_H_GUARD
#define COVARRYING_VECTORS_H_H_GUARD

#include "itpp/itbase.h"
#include "semble/semble_vector.h"
#include "radmat/utils/pow2assert.h"
#include <complex>

namespace radmat
{

  // wrap the itpp random deviate call into a template
  template<typename T> 
  itpp::Vec<T> wrapRandn(const int length)
  {
    return itpp::Vec<T>(length);
  }

  template<>
  itpp::Vec<double> wrapRandn(const int length)
  {
    return itpp::randn(length);
  }

  template<>
  itpp::Vec<std::complex<double> > wrapRandn(const int length)
  {
    return itpp::randn_c(length);
  }

  
  // generate an ensemble with mean mean and whose correlation matrix (the unit normalized covariance matrix) 
  // approaches the one passed as ncfgs runs off to inf
  template<typename T>
  SEMBLE::SembleVector<T> genCovVec(const itpp::Vec<T> &mean, 
				    const itpp::Vec<T> &correlation_matrix, 
				    const int ncfg)
  {
    int dim = mean.size();
    POW2_ASSERT((dim == correlation_matrix.rows()) && 
		(correlation_matrix.rows() == correlation_matrix.cols()));

    SEMBLE::SembleVector<T> V(ncfg,mean.size());
    V = mean;

    itpp::Mat<T> U,V,C(correlation_matrix),rootCdag;
    itpp::Vec<double> s, rootS;

    if(correlation_matrix != itpp::hermitian_transpose(correlation_matrix) )
      {
	std::cout << "Warning, non-symmetric (hermitian) correlation matrix is evil, symmeterizing" << endl;
	C = 0.5*(C + itpp::hermitian_transpose(C));
      }

    itpp::svd(C,U,s,V);
    rootS = s;
    for(int i = 0; i < dim; i++)
      rootS[i] = sqrt(s);

    rootCdag = U*rootS;

    for(int i = 0; i < ncfg; i++)
      V[bin] += itpp::conj(rootCdag * wrapRandn<T>(dim) );
    
    return V;
  }
  


} // close radmat 

#endif
