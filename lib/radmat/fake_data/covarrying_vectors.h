#ifndef COVARRYING_VECTORS_H_H_GUARD
#define COVARRYING_VECTORS_H_H_GUARD

#include "itpp/itbase.h"
#include <vector>
#include <complex>

namespace radmat
{

  struct cd
  {
    typedef std::complex<double> dc;
  };

  // return a vector/matrix of random numbers
  struct wrapRandN
  {
  public:
    template<typename T>
    itpp::Vec<T> wrap(const int length) const;

    template<typename T>
    itpp::Mat<T> wrap(const int row, const int col) const;
  };

  template<> itpp::Vec<double> wrapRandN::wrap<double>(const int) const;
  template<> itpp::Vec<cd::dc> wrapRandN::wrap<cd::dc>(const int) const;
  template<> itpp::Mat<double> wrapRandN::wrap<double>(const int , const int) const;
  template<> itpp::Mat<cd::dc> wrapRandN::wrap<cd::dc>(const int , const int) const;

  //
  //


  // calculate the covariance of the distribution 
  struct wrapCov
  {
  public: 
    template<typename T>
    itpp::Mat<T> cov(const typename std::vector<itpp::Vec<T> > &V) const;
  };

  template<> itpp::Mat<double> wrapCov::cov<double>(const std::vector<itpp::Vec<double> > &) const;
  template<> itpp::Mat<cd::dc> wrapCov::cov<cd::dc>(const std::vector<itpp::Vec<cd::dc> > &) const;
  
  //
  //

  // calculate the correlation matrix of the distribution
  struct wrapCor
  {
  public:
    template<typename T>
    itpp::Mat<T> cor(const typename std::vector<itpp::Vec<T> > & V) const;
  };

  template<> itpp::Mat<double> wrapCor::cor<double>(const std::vector<itpp::Vec<double> > &) const;
  template<> itpp::Mat<cd::dc> wrapCor::cor<cd::dc>(const std::vector<itpp::Vec<cd::dc> > &) const;

  //
  //

  // generate a covarrying distribution
  template<typename T>
  std::vector<itpp::Vec<T> > genCovarryingDist(const itpp::Vec<T> &mean,
					       const itpp::Vec<double> &variance,       // (y - ybar)* x (y - ybar) -- real!
					       const itpp::Mat<T> &correlation_matrix,  // unit normalized guy
					       const int ncfg) 
  {
    int sz = mean.size();
    
    if( (variance.size() != sz) 
	|| (correlation_matrix.rows() != sz) 
	|| (correlation_matrix.rows() != correlation_matrix.cols()) )
      __builtin_trap();
    
    itpp::Mat<T> CovarianceMatrix;

    if(correlation_matrix != itpp::hermitian_transpose(correlation_matrix) )
      {
	std::cout <<__PRETTY_FUNCTION__ << std::endl;
	std::cout << "Warning, non symmetric (hermitian) correlation matricies" 
		  << " are evil, symmetrizing (hermiterizing)." << std::endl;
	CovarianceMatrix = 0.5 * (correlation_matrix + itpp::hermitian_transpose(correlation_matrix) );
     }
    else
      CovarianceMatrix = correlation_matrix;
    
    itpp::Vec<double> d(variance);
    for(int i =0; i < sz; i++)
      d(i) = std::sqrt(d(i));

    itpp::Mat<double> D = itpp::diag(d);

    CovarianceMatrix = D * CovarianceMatrix * D;

    /*
    const itpp::Mat<T> Chalf = itpp::chol(CovarianceMatrix).H();
    std::vector<itpp::Vec<T> > cv(ncfg,mean);
    typename std::vector<itpp::Vec<T> >::iterator it;
    wrapRandN wrapper;
    */

    /* wrapRandN wrapper;
    typename std::vector<itpp::Vec<T> >::iterator it;
    std::vector<itpp::Vec<T> > cv(ncfg,mean);
    itpp::Mat<T> UU,VV;
    itpp::Vec<double> ss;
    itpp::svd(CovarianceMatrix,UU,ss,VV);
    
    for(int i = 0 ; i < sz; i++)
      ss(i) = std::sqrt(ss(i));

    const itpp::Mat<T> Chalf = UU*itpp::diag(ss);
    */

    // this "converged" the fastest
    wrapRandN wrapper;
    typename std::vector<itpp::Vec<T> >::iterator it;
    std::vector<itpp::Vec<T> > cv(ncfg,mean);
    itpp::Mat<T> V;
    itpp::Vec<double> lambda;
    itpp::eig_sym(CovarianceMatrix,lambda,V);
    
    for(int i = 0; i < sz; i++)
      lambda(i) = std::sqrt(lambda(i));

    const itpp::Mat<T> Chalf = V*itpp::diag(lambda);

    for(it = cv.begin(); it != cv.end(); it++)
      *it += Chalf * wrapper.wrap<T>(sz);

    return cv;
  }

  //
  //

  // generate a valid correlation Matrix in R^n or C^n
  struct corMat
  {
  public:
    template<typename T>
    itpp::Mat<T> genMat(const int rank) const;
  };


  template<> itpp::Mat<double> corMat::genMat<double>(const int) const;
  template<> itpp::Mat<cd::dc> corMat::genMat<cd::dc>(const int) const;

} // close namespace



#endif
