// covarrying_vectors.cc -
//
// Tuesday, June 12 2012
//

#include "covarrying_vectors.h"
#include "itpp/itbase.h"
#include <complex>
#include <math.h>

namespace radmat
{
  // return a vector of random normal numbers
  template<>
  itpp::Vec<double> wrapRandN::wrap<double>(const int length) const
  {
    return itpp::randn(length);
  }

  template<>
  itpp::Vec<cd::dc> wrapRandN::wrap<cd::dc>(const int length) const
  {
    return itpp::randn_c(length);
  }

  template<>
  itpp::Mat<double> wrapRandN::wrap<double>(const int a, const int b) const
  {
    return itpp::randn(a,b);
  }

  template<>
  itpp::Mat<cd::dc> wrapRandN::wrap<cd::dc>(const int a, const int b) const 
  {
    return itpp::randn_c(a,b);
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////



  // calculate covariance of the distribution
  template<>
  itpp::Mat<double> wrapCov::cov<double>(const std::vector<itpp::Vec<double> > &V) const
  {
    int ncfg = V.size();
    int sz = V.begin()->size();
    itpp::Mat<double> C(sz,sz);
    itpp::Vec<double> mean(sz);
    
    C.zeros();
    mean.zeros();

    std::vector<itpp::Vec<double> >::const_iterator it;
    for(it = V.begin(); it != V.end(); it++)
      mean += *it;

    mean /= double(ncfg);
    
    std::vector<itpp::Vec<double> > half_var(ncfg,mean);
    for(int i = 0; i < ncfg; i++)
      half_var[i] -= V[i];

    for(int row = 0; row < sz; row++)
      for(int col = row; col < sz; col++) // symmetric
	{
	  double val = 0;
	  for(it = half_var.begin(); it != half_var.end(); it++)
	    val += (*it)[row]*(*it)[col];
	  val /= double(ncfg);
	  C(row,col) = val;
	  C(col,row) = val;
	}

    return C;
  }

  
  template<>
  itpp::Mat<cd::dc> wrapCov::cov<cd::dc>(const std::vector<itpp::Vec<cd::dc> > &V) const
  {

    int ncfg = V.size();
    int sz = V.begin()->size();
    itpp::Mat<cd::dc> C(sz,sz);
    itpp::Vec<cd::dc> mean(sz);

    C.zeros();
    mean.zeros();
    
    std::vector<itpp::Vec<cd::dc> >::const_iterator it;
    for(it = V.begin(); it != V.end(); it++)
      mean += *it;

    mean /= double(ncfg);

    std::vector<itpp::Vec<cd::dc> > half_var(ncfg,mean);
    for(int i = 0; i < ncfg; i++)
      half_var[i] -= V[i];

    for(int row = 0; row < sz; row++)
      for(int col = row; col < sz; col++) // hermitian
	{
	  cd::dc val = 0;
	  
	  for(it = half_var.begin(); it != half_var.end(); it++)
	    val += (*it)[row]*std::conj((*it)[col]);    // there is an extra dagger in this formula -- think about it
                                                        // as v cross vdag 
	  val /= double(ncfg);
	  
	  C(row,col) = val;
	  C(col,row) = std::conj(val);
	}
    return C;
  }
  


  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////


  // calculate the correlation matrix  
  template<> 
  itpp::Mat<double> wrapCor::cor<double>(const std::vector<itpp::Vec<double> > &V) const
  {
    wrapCov wrapper;
    itpp::Mat<double> cov = wrapper.cov<double>(V);
    itpp::Vec<double> d = itpp::diag(cov);
    int sz = d.size();

    for(int i = 0; i < sz; i++)
      d(i) = 1./std::sqrt(d(i));

    return itpp::diag(d) * cov * itpp::diag(d);
  }

  template<>
  itpp::Mat<cd::dc> wrapCor::cor<cd::dc>(const std::vector<itpp::Vec<cd::dc> > &V) const
  {
    wrapCov wrapper;
    itpp::Mat<cd::dc> cov = wrapper.cov<cd::dc>(V);
    itpp::Vec<cd::dc> cd = itpp::diag(cov);
    int sz = cd.size();
    itpp::Vec<double> d(sz);


    for(int i = 0; i < sz; i++)
      d(i) = 1./pow(std::norm(cd(i)),0.25);  // 4th root of the norm is the sqrt of the variance

    return itpp::diag(d) * cov * itpp::diag(d);
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////


  // holmes correlation matrix -- NB these aren't truly random since
  // the points on the n-ball aren't truly random, the diagonals have 
  // a higher probability density
  template<>
  itpp::Mat<double> corMat::genMat<double>(const int rank) const
  {
    itpp::Mat<double> foo = itpp::randn(rank,rank);

    for(int i = 0; i < rank; i++)
      foo.set_row(i,foo.get_row(i)/std::sqrt(itpp::dot(foo.get_row(i),foo.get_row(i))));

    foo = foo*itpp::hermitian_transpose(foo);
    
    return 0.5*(foo + itpp::hermitian_transpose(foo));
  };


  template<>
  itpp::Mat<cd::dc> corMat::genMat<cd::dc>(const int rank) const
  {
    itpp::Mat<cd::dc> foo = itpp::randn_c(rank,rank);
    
    for(int i =0; i < rank; i++)
      {
	double bar = std::real(itpp::elem_mult_sum(foo.get_row(i),itpp::conj(foo.get_row(i))));
	foo.set_row(i,foo.get_row(i)/std::sqrt(bar));
      }

    return foo * itpp::hermitian_transpose(foo);
  }

}
