// fake_overlaps.cc -
//
// Thursday, June 28 2012
//

#include "fake_overlaps.h"
#include "fake_data_ini.h"
#include "covarrying_vectors.h"
#include "itpp/itbase.h"
#include <complex>
#include <vector>
#include "semble/semble_matrix.h"
#include "radmat/utils/splash.h"

namespace  // a few useful functions that don't need to exist elsewhere
{
  itpp::Mat<double> expSkewMat(const int dim, const double tol=1e-9)
  {

    //generate a random skew symmetric matrix
    itpp::Mat<double> skew(dim,dim);
    skew.zeros();

    for(int i = 0; i < dim; ++i)
      for(int j = i+1; j < dim; ++j)
	{
	  skew(i,j) = itpp::randn();
	  skew(j,i) = -skew(i,j);
	}
  
    //perform a schur decomposition
    itpp::Mat<double> U,S;
    itpp::schur(skew,U,S);
  
    //drop all the elements comparable to zero left over from the decomposition
    S = itpp::round_to_zero(S,tol);

    itpp::Mat<double> exp_skew(dim,dim);
    int bound;
    if(dim % 2 == 0)
      bound = dim;
    else
      bound = dim -1;

    exp_skew.zeros();

    for(int block = 0; block < exp_skew.rows(); block+=2)
      {
	itpp::Mat<double> dum(2,2);
	double sigma = S(block,block+1);
	dum(0,0) = cos(sigma);
	dum(0,1) = -sin(sigma);
	dum(1,0) = -dum(0,1);
	dum(1,1) = dum(0,0);
	exp_skew.set_submatrix(block,block,dum);
      }

    //cant forget the one in the case of odd dimension
    if(dim % 2 != 0)
      exp_skew.set(dim-1,dim-1,1.);

    return U*exp_skew*U.T();
  }

  // exp(i*M)
  itpp::Mat<std::complex<double> > expHermMat(const int dim)
  {
    itpp::Mat<std::complex<double> > foo = itpp::randn_c(dim,dim);
    foo += itpp::hermitian_transpose(foo);

    itpp::Mat<std::complex<double> > V;
    itpp::Vec<double> s;
    itpp::Vec<std::complex<double> > exp_s(dim);

    itpp::eig_sym(foo,s,V);
  
    for(int i =0; i < dim; i++)
      exp_s(i) = std::complex<double>(cos(s(i)),sin(s(i)));

    return V * itpp::diag(exp_s) * itpp::hermitian_transpose(V);
  }

  template<typename T>
  itpp::Mat<T> unitary(const int dim);

  template<>
  itpp::Mat<double> unitary(const int dim)
  {
    return expSkewMat(dim);
  }

  template<>
  itpp::Mat<std::complex<double> > unitary(const int dim)
  {
    return expHermMat(dim);
  }

  template<typename T>
  itpp::Mat<T> genZ(const int dim, const std::string &s)
  {
    if(s == std::string("unitary") )
      return unitary<T>(dim);
    SPLASH("unknown overlap generator type");
    exit(1);
  }

  // hide the real work here so that the template in the headder is only instantiated for
  // the types we care about and won't compile otherwise
  template<typename T> 
  struct genLap
  {
    typedef typename std::vector<SEMBLE::SembleMatrix<T> > Zt_t;
    typedef typename std::vector<itpp::Vec<T> > v_itpp_t;

    Zt_t gen(const int dim, const radmat::FakeDataIni_t &ini) const
    {
      const int ncfg = ini.dataProps.ncfg;
      const bool uCov = ini.matElemProps.updateCovariance;
      const bool uVar = ini.matElemProps.updateVariance;
      const double varO = ini.matElemProps.varianceOrder;
      const int Lt = abs(ini.timeProps.tsink - ini.timeProps.tsource);
      
      typename Zt_t::value_type Zero(ncfg,dim,dim);
      Zero.zeros();
      Zt_t Zt(Lt,Zero);
      itpp::Mat<T> z = genZ<T>(dim,ini.stateProps.overlapGenerator);
      itpp::Vec<T> mean(Lt);
      mean.zeros();
      itpp::Vec<double> var = varO*itpp::randn(Lt);
      radmat::corMat genCor;
      itpp::Mat<T> cor = genCor.genMat<T>(Lt);
      v_itpp_t work;
     
      for(int row = 0; row < dim; row ++)
	for(int col = 0; col < dim; col++)
	  {
	    if(uCov)
	      cor = genCor.genMat<T>(Lt);
	    if(uVar)
	      var = varO*itpp::randn(Lt);
	    
	    mean = z(row,col);
	    work = radmat::genCovarryingDist(mean,var,cor,ncfg);

	    for(int t = 0; t < Lt; t++)
	      for(int bin = 0; bin < ncfg; bin++)
		Zt[t].setElement(bin,row,col,work[bin][t]);
		
	  }

      return Zt;
    }

  };

} // namespace 




namespace radmat
{

  template<> 
  std::vector<SEMBLE::SembleMatrix<std::complex<double> > > 
  FakeOverlaps::generate<std::complex<double> >(const int dim) const
  {
    genLap<std::complex<double> > gen;
    return gen.gen(dim,m_ini);
  }
 
  template<> std::vector<SEMBLE::SembleMatrix<double> > FakeOverlaps::generate<double>(const int dim) const
  {
    genLap<double> gen;
    return gen.gen(dim,m_ini);
  }



}
