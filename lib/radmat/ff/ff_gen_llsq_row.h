#ifndef FF_GEN_LLSQ_ROW_H_H_GUARD
#define FF_GEN_LLSQ_ROW_H_H_GUARD


#include "radmat/utils/pow2assert.h"
#include "radmat/utils/tensor.h"
#include "adat/handle.h"
#include "itpp/itbase.h"
#include "io/adat_xmlio.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include <utility>
#include <sstream>
#include <string>
#include <list>


/*
  BRIEF NOTES  -- I think this may work -- the inner tex in the document is at 
  least correct in that it compiles.  Just show that q_mu j^mu = 0 since
  it is a Noether current. 

  {\color[rgb]{0.500000,0.500000,0.500000}\documentclass[10pt]{article}
  \usepackage[usenames]{color} %used for font color
  \usepackage{amssymb} %maths
  \usepackage{amsmath} %maths
  \usepackage[utf8]{inputenc} }

  \begin{document}

  \begin{align*}
  j^{\mu}(x) &\equiv \langle x, p^{\prime}, J^{\prime P^{\prime}} | \left[\bar{\Psi}\gamma^{\mu}\Psi\right] |x, p, J^P\rangle \\
  &= \langle p^{\prime}, J^{\prime P^{\prime}} | e^{i\hat{P}_{\alpha}\hat{X}^{\alpha}} \left[\bar{\Psi}\gamma^{\mu}\Psi\right] e^{-i\hat{P}_{\alpha}\hat{X}^{\alpha}}|p, J^P\rangle \\
  &= e^{i \left(p^{\prime} -p \right)_{\alpha} x^{\alpha}} \langle p^{\prime}, J^{\prime P^{\prime}} | \left[\bar{\Psi}\gamma^{\mu}\Psi\right] |p, J^P\rangle 
  \end{align*}

  $j^{\mu}$ is a Noether current so it has a zero 4-divergence, $\partial_{\mu}j^{\mu} = 0$

  \begin{align*}
  0 = \partial_{\mu}j^{\mu} &= \partial_{\mu} e^{i \left(p^{\prime} -p \right)_{\alpha} x^{\alpha}} \langle p^{\prime}, J^{\prime P^{\prime}} | \left[\bar{\Psi}\gamma^{\mu}\Psi\right] |p, J^P\rangle \\
  &\propto \left(p^{\prime} -p \right)_{\mu} \langle p^{\prime}, J^{\prime P^{\prime}} | \left[\bar{\Psi}\gamma^{\mu}\Psi\right] |p, J^P\rangle \\
  &\rightarrow q_{\mu}\langle p^{\prime}, J^{\prime P^{\prime}} | \left[\bar{\Psi}\gamma^{\mu}\Psi\right] |p, J^P\rangle = 0 \\
  &\rightarrow q_{\mu} j^{\mu} = 0 \qquad q_{\mu} \equiv \left(p^{\prime} -p \right)_{\mu}
  \end{align*}



  \end{document}

*/

using namespace ADATXML;
//using namespace ADATIO;

namespace radmat
{

  // basically hold two momenta in the tensor form
  //                     p_f               p_i
  typedef std::pair<Tensor<double,1>, Tensor<double,1> > PInv_t;
 
  
  Tensor<double,1> pPlus(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i);   // p_f + p_i
  Tensor<double,1> pMinus(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i);  // p_f - p_i


  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  //! copy a rank one tensor into an itpp vector for easy manipulation of lorentz 4-vectors
  template<typename T>
  itpp::Vec<T> toItpp(const Tensor<T,1> &v)
  {
    idx_t dim = v.getDim(0);
    itpp::Vec<T> foo(dim);
    for(idx_t i = 0; i < dim; i++)
      foo[i] = v[i];

    return foo;
  }

  template<typename T>
  Tensor<T,1> toTensor(const itpp::Vec<T> &v)
  {
    idx_t dim = v.size();
    Tensor<T,1> foo(TensorShape<1>()[dim],0.);
    for(idx_t i = 0; i < dim; i++)
      foo[i] = v[i];

    return foo;
  }
  

  //! copy a rank 2 tensor into an itpp mat 
  //  -- don't know why I bothered writing this but it may prove to be useful
  template<typename T>
  itpp::Mat<T> toItpp(const Tensor<T,2> &m)
  {
    idx_t row = m.getDim(0);
    idx_t col = m.getDim(1);
    itpp::Mat<T> foo(row,col);

    for(idx_t i = 0; i < row; i++)
      for(idx_t j = 0; j < col; j++)
	foo(i,j) = m[i][j];

    return foo;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  
  // F(Q^2)*v^{\mu}
  // this is a base class to set up list of handles to polymorphic functors to construct a general 
  // linear system
  template<typename T>
  struct ffBlockBase_t
  {
    virtual std::string ff(void) const {return std::string("unknown");}
    virtual Tensor<T,1> operator()(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i) const = 0;
    virtual ~ffBlockBase_t(void) {}
  };


  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////


  // generate some linear least squares system w/o having to know very much about the 
  // matrix element invariants. These are still at the configuration level. Just 
  // shove this into a SembleMatrix bin by bin to deal with the ensemble stats.

  // this will be polymorphic with different named classes corresponding to 
  // different quantum numbers, ie: PiPi  <-->   < 0+ | j_mu | 0+ >
  template<typename T> 
  struct ffBase_t
  {
    // save some typing
    typedef  ffBlockBase_t<T> BBType;
    typedef typename ADAT::Handle< ffBlockBase_t< T > > BBHandle_t;
    typedef std::list< BBHandle_t > ff_list;

    // this will be useful when we derive
    // a PiPi mat elem constructor will be something like  |    PiPi(void)
    //                                                     |     : ffBase_t(getPiPi()) {}
    ffBase_t(const ff_list& list)
    : m_list(list) 
    {  }

    ffBase_t& operator=(const ffBase_t &o)
    {
      if(this != &o)
	{
	  m_list = o.m_list;
	}
      return *this;
    }

    ffBase_t(const ffBase_t &o)
    : m_list(o.m_list)
    {  }
  
    // needs to be present and virtual b/c we are using pointers to derived
    virtual ~ffBase_t(void) {}

    // useful higher up
    virtual int nFacs(void) {return m_list.size();}

    // generate some tex code corresponding to the 
    //  string of stuff we think this is making
    std::string ff(void) const 
    {
      std::stringstream ss;
      typename ff_list::const_iterator it;
      for (it = m_list.begin(); it != m_list.end(); it++)
	ss << (*it)->ff() << "  ";
      return ss.str();
    }

    // generate the linear system based on the available set of kinematic factors
    // eventually wanna pass this through some type of isZero() filter?  -- do it at the ensemble level 
    // -- sometimes things will obviously be zero eg eps(z)_0,3 == 0
    virtual itpp::Mat<T> operator()(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i) const
    {
      itpp::Mat<T> ret;
      typename ff_list::const_iterator it;
      for (it = m_list.begin(); it != m_list.end(); it++)
	ret.append_col(toItpp<T>((**it)(p_f,p_i)));

      return ret;
    }

  protected:  // hide ctor
    ffBase_t(void);
    
    // data store
    ff_list m_list;
  };


  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  // save some typing
  typedef std::pair<SEMBLE::SembleVector<double> , 
		    SEMBLE::SembleVector<double> > SemblePInv_t;


SemblePInv_t makeMomInvariants(const EnsemReal &E_f, 
			       const EnsemReal &E_i,
			       const Array<int> &p_f,
			       const Array<int> &p_i,
			       const double mom_factor);  // 1/xi * 2pi /L_s -- the "unit" size

/**
   @brief generate the kinematic factor matrix for one measurement,
     
   @details ie: a 4 x n form factors matrix where the row index is the lorentz index.
   These then need to get stitched together higher up 
     
*/
template<typename T>
struct ffKinematicFactors_t
{
  // save some typing
  typedef typename ADAT::Handle<ffBase_t<T> > ffBase_h;
  typedef typename SEMBLE::SembleMatrix<T> KinematicFactorMatrix;

  // the only available constructor
  ffKinematicFactors_t(const ffBase_h &KFacGen)
  : m_KFacGen(KFacGen) 
  {  }

  ~ffKinematicFactors_t(void) {} // handle cleans itself up

  // basically generate the 4 X (n multipole) matrix of kinematic factors, 
  // the row index is the lorentz index of the lattice matrix element
  KinematicFactorMatrix genFactors(const SemblePInv_t &moms)
  {
    SEMBLE::SembleVector<double> p_f(moms.first), p_i(moms.second);
    int nfacs = m_KFacGen->nFacs();
    int nbins = p_f.getB();

    POW2_ASSERT_DEBUG( (nbins == p_i.getB()) && (nfacs > 0) );

    KinematicFactorMatrix KF(nbins,4,nfacs);
    KF.zeros();

    // scale down
    p_f.rescaleSembleDown();
    p_i.rescaleSembleDown();

    for(int bin = 0; bin < nbins; bin++)
      KF[bin] = (*m_KFacGen)(toTensor<double>(p_f[bin]),toTensor<double>(p_i[bin]));

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /*
      TO DO -- put the isZero()? type check in here before the return statement
    */
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    // scale up
    KF.rescaleSembleUp();
    
    return KF;
  }

  int nFacs(void) {return m_KFacGen->nFacs();}

  // hide ctor    
private:
  ffKinematicFactors_t(void);
  ffKinematicFactors_t(const ffKinematicFactors_t &o);
  ffKinematicFactors_t& operator=(const ffKinematicFactors_t<T> &o);

private:
  ffBase_h m_KFacGen;
};

}

#endif
