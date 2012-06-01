#ifndef COMMON_H_H_GUARD
#define COMMON_H_H_GUARD


#include "radmat/utils/pow2assert.h"
#include "radmat/utils/tensor.h"
#include "adat/handle.h"
#include "itpp/itbase.h"
#include <utility>
#include <sstream>
#include <string>
#include <list>


/**
   @file common.h common definitions to lorentz decomposition ffs and multipole ffs
   @brief contains some common definitions 
*/


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

namespace radmat
{

  // basically hold two momenta in the tensor form
  //                     p_f               p_i
  typedef std::pair<Tensor<double,1>, Tensor<double,1> > PInv_t;
 
  
  Tensor<double,1> pPlus(const PInv_t &moms);   // p_f + p_i
  Tensor<double,1> pMinus(const PInv_t &moms);  // p_f - p_i


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
    virtual Tensor<T,1> operator()(const PInv_t &moms) const = 0;
    virtual ~ffBlockBase_t(void)=0;
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
    typedef typename ADAT::Handle< ffBlockBase_t< T > > BBHandle_t;
    typedef std::list< BBHandle_t > ff_list;

    // this will be useful when we derive
    // a PiPi mat elem constructor will be something like  |    PiPi(void)
    //                                                     |     : ffBase_t(getPiPi()) {}
    ffBase_t(const ff_list& list)
    : m_list(list) 
    {  }

  
    // needs to be present and virtual b/c we are using pointers to derived
    virtual ~ffBase_t(void) {}


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

    // generate the linear system base on the available set of kinematic factors
    // eventually wanna pass this through some type of isZero() filter?  
    // -- sometimes things will obviously be zero eg eps(z)_0,3 == 0
    virtual itpp::Mat<T> operator()(const PInv_t &moms) const
    {
      itpp::Mat<T> ret;
      typename ff_list::const_iterator it;
      for (it = m_list.begin(); it != m_list.end(); it++)
	ret.append_col(toItpp<T>((**it)(moms)));

      return ret;
    }

  protected:  // hide ctor
    ffBase_t(void);
    ffBase_t& operator=(const ffBase_t &);
    ffBase_t(const ffBase_t &);
    
    // data store
    ff_list m_list;
  };


  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  

  /*
    There needs to be something here that takes an ensemble of pairs of energies
    and two lattice momentum directions and spits back out an ensemble of 
    kinematic factors. 
  */

  /*
    ROTATIONS ? -- could we set up some driver to rip things out of adat and try to get the phases 
    right by starting with a meson op and adding some extra bit of phase?
  */

}

#endif
