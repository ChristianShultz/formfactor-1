#ifndef LORENTZFF_PIRHO_H
#define LORENTZFF_PIRHO_H

#include "ff_gen_llsq_row.h"
#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"


namespace radmat
{

  namespace PiRho
  {
    template<short lambda>
      struct F1 
      : public ffBlockBase_t<std::complex<double> > , 
      public rightPTensor<1,lambda>
    {

      std::string ff(void) const
      {
        std::string s; 
        s = "F_1(Q^2) \\epsilon^{\\mu,\\nu,\\rho,\\sigma}\\epsilon_{\\nu}";
        s += "(p,\\lambda)p_{\\rho}^{+}p_{\\sigma}^{-}";
        return s; 
      }

      Tensor<std::complex<double> , 1>  operator()(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, const double mom_fac) const
      {
        // come up with the ingredient list
        Tensor<std::complex<double>, 1> epsilon = this->right_p_tensor(p_i,mom_fac,p_f); 
        Tensor<std::complex<double>, 1> pplus, pminus;
        pplus = convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
        pminus = convertTensorUnderlyingType<std::complex<double>,double,1>( pMinus(p_f,p_i) );
        Tensor<std::complex<double>, 4> levi = levi_civita<std::complex<double> , 4>(); 
        Tensor<std::complex<double>, 2> gdd;
        gdd = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd());

        // the intermediary steps.. since we are contracting w/ vectors
        // we go down by one rank at each step
        Tensor<std::complex<double>, 3> foo;
        Tensor<std::complex<double>, 2> bar;
        Tensor<std::complex<double>, 1> baz; 

        // do contractions 
        foo = contract(levi,applyMetric(pminus,gdd,0),3,0);
        bar = contract(foo,applyMetric(pplus,gdd,0),2,0);
        baz = contract(bar,applyMetric(epsilon,gdd,0),1,0);

#if 1
        // sanity
        Tensor<std::complex<double>, 0> inner_prod; 
        inner_prod = contract(epsilon,applyMetric(p_i,gdd,0),0,0); 
        std::stringstream ss; 
        ss << "momentum dotted into polarization was " << inner_prod.value() << std::endl ;
        if(std::norm(inner_prod.value()) > 0.000001 ) 
          std::cout << ss.str() << std::endl;  

#endif
        // the kinematic factor carries one lorentz index
        return baz;
      }
    };



    // generate a list for the PiPi constructor
    template<short lambda>
      ffBase_t<std::complex<double> >::ff_list genList(void)
      {
        ffBase_t<std::complex<double> >::ff_list retPiRho;
        ffBase_t<std::complex<double> >::BBType *blockPtr;
        blockPtr = new F1<lambda>();
        POW2_ASSERT(blockPtr);
        retPiRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(blockPtr));
        return retPiRho;
      }

    template<short lambda>
      struct PiRho : public ffBase_t<std::complex<double> >
    {
      PiRho(void)
        : ffBase_t<std::complex<double> >(radmat::PiRho::genList<lambda>())
      {  }

      PiRho& operator=(const PiRho &o)
      {
        if(this != &o)
          ffBase_t<std::complex<double> >::operator=(o);
        return *this;
      }

      PiRho(const PiRho &o)
        : ffBase_t<std::complex<double> >(o)
      { }

      private:
      PiRho(const ffBase_t<std::complex<double> >::ff_list &);
      PiRho(const ffBase_t<std::complex<double> >::ff_list); 

    };

  }
}







#endif /* LORENTZFF_PIRHO_H */
