#ifndef LORENTZFF_RHOPI_H
#define LORENTZFF_RHOPI_H

#include "ff_gen_llsq_row.h"
#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"


namespace radmat
{

  namespace RhoPi
  {
    template <short lambda> 
      struct F1 : public ffBlockBase_t<std::complex<double> > ,
      public leftPTensor<1,lambda>
    {

      std::string ff(void) const
      {
        std::string s;
        s += "F_1(Q^2) \\epsilon^{\\mu,\\nu,\\rho,\\sigma}\\epsilon^{*}_{\\nu}";
        s += "(p,\\lambda)p_{\\rho}^{+}p_{\\sigma}^{-}";
        return s; 
      }


      Tensor<std::complex<double> , 1>  operator()(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, const double mom_fac) const
      {

        // come up with the ingredient list
        Tensor<std::complex<double>, 1> epsilon = this->left_p_tensor(p_f,p_i,mom_fac);
        Tensor<std::complex<double>, 1> pplus, pminus;
        pplus = convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
        pminus = convertTensorUnderlyingType<std::complex<double>,double,1>( pMinus(p_f,p_i) );
        Tensor<std::complex<double>, 4>  levi = levi_civita<std::complex<double>,4>(); 
        Tensor<std::complex<double>, 2> gdd;
        gdd = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd());

        pminus = applyMetric(pminus,gdd,0); 
        pplus = applyMetric(pplus,gdd,0); 
        epsilon = applyMetric(epsilon,gdd,0); 

#if 1
        Tensor<std::complex<double> , 0> inner_prod = contract( epsilon, p_f , 0 , 0 ) ; 
        if ( std::norm ( inner_prod.value() ) >  1e-6 ) 
          std::cout << "mom dotted into polarization was " << inner_prod.value() << std::endl; 
#endif 

        return contract(
            contract(
                contract(levi,
                    pminus , 3 , 0),
                pplus , 2 , 0 ),
            epsilon , 1 , 0 );
      }

    };


    template<short lambda>
      ffBase_t<std::complex<double> >::ff_list genList(void)
      {
        ffBase_t<std::complex<double> >::ff_list retRhoPi;
        ffBase_t<std::complex<double> >::BBType *blockPtr;
        blockPtr = new radmat::RhoPi::F1<lambda>();
        POW2_ASSERT(blockPtr);
        retRhoPi.push_back(ffBase_t<std::complex<double> >::BBHandle_t(blockPtr));
        return retRhoPi;
      }


    template<short lambda>
      struct RhoPi : public ffBase_t<std::complex<double> >
    {
      RhoPi(void) 
        : ffBase_t<std::complex<double> >(radmat::RhoPi::genList<lambda>())
      {   }

      RhoPi& operator=(const RhoPi &o)
      {
        if(this != &o)
          ffBase_t<std::complex<double> >::operator=(o);
        return *this; 
      }

      RhoPi(const RhoPi &o)
        : ffBase_t<std::complex<double> >(o)
      {  }

      private:
      RhoPi(const ffBase_t<std::complex<double> >::ff_list &);
      RhoPi(const ffBase_t<std::complex<double> >::ff_list); 

    };

  } // namespace RhoPi



} // namespace anonomyous


#endif /* LORENTZFF_RHOPI_H */
