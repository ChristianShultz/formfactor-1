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
      struct F1 : public ffBlockBase_t<std::complex<double> > , public embedHelicityPolarizationTensor<1,lambda>
    {

      std::string ff(void) const
      {
        return std::string("F_1(Q^2) \\epsilon^{\\mu,\\nu,\\rho,\\sigma}\\epsilon^{*}_{\\nu}(p,\\lambda)p_{\\rho}^{+}p_{\\sigma}^{-}");
      }


      Tensor<std::complex<double> , 1>  operator()(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, const double mom_fac) const
      {

        // come up with the ingredient list
        Tensor<std::complex<double>, 1> epsilon = this->conjugate((this->ptensor(p_f,mom_fac))); // nb final state conjugation 
        //   Tensor<std::complex<double>, 1 > epsilon = this->ptensor(p_f,mom_fac); 
        Tensor<std::complex<double>, 1> pplus, pminus;
        pplus = convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
        pminus = convertTensorUnderlyingType<std::complex<double>,double,1>( pMinus(p_f,p_i) );
        Tensor<std::complex<double>, 4> levi = levi_civita<std::complex<double> , 4>(); 
        Tensor<std::complex<double>, 2> gdd;
        gdd = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd());

        // the intermediary steps.. since we are contracting w/ vectors we go down by one
        // rank at each step.. duh
        Tensor<std::complex<double>, 3> foo;
        Tensor<std::complex<double>, 2> bar;
        Tensor<std::complex<double>, 1> baz; 

        pminus = applyMetric(pminus,gdd,0); 
        pplus = applyMetric(pplus,gdd,0); 
        epsilon = applyMetric(epsilon,gdd,0); 

#if 0
        std::cout << "pplus_mu = " << pplus << std::endl;
        std::cout << "pminus_mu = " << pminus << std::endl;
        std::cout << "eps_mu = " << epsilon << std::endl;


        std::cout << "levi non-zero " << std::endl;
        Tensor<std::complex<double>, 1> tt((TensorShape<1>())[4],0.); 


        for(int i = 0; i < 4; ++i)
          for(int j = 0; j < 4; ++j)
            for(int k = 0; k < 4; ++k)
              for(int l = 0; l < 4; ++l)
                if(levi[i][j][k][l] != 0.)
                {
                  tt[i] += levi[i][j][k][l]*epsilon[j]*pplus[k]*pminus[l];
                 // std::cout << i << " " << j << " " << k << " " << l << "   " << levi[i][j][k][l] << std::endl;
                }
        std::cout << __func__ << " the answer should be " << std::endl;
        std::cout << tt << std::endl; 
#endif 

        // do contractions 
        foo = contract(levi,pminus,3,0);
        bar = contract(foo,pplus,2,0);
        baz = contract(bar,epsilon,1,0);

        // the kinematic factor carries one lorentz index
        return baz;
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
