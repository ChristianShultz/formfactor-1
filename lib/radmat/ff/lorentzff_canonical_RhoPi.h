#ifndef LORENTZFF_CANONICAL_RHOPI_H
#define LORENTZFF_CANONICAL_RHOPI_H 



#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_polarization_embedding.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"

namespace radmat
{

  namespace CanonicalRhoPi
  {
    struct F1impl
      : public FormFacRotationManager<F1impl>
    {
      virtual std::string ff(void) const
      {
        std::string s;
        s += "F_1(Q^2) \\epsilon^{\\mu,\\nu,\\rho,\\sigma}\\epsilon^{*}_{\\nu}";
        s += "(p,\\lambda)p_{\\rho}^{+}p_{\\sigma}^{-}";
        return s; 
      }

      virtual Tensor<std::complex<double> , 1>  
        impl(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i,
          const double mom_fac, 
          const int hel, 
          const int zzero) const
      {

        // come up with the ingredient list
        Tensor<std::complex<double>, 1> epsilon = this->left_p_tensor(p_f,p_i,mom_fac,hel);
        Tensor<std::complex<double>, 1> pplus, pminus;
        pplus = convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
        pminus = convertTensorUnderlyingType<std::complex<double>,double,1>( pMinus(p_f,p_i) );
        Tensor<std::complex<double>, 4>  levi = levi_civita<std::complex<double>,4>(); 
        Tensor<std::complex<double>, 2> gdd;
        gdd = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd());

        pminus = applyMetric(pminus,gdd,0); 
        pplus = applyMetric(pplus,gdd,0); 
        epsilon = applyMetric(epsilon,gdd,0); 

        return contract(
            contract(
              contract(levi,
                pminus , 3 , 0),
              pplus , 2 , 0 ),
            epsilon , 1 , 0 );
      }
    }

    template<int lambda_l>
      struct F1
      : public embedRotaionManager<1,1,lambda_l,0,F1impl>
      {
        virtual ~F1() {}
      };

    };

    template<short lambda>
      ffBase_t<std::complex<double> >::ff_list genList(void)
      {
        ffBase_t<std::complex<double> >::ff_list retCanonicalRhoPi;
        ffBase_t<std::complex<double> >::BBType *blockPtr;
        blockPtr = new radmat::CanonicalRhoPi::F1<lambda>();
        POW2_ASSERT(blockPtr);
        retCanonicalRhoPi.push_back(ffBase_t<std::complex<double> >::BBHandle_t(blockPtr));
        return retCanonicalRhoPi;
      }


    template<short lambda>
      struct CanonicalRhoPi : public ffBase_t<std::complex<double> >
    {
      CanonicalRhoPi(void) 
        : ffBase_t<std::complex<double> >(radmat::CanonicalRhoPi::genList<lambda>())
      {   }

      CanonicalRhoPi& operator=(const CanonicalRhoPi &o)
      {
        if(this != &o)
          ffBase_t<std::complex<double> >::operator=(o);
        return *this; 
      }

      CanonicalRhoPi(const CanonicalRhoPi &o)
        : ffBase_t<std::complex<double> >(o)
      {  }

      private:
      CanonicalRhoPi(const ffBase_t<std::complex<double> >::ff_list &);
      CanonicalRhoPi(const ffBase_t<std::complex<double> >::ff_list); 

    };

  } // CanonicalRhoPi

} // radmat



#endif /* LORENTZFF_CANONICAL_RHOPI_H */
