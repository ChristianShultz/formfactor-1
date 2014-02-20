#ifndef LORENTZFF_CANONICAL_RHOPI_H
#define LORENTZFF_CANONICAL_RHOPI_H 



#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_polarization_embedding.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_frame_formfactor.h"

namespace radmat
{

  namespace CanonicalRhoPi
  {
    struct F1impl
      : public FormFacRotationManager<F1impl,std::complex<double> >,
      public leftSpinPTensor<1>
    {
      typedef std::complex<double> Data_t; 

      virtual ~F1impl() {}

      virtual std::string ff_impl(void) const
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

          Tensor<double, 0> m_left, m_right;
          std::complex<double> norm; 
          m_left = contract(p_f,applyMetric(p_f,g_dd(),0),0,0);
          m_right = contract(p_i,applyMetric(p_i,g_dd(),0),0,0);
          norm = std::complex<double>( 2./( sqrt(m_left.value()) + sqrt(m_right.value()) ), 0.); 

     //     std::cout << __func__ << "hl " << hel << " hr " << zzero << std::endl;
     //
     //     std::cout 
     //      << __func__ << ": pleft " << p_f  
     //      << __func__ << ": pright " << p_i 
     //      << __func__ << ": pplus " << pplus 
     //      << __func__ << ": minus " << pminus 
     //      << __func__ << ": epsilon " << epsilon 
     //      << std::endl;

          return norm * contract(
              contract(
                contract(levi,
                  pminus , 3 , 0),
                pplus , 2 , 0 ),
              epsilon , 1 , 0 );
        }
    };

      struct F1
      : public canonicalFrameFormFactor<1,0,F1impl>
      {
        virtual ~F1() {}
      };


    template<int embed>
      FFAbsBase_t::FFAbs_list genList(void)
      {
        FFAbsBase_t::FFAbs_list retCanonicalRhoPi;
        FFAbsBase_t::BBType *blockPtr;
        blockPtr = new radmat::CanonicalRhoPi::F1<embed>();
        POW2_ASSERT(blockPtr);
        retCanonicalRhoPi.push_back(FFAbsBase_t::BBHandle_t(blockPtr));
        return retCanonicalRhoPi;
      }


    template<int embed>
      struct CanonicalRhoPi : public FFAbsBase_t
    {
      CanonicalRhoPi(void) 
        : FFAbsBase_t(radmat::CanonicalRhoPi::genList<embed>())
      {   }

      CanonicalRhoPi& operator=(const CanonicalRhoPi &o)
      {
        if(this != &o)
          FFAbsBase_t::operator=(o);
        return *this; 
      }

      CanonicalRhoPi(const CanonicalRhoPi &o)
        : FFAbsBase_t(o)
      {  }

      private:
      CanonicalRhoPi(const FFAbsBase_t::FFAbs_list &);
      CanonicalRhoPi(const FFAbsBase_t::FFAbs_list); 

    };

  } // CanonicalRhoPi

} // radmat



#endif /* LORENTZFF_CANONICAL_RHOPI_H */
