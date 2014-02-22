#ifndef LORENTZFF_CANONICAL_RHOPI_H
#define LORENTZFF_CANONICAL_RHOPI_H 



#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"
#include "lorentzff_polarization_embedding.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_frame_formfactor.h"

namespace radmat
{


    struct RhoPiF1impl
      : public FormFacRotationManager<RhoPiF1impl,std::complex<double> >,
      public leftSpinPTensor<1>
    {
      typedef std::complex<double> Data_t; 

      virtual ~RhoPiF1impl() {}

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

    struct RhoPiF1;
    REGISTER_STRINGIFY_TYPE( RhoPiF1 ); 

    struct RhoPiF1
      : public canonicalFrameFormFactor<1,0,RhoPiF1impl>
    {
      virtual ~RhoPiF1() {}
      virtual std::string id() const { return Stringify<RhoPiF1>(); }
    };


    template<int embedl, int embedr>
      LorentzFFAbsBase_t::LorentzFFAbs_list RhoPiGenList(void)
      {
        LorentzFFAbsBase_t::LorentzFFAbs_list retRhoPi;
        LorentzFFAbsBase_t::BBType *blockPtr;
        blockPtr = new radmat::RhoPiF1();
        POW2_ASSERT(blockPtr);
        retRhoPi.push_back(LorentzFFAbsBase_t::BBHandle_t(blockPtr));
        return retRhoPi;
      }


  template<int embedl, int embedr> struct RhoPi;
  REGISTER_STRINGIFY_TYPE2( RhoPi<1,0> ); 



  template<int embedl, int embedr>
    struct RhoPi : public LorentzFFAbsBase_t
  {
    RhoPi(void) 
      : LorentzFFAbsBase_t(radmat::RhoPiGenList<embedl,embedr>())
    {   }

    RhoPi& operator=(const RhoPi &o)
    {
      if(this != &o)
        LorentzFFAbsBase_t::operator=(o);
      return *this; 
    }

    RhoPi(const RhoPi &o)
      : LorentzFFAbsBase_t(o)
    {  }

    virtual ~RhoPi() {}

    virtual std::string id(void) { return Stringify< RhoPi<embedl,embedr> >(); }
    virtual int left_spin(void) const { return embedl; }
    virtual int right_spin(void) const { return embedr; }

    private:
    RhoPi(const LorentzFFAbsBase_t::LorentzFFAbs_list &);
    RhoPi(const LorentzFFAbsBase_t::LorentzFFAbs_list); 

  };



} // radmat





#endif /* LORENTZFF_CANONICAL_RHOPI_H */
