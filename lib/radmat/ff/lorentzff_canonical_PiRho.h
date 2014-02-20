#ifndef LORENTZFF_CANONICAL_PIRHO_H
#define LORENTZFF_CANONICAL_PIRHO_H 


#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_frame_formfactor.h"


namespace radmat
{

  namespace CanonicalPiRho
  {

    // actual implementation of the kinematic factor
    struct F1impl
      : public FormFacRotationManager<F1impl, std::complex<double> > , 
      public rightSpinPTensor<1>
    {
      typedef std::complex<double> Data_t; 
      virtual ~F1impl() {}

      virtual std::string
        ff_impl(void) const
        {
          std::string s; 
          s = "F_1(Q^2) \\epsilon^{\\mu,\\nu,\\rho,\\sigma}\\epsilon_{\\nu}";
          s += "(p,\\lambda)p_{\\rho}^{+}p_{\\sigma}^{-}";
          return s; 
        }

      virtual Tensor<std::complex<double> , 1> 
        impl(const Tensor<double,1> &p_f, 
            const Tensor<double,1> &p_i, 
            const double mom_fac,
            int pihel,
            int rhohel) const
        {
          // come up with the ingredient list
          Tensor<std::complex<double>, 1> epsilon = this->right_p_tensor(p_f,p_i,mom_fac,rhohel); 
          Tensor<std::complex<double>, 1> pplus, pminus;
          pplus = convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
          pminus = convertTensorUnderlyingType<std::complex<double>,double,1>( pMinus(p_f,p_i) );
          Tensor<std::complex<double>,4> levi = levi_civita<std::complex<double>,4>(); ; 
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

#if 1
          Tensor<std::complex<double> , 0> inner_prod = contract( epsilon, p_i , 0 , 0 ) ; 
          if ( std::norm ( inner_prod.value() ) >  1e-6 ) 
            std::cout << "mom dotted into polarization was " << inner_prod.value() << std::endl; 
#endif 

          return  norm * contract(
              contract(
                contract(levi,
                  pminus , 3 , 0),
                pplus , 2 , 0 ),
              epsilon , 1 , 0 );
        }
    };

    // holder class for canonical frame interface
    struct F1
      : public canonicalFrameFormFactor<0,1,F1impl>
    {
      virtual ~F1() {}
    };



    // generate a list for the PiPi constructor
    //
    //  use an embedding so we can play with subduction later
    //
    template< int embed > 
      FFAbsBase_t::FFAbs_list genList(void)
      {
        FFAbsBase_t::FFAbs_list retCanonicalPiRho;
        FFAbsBase_t::BBType *blockPtr;
        blockPtr = new F1();
        POW2_ASSERT(blockPtr);
        retCanonicalPiRho.push_back(FFAbsBase_t::BBHandle_t(blockPtr));
        return retCanonicalPiRho;
      }


    //  use an embedding so we can play with subduction later
    template<int embed>
      struct CanonicalPiRho : public FFAbsBase_t
    {
      CanonicalPiRho(void)
        : FFAbsBase_t(radmat::CanonicalPiRho::genList<embed>())
      {  }

      CanonicalPiRho& operator=(const CanonicalPiRho &o)
      {
        if(this != &o)
          FFAbsBase_t::operator=(o);
        return *this;
      }

      CanonicalPiRho(const CanonicalPiRho &o)
        : FFAbsBase_t(o)
      { }

      private:
      CanonicalPiRho(const FFAbsBase_t::FFAbs_list &);
      CanonicalPiRho(const FFAbsBase_t::FFAbs_list); 

    };

  }
}








#endif /* LORENTZFF_CANONICAL_PIRHO_H */
