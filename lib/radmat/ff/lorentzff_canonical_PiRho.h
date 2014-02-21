#ifndef LORENTZFF_CANONICAL_PIRHO_H
#define LORENTZFF_CANONICAL_PIRHO_H 


#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_frame_formfactor.h"


namespace radmat
{


  // actual implementation of the kinematic factor
  struct PiRhoF1impl
    : public FormFacRotationManager<PiRhoF1impl, std::complex<double> > , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t; 
    virtual ~PiRhoF1impl() {}

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

  struct PiRhoF1; 
  REGISTER_STRINGIFY_TYPE( PiRhoF1 ); 

  // holder class for canonical frame interface
  struct PiRhoF1
    : public canonicalFrameFormFactor<0,1,PiRhoF1impl>
  {
    virtual ~PiRhoF1() {}
    virtual std::string id() const { return Stringify<PiRhoF1>(); }
  };


  // generate a list for the PiPi constructor
  //
  //  use an embedding so we can play with subduction later
  //
  template< int embedl , int embedr  > 
    FFAbsBase_t::FFAbs_list PiRhoGenList(void)
    {
      FFAbsBase_t::FFAbs_list retCanonicalPiRho;
      FFAbsBase_t::BBType *blockPtr;
      blockPtr = new PiRhoF1();
      POW2_ASSERT(blockPtr);
      retCanonicalPiRho.push_back(FFAbsBase_t::BBHandle_t(blockPtr));
      return retCanonicalPiRho;
    }


  ////////////////////////////
  ////////////////////////////

  template<int embedl, int embedr> struct PiRho; 
  REGISTER_STRINGIFY_TYPE2( PiRho<0,1> ); 


  //  use an embedding so we can play with subduction later
  template<int embedl, int embedr>
    struct PiRho : public FFAbsBase_t
  {
    PiRho(void)
      : FFAbsBase_t(radmat::PiRhoGenList<embedl,embedr>())
    {  }

    PiRho& operator=(const PiRho &o)
    {
      if(this != &o)
        FFAbsBase_t::operator=(o);
      return *this;
    }

    PiRho(const PiRho &o)
      : FFAbsBase_t(o)
    { }

    virtual ~PiRho() {}

    virtual std::string id(void) const { return Stringify<PiRho<embedl,embedr> >(); }
    virtual int left_spin(void) const { return embedl; }
    virtual int right_spin(void) const { return embedr; }

    private:
    PiRho(const FFAbsBase_t::FFAbs_list &);
    PiRho(const FFAbsBase_t::FFAbs_list); 

  };

}








#endif /* LORENTZFF_CANONICAL_PIRHO_H */
