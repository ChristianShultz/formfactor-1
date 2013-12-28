#ifndef LORENTZFF_CANONICAL_RHORHO_H
#define LORENTZFF_CANONICAL_RHORHO_H 


#include "lorentzff_canonical_frame_formfactor.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include <exception>
#include <sstream>

namespace radmat
{

  namespace RhoRho
  {

    // PRD 73, 074507 (2006) 
    //

    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////

    struct G1impl
      : public FormFacRotationManager<G1impl, std::complex<double> >,
      public leftSpinPTensor<1> , 
      public rightSpinPTensor<1>
    {
      typedef std::complex<double> Data_t;
      virtual ~G1impl() {}

      virtual std::string
        ff_impl() const
        {
          return "-G_1(Q^2)(p_f + p_i)^\\mu \\epsilon^*_\\nu(p_f,\\lambda_f)\\epsilon^\\nu(p_i,\\lambda_i) \\\\";
        }

      virtual Tensor<std::complex<double> , 1> 
        impl(const Tensor<double,1> &p_f, 
            const Tensor<double,1> &p_i, 
            const double mom_fac,
            const int lh,
            const int rh)  const
        {
          Tensor<std::complex<double> , 1> ret( (TensorShape<1>())[4], std::complex<double>(0.,0.)); 
          Tensor<std::complex<double> , 1> eps_left, eps_right; 
          eps_left = left_p_tensor(p_f,p_i,mom_fac,lh);
          eps_right = right_p_tensor(p_f,p_i,mom_fac,rh);
          Tensor<std::complex<double>,0> val; 
          Tensor<std::complex<double>,2> metric; 
          metric = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd()); 
          val = contract(eps_left,applyMetric(eps_right,metric,0),0,0); 

          ret = convertTensorUnderlyingType<std::complex<double>, double,1>(p_f + p_i); 
        
          return (- val.value() * ret ); 
        }
    };

    template<int lambda_left, int lambda_right>
      struct G1 : public canonicalFrameFormFactor<1,1,lambda_left,lambda_right, G1impl >
    {
      virtual ~G1() {}
    };


    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    struct G2impl : public FormFacRotationManager<G2impl, std::complex<double> >,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
    {
      typedef std::complex<double> Data_t;
      virtual ~G2impl() {}
      virtual std::string 
        ff_impl() const
        {
          std::string s = "G_2(Q^2)\\left[ \\epsilon^\\mu(p_i,\\lambda_i)\\epsilon^*_\\nu(p_f,\\lambda_f)p_i^\\nu";
          s += "+ \\epsilon^{*\\mu}(p_f,\\lambda_f)\\epsilon_\\nu(p_i,\\lambda_i)p_f^\\nu \\right] \\\\";
          return s;
        }

      virtual Tensor<std::complex<double> , 1> 
        impl(const Tensor<double,1> &p_f, 
            const Tensor<double,1> &p_i, 
            const double mom_fac,
            const int lh, 
            const int rh)  const
        {
          Tensor<std::complex<double> , 1> ret( (TensorShape<1>())[4], std::complex<double>(0.,0.)); 
          Tensor<std::complex<double> , 1> p_left,p_right; 
          Tensor<std::complex<double> , 1> eps_left, eps_right; 
          Tensor<std::complex<double> , 0> val_a, val_b; 
          Tensor<std::complex<double> , 2> metric; 
          metric = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd()); 
          eps_left = left_p_tensor(p_f,p_i,mom_fac,lh);
          eps_right = right_p_tensor(p_f,p_i,mom_fac,rh);
          p_left = convertTensorUnderlyingType< Data_t , double , 1>(p_f); 
          p_right = convertTensorUnderlyingType< Data_t , double , 1>(p_i); 
          val_a = contract(eps_left, applyMetric(p_right,metric,0),0,0);  
          val_b = contract(eps_right, applyMetric(p_left,metric,0),0,0);  

          ret = val_a.value() * eps_right + val_b.value() * eps_left; 

          return ret; 
        }
    };


    template<int lambda_left, int lambda_right>
      struct G2 : public canonicalFrameFormFactor<1,1,lambda_left,lambda_right, G2impl >
    {
      virtual ~G2() {}
    };

    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    struct G3impl : public FormFacRotationManager<G3impl, std::complex<double> >,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
    {
      typedef std::complex<double> Data_t;
      virtual ~G3impl() {}
      virtual std::string
        ff_impl() const
        {
          std::string s =  "-\\frac{G_3(Q^2)}{2m_v^2}(p_f + p_i)^\\mu";
          s += " \\epsilon^*_\\nu(p_f,\\lambda_f)p_i^\\nu";
          s += " \\epsilon_\\alpha(p_i,\\lambda_f)p_f^\\alpha \\\\ ";
          return s; 
        }
      virtual Tensor<std::complex<double> , 1> 
        impl(const Tensor<double,1> &p_f, 
            const Tensor<double,1> &p_i, 
            const double mom_fac,
            const int lh,
            const int rh)  const
        {
          Tensor<std::complex<double> , 1> ret( (TensorShape<1>())[4], std::complex<double>(0.,0.)); 
          Tensor<std::complex<double> , 1> p_left,p_right; 
          Tensor<std::complex<double> , 1> eps_left, eps_right; 
          Tensor<std::complex<double> , 0> val_a, val_b, mass; 
          Tensor<std::complex<double> , 2> metric; 
          metric = convertTensorUnderlyingType<std::complex<double>,double,2>(g_dd()); 
          eps_left = left_p_tensor(p_f,p_i,mom_fac,lh);
          eps_right = right_p_tensor(p_f,p_i,mom_fac,rh);
          p_left = convertTensorUnderlyingType< Data_t , double , 1>(p_f); 
          p_right = convertTensorUnderlyingType< Data_t , double , 1>(p_i); 
          val_a = contract(eps_left, applyMetric(p_right,metric,0),0,0);  
          val_b = contract(eps_right, applyMetric(p_left,metric,0),0,0);  
          mass = contract(p_right, applyMetric(p_right,metric,0),0,0); 
    
          ret = convertTensorUnderlyingType<std::complex<double>, double,1>(p_f + p_i); 
          return  - (val_a.value() * val_b.value() / 2.*mass.value()) * ret; 
        }
    };

    template<int lambda_left, int lambda_right>
      struct G3 : public canonicalFrameFormFactor<1,1,lambda_left,lambda_right,G3impl>
    {
      virtual ~G3() {}
    };


    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    template<int lambda_left, int lambda_right> 
      ffBase_t<std::complex<double> >::ff_list genList(void)
      {
        ffBase_t<std::complex<double> >::ff_list retRhoRho; 
        ffBase_t<std::complex<double> >::BBType *g1 , *g2, *g3; 

        try
        {
          g1 = new radmat::RhoRho::G1<lambda_left,lambda_right>();
          g2 = new radmat::RhoRho::G2<lambda_left,lambda_right>();
          g3 = new radmat::RhoRho::G3<lambda_left,lambda_right>();

          // POW2_ASSERT(g1 && g2 && g3);
          POW2_ASSERT( g1 );
          POW2_ASSERT( g2 );
          POW2_ASSERT( g3 );

          retRhoRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(g1)); 
          retRhoRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(g2)); 
          retRhoRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(g3)); 
        }
        catch(...)
        {
          POW2_ASSERT(false); 
        } 

        return retRhoRho;
      }


    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    template<int lambda_left, int lambda_right>
      struct RhoRho : public ffBase_t<std::complex<double> >
    {
      RhoRho(void)
        : ffBase_t<std::complex<double> >(radmat::RhoRho::genList<lambda_left,lambda_right>())
      { }

      RhoRho& operator=(const RhoRho &o)
      {
        if (this != &o)
          ffBase_t<std::complex<double> >::operator=(o);

        return *this; 
      }

      RhoRho(const RhoRho &o)
        : ffBase_t<std::complex<double> >(o)
      {  }

      private: 
      RhoRho(const ffBase_t<std::complex<double> >::ff_list &); 
      RhoRho(const ffBase_t<std::complex<double> >::ff_list ); 
    };

  } // namespace RhoRho


} // namespace radmat



#endif /* LORENTZFF_CANONICAL_RHORHO_H */
