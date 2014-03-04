#ifndef LORENTZFF_CANONICAL_RHORHO_H
#define LORENTZFF_CANONICAL_RHORHO_H 


#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"
#include <exception>
#include <sstream>
#include <complex>

namespace radmat
{


  // PRD 73, 074507 (2006) 
  //

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  struct RhoRhoG1; 
  REGISTER_STRINGIFY_TYPE( RhoRhoG1 ); 


  struct RhoRhoG1
    : public FormFacRotationManager<RhoRhoG1, std::complex<double> , 1, 1>,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoG1() {}

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

        return -( val.value() * ret ); 
      }
  };

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  struct RhoRhoG2; 
  REGISTER_STRINGIFY_TYPE( RhoRhoG2 ); 

  struct RhoRhoG2 
    : public FormFacRotationManager<RhoRhoG2, std::complex<double> , 1 , 1>,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoG2() {}
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

#if 0 
        std::cout << __func__ << ": pleft = " << p_left 
          << " pright = " << p_right 
          << " epsl = " << eps_left 
          << " epsr = " << eps_right 
          << " el.pr = " << val_a.value() << "     er.pl = " << val_b.value() 
          << std::endl;  
#endif 

        ret = val_a.value() * eps_right + val_b.value() * eps_left; 

        return ret; 
      }
  };


  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  struct RhoRhoG3; 
  REGISTER_STRINGIFY_TYPE(RhoRhoG3); 

  struct RhoRhoG3 
    : public FormFacRotationManager<RhoRhoG3, std::complex<double> , 1, 1>,
    public leftSpinPTensor<1> , 
    public rightSpinPTensor<1>
  {
    typedef std::complex<double> Data_t;
    virtual ~RhoRhoG3() {}
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
        // return - val_a.value() * val_b.value() * ret;
        return  - ( (val_a.value() * val_b.value() ) / (2.*mass.value()) ) * ret; 
      }
  };

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  template<int embedl, int embedr>
    LorentzFFAbsBase_t::LorentzFFAbs_list RhoRhoGenList(void)
    {
      LorentzFFAbsBase_t::LorentzFFAbs_list retRhoRho; 
      LorentzFFAbsBase_t::BBType *g1 , *g2, *g3; 

      try
      {
        g1 = new radmat::RhoRhoG1();
        g2 = new radmat::RhoRhoG2();
        g3 = new radmat::RhoRhoG3();

        // POW2_ASSERT(g1 && g2 && g3);
        POW2_ASSERT( g1 );
        POW2_ASSERT( g2 );
        POW2_ASSERT( g3 );

        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(g1)); 
        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(g2)); 
        retRhoRho.push_back(LorentzFFAbsBase_t::BBHandle_t(g3)); 
      }
      catch(...)
      {
        POW2_ASSERT(false); 
      } 

      return retRhoRho;
    }



  template<int embedl,int embedr> struct RhoRho; 
  REGISTER_STRINGIFY_TYPE2( RhoRho<1,1> ); 


  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  template<int embedl, int embedr>
    struct RhoRho : public LorentzFFAbsBase_t
  {
    RhoRho(void)
      : LorentzFFAbsBase_t(radmat::RhoRhoGenList<embedl,embedr>())
    { }

    RhoRho& operator=(const RhoRho &o)
    {
      if (this != &o)
        LorentzFFAbsBase_t::operator=(o);

      return *this; 
    }

    RhoRho(const RhoRho &o)
      : LorentzFFAbsBase_t(o)
    {  }

    virtual ~RhoRho() {}

    virtual std::string id(void) const { return Stringify< RhoRho<embedl,embedr> >(); }
    virtual int left_spin(void) const { return embedl; }
    virtual int right_spin(void) const { return embedr; }

    private: 
    RhoRho(const LorentzFFAbsBase_t::LorentzFFAbs_list &); 
    RhoRho(const LorentzFFAbsBase_t::LorentzFFAbs_list ); 
  };


} // namespace radmat



#endif /* LORENTZFF_CANONICAL_RHORHO_H */
