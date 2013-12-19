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

    // will bug if indicies arent correct
    template<typename L, typename R>
      typename Promote<L,R>::Type_t dot_out(const Tensor<L,1> &l, const Tensor<R,1> &r)
      {
        return contract_b(l,r); 
      }


    struct Ingredients
      : public leftSpinPTensor<1>, 
      public rightSpinPTensor<1>
    {
      Ingredients() {}

      // ingredient list
      Ingredients(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i, 
          const double mom_fac,
          const int lambda_left,
          const int lambda_right) 
      {
        // polarization tensors (up indicies)
        eps_left = this->left_p_tensor(p_f,p_i,mom_fac,lambda_left); 
        eps_right = this->right_p_tensor(p_f,p_i,mom_fac,lambda_right); 

        // momentum tensors (up indicies) 
        p_left = convertTensorUnderlyingType<std::complex<double> , double, 1>(p_f);
        p_right = convertTensorUnderlyingType<std::complex<double> , double, 1>(p_i);

        // sum and difference (up indicies)
        pplus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pPlus(p_f,p_i)); 
        pminus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pMinus(p_f,p_i)); 

        // metric tensor
        gdd = convertTensorUnderlyingType<std::complex<double>, double, 2>(g_dd()); 

        // lorentz scalars
        mm_left = dot_out(p_left,applyMetric(p_left,gdd,0)).real(); 
        mm_right = dot_out(p_right,applyMetric(p_right,gdd,0)).real(); 
        qq = dot_out(pminus,applyMetric(pminus,gdd,0)).real(); 
        Q2 = -qq; 

        // sanity
        enforce_positive(mm_left,"negative mass left");
        enforce_positive(mm_right,"negative mass right"); 
        enforce_equal(mm_right,mm_left,"m_left != m_right"); 
      }

      void print_ingredients(void)
      {
        std::cout << __func__ << ":\n"
          << "p_left " << p_left << "" 
          << "eps_left " << eps_left << ""
          << "mm_left " << mm_left << "\n" 
          << "\np_right " << p_right << ""
          << "eps_right " << eps_right << ""
          << "mm_right " << mm_right << std::endl; 
      }

      void 
        enforce(const bool &b, const std::string &msg) const
        {
          if ( !!! b ) 
          {
            std::cout << __PRETTY_FUNCTION__ << ":error " << msg << std::endl; 
            exit(1); 
          }
        }

      // 1 part in 20
      bool 
        is_equal(const double &l, const double &r) const
        {
          return (fabs( (l - r) / (l+r)) < 5e-2); 
        }

      void 
        enforce_positive(const double &d, const std::string &msg) const
        {
          enforce( d > 1e-6 , msg); 
        }

      void 
        enforce_equal(const double d1 , const double &d2, const std::string &msg) const
        {
          std::stringstream ss;
          ss << "d1=" << d1 << " d2=" << d2;
          enforce(is_equal(d1,d2) , msg + ss.str()); 
        }

      Tensor<std::complex<double> , 1> 
        init_4_tens() const
        {
          return Tensor<std::complex<double> , 1> ((TensorShape<1>())[4]); 
        }

      std::complex<double>
        eps_left_dot_p_right() const
        {
          return dot_out(eps_left , applyMetric(p_right,gdd,0)); 
        }


      std::complex<double> 
        eps_right_dot_p_left() const
        {
          return dot_out(eps_right, applyMetric(p_left,gdd,0)); 
        }

      std::complex<double> 
        eps_left_dot_eps_right() const
        {
          return dot_out(eps_left, applyMetric(eps_right,gdd,0)); 
        }

      Tensor<std::complex<double>,1> eps_left; 
      Tensor<std::complex<double>,1> eps_right; 
      Tensor<std::complex<double>,1> pplus; 
      Tensor<std::complex<double>,1> pminus; 
      Tensor<std::complex<double>,2> gdd; 
      Tensor<std::complex<double>,1> p_left; 
      Tensor<std::complex<double>,1> p_right; 
      double mm_left;
      double mm_right; 
      double qq; 
      double Q2; 
    };



    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    
    struct G1impl
      : public FormFacRotationManager<G1impl>
    {
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
          Ingredients i(p_f,p_i,mom_fac,lh,rh); 
          return  -i.pplus * i.eps_left_dot_eps_right(); 
        }
    };

    template<int lambda_left, int lambda_right>
      struct G1 : public canonicalFrameFormFactor<1,1,lambda_left,lambda_right, G1impl >
    {
      virtual ~G1() {}
    };


    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
      struct G2impl : public FormFacRotationManager<G2impl>
    {
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
            const int l, 
            const int r)  const
        {
          Ingredients i(p_f,p_i,mom_fac,l,r); 
          Tensor<std::complex<double> , 1> ret(i.init_4_tens()); 

          ret = i.eps_right * i.eps_left_dot_p_right(); 
          ret += i.eps_left * i.eps_right_dot_p_left(); 

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
      struct G3impl : public FormFacRotationManager<G3impl>
    {
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
            const int l,
            const int r)  const
        {
          Ingredients i(p_f,p_i,mom_fac,l,r); 
          return - i.pplus * ( i.eps_left_dot_p_right() * i.eps_right_dot_p_left() / (2. * i.mm_left)); 
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
