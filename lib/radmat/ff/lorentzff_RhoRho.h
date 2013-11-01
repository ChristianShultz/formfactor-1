#ifndef LORENTZFF_RHORHO_H
#define LORENTZFF_RHORHO_H 

#include "ff_gen_llsq_row.h"
#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "radmat/utils/pow2assert.h"
#include <exception>

namespace radmat
{

  namespace RhoRho
  {

    // arXiv:0902.2241v1 
    //
    //  F1 -> M1
    //
    //  F2 -> E2
    //
    //  F3 -> C0
    //
    //  F4 -> C2 
    //
    //


    namespace util
    {
      Tensor<std::complex<double> > PiMu(const Tensor<double,1>, &p_f, const Tensor<double,1> &p_i)
    }

    template<short lambda_left, short lambda_right>
      struct F1 : public ffBlockBase_t<std::complex<double> > , 
        public leftPTensor<1,lambda_left>, 
        public rightPTensor<1,lambda_right>
    {
      Tensor<std::complex<double> , 1> operator()(const Tensor<double,1> &p_f, 
          const Tensor<double> &p_i, 
          const double mom_fac)  const
      {
        // ingredient list
        Tensor<std::complex<double> , 1> eps_left, eps_right, pplus, pminus, gdd;
        eps_left = this->left_p_tensor(p_f,mom_fac); 
        eps_right = this->right_p_tensor(p_i,mom_fac); 
        pplus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pPlus(p_f,p_i)); 
        pminus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pMinus(p_f,p_i)); 
        gdd = convertTensorUnderlyingType<std::complex<double>, double, 1>(g_dd()); 
      }
    };


    template<short lambda_left, short lambda_right>
      struct F2 : public ffBlockBase_t<std::complex<double> > , 
        public leftPTensor<1,lambda_left>, 
        public rightPTensor<1,lambda_right>
    {
      Tensor<std::complex<double> , 1> operator()(const Tensor<double,1> &p_f, 
          const Tensor<double> &p_i, 
          const double mom_fac)  const
      {
        // ingredient list
        Tensor<std::complex<double> , 1> eps_left, eps_right, pplus, pminus, gdd;
        eps_left = this->left_p_tensor(p_f,mom_fac); 
        eps_right = this->right_p_tensor(p_i,mom_fac); 
        pplus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pPlus(p_f,p_i)); 
        pminus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pMinus(p_f,p_i)); 
        gdd = convertTensorUnderlyingType<std::complex<double>, double, 1>(g_dd()); 
      }
    };

    template<short lambda_left, short lambda_right>
      struct F3 : public ffBlockBase_t<std::complex<double> > , 
        public leftPTensor<1,lambda_left>, 
        public rightPTensor<1,lambda_right>
    {
      Tensor<std::complex<double> , 1> operator()(const Tensor<double,1> &p_f, 
          const Tensor<double> &p_i, 
          const double mom_fac)  const
      {
        // ingredient list
        Tensor<std::complex<double> , 1> eps_left, eps_right, pplus, pminus, gdd;
        eps_left = this->left_p_tensor(p_f,mom_fac); 
        eps_right = this->right_p_tensor(p_i,mom_fac); 
        pplus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pPlus(p_f,p_i)); 
        pminus = convertTensorUnderlyingType<std::complex<double>, double, 1>(pMinus(p_f,p_i)); 
        gdd = convertTensorUnderlyingType<std::complex<double>, double, 1>(g_dd()); 
      }
    };


    template<short lambda_left, short lambda_right> 
      ffBase_t<std::complex<double> >::ff_list genList(void)
      {
        ffBase_t<std::complex<double> >::ff_list retRhoRho; 
        ffBase_t<std::complex<double> >::BBType *f1 , *f2, *f3, *f4; 
       
        try
        {
          f1 = new radmat::RhoRho::F1<lambda_left,lambda_right>();
          f2 = new radmat::RhoRho::F2<lambda_left,lambda_right>();
          f3 = new radmat::RhoRho::F3<lambda_left,lambda_right>();
          f4 = new radmat::RhoRho::F3<lambda_left,lambda_right>();

          // POW2_ASSERT(f1 && f2 && f3);
          POW2_ASSERT( f1 );
          POW2_ASSERT( f2 );
          POW2_ASSERT( f3 );

          retRhoRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(f1)); 
          retRhoRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(f2)); 
          retRhoRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(f3)); 
          retRhoRho.push_back(ffBase_t<std::complex<double> >::BBHandle_t(f4)); 
        }
        catch(...)
        {
          POW2_ASSERT(false); 
        } 
      
        return retRhoRho;
      }


    template<short lambda_left, short lambda_right>
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






#endif /* LORENTZFF_RHORHO_H */
