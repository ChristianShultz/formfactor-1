// lorentzff_PiPiStar.cc -
//
// Tuesday, July 17 2012
//

#include"lorentzff_PiPiStar.h"
#include "radmat/utils/pow2assert.h"

namespace radmat
{

  namespace PiPiStar
  {

    std::string F1::ff(void) const
    {
      return std::string("F_1(Q^2)\\left( p_+^{\\mu}\frac{Q^2}{m_{\\pi*}^2 -m_{\\pi}^2} + p_-\\right)");
    }

    Tensor<std::complex<double> , 1> F1::operator()(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i) const
    {
      Tensor<std::complex<double>,1> pp,pm;
      pp = convertTensorUnderlyingType<std::complex<double>,double,1>(pPlus(p_f,p_i));
      pm = convertTensorUnderlyingType<std::complex<double>,double,1>(pMinus(p_f,p_i));

     // double num = (p_f - p_i) * (g_dd() * (p_f - p_i));
     // double denom = (p_f - p_i) * (g_dd() * (p_f + p_i));

      double num = value(contract(p_f-p_i,p_f-p_i,g_dd(),0,0));
      double denom = value(contract(p_f-p_i,p_f+p_i,g_dd(),0,0));

      POW2_ASSERT_DEBUG(denom != 0.);

      return ( (num/denom)*pp +  pm);
    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ffBase_t<std::complex<double> >::ff_list genList(void)
    {
      ffBase_t<std::complex<double> >::ff_list retPiPi;
      ffBase_t<std::complex<double> >::BBType *blockPtr;
      blockPtr = new F1();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retPiPi.push_back(ffBase_t<std::complex<double> >::BBHandle_t(blockPtr));
      return retPiPi;
    }

  } // PiPiStar

} // radmat
