// lorentzff_PiPi.cc -
//
// Wednesday, May 30 2012
//

#include"lorentzff_PiPi.h"

namespace radmat
{

  namespace PiPi
  {

    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    std::string F1::ff(void) const
    {
      return std::string(" F_1(Q^2) p_+^{\\mu} ");
    }

    // return a complex version of p_+
    Tensor<std::complex<double> , 1> F1::operator()(const Tensor<double,1> &p_f,
						    const Tensor<double,1> &p_i, const double mom_fac) const
    {
      return convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
    }

    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    // generate a list for the PiPi constructor
    ffBase_t<std::complex<double> >::ff_list genList(void)
    {
      ffBase_t<std::complex<double> >::ff_list retPiPi;
      ffBase_t<std::complex<double> >::BBType *blockPtr;
      blockPtr = new F1();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retPiPi.push_back(ffBase_t<std::complex<double> >::BBHandle_t(blockPtr));
      return retPiPi;
    }
    
  }

}
