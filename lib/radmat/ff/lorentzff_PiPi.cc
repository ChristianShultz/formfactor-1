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
    Tensor<std::complex<double> , 1> F1::operator()(const PInv_t &mom) const
    {
      return convertTensorUnderlyingType<std::complex<double>,double,1>(pPlus(mom) );
    }

    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    // generate a list for the PiPi constructor
    ffBase_t<std::complex<double> >::ff_list genList(void)
    {
      ffBase_t<std::complex<double> >::ff_list retPiPi;
      retPiPi.push_back(ffBase_t<std::complex<double> >::BBHandle_t(new F1()));
      return retPiPi;
    }
    
  }

}
