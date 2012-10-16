#ifndef LORENTZFF_PIPI_H_H_GUARD
#define LORENTZFF_PIPI_H_H_GUARD

#include "ff_gen_llsq_row.h"
#include <complex>

namespace radmat
{

  namespace PiPi
  {

    /*
DECOMPOSITION:

\begin{align*}
p_+ &\equiv p_f + p_i  \qquad p_- \equiv p_f - p_i \\
\langle M^{\mu} \rangle &= F_+ p_+^{\mu} + F_- p_-^{\mu} \\
q_{\mu}\langle M^{\mu} \rangle = 0 &\rightarrow F_+p_+^{\mu}p_{-\mu} - F_-Q^2 = 0 \\
&F_+\left(m_{\pi}^2 -m_{\pi}^2 \right) -F_-Q^2 = 0 \rightarrow F_- = 0 \\
\\
\langle M^{\mu} \rangle &= F_1(Q^2) p_+^{\mu}
\end{align*}
     */


    // only one ff
    struct F1 : public ffBlockBase_t<std::complex<double> >
    {
      std::string ff(void) const;
      Tensor<std::complex<double> , 1> operator()(const Tensor<double,1> &p_f, 
          const Tensor<double,1> &p_i) const;
    };


    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////


    // generate a list for the PiPi constructor
    ffBase_t<std::complex<double> >::ff_list genList(void);


    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    // only need to derive the constructor, everything else is in the 
    // base class, do this polymorphically, make some function that 
    // returns the appropriate handle based on the requested matrix element type
    struct PiPi : public ffBase_t<std::complex<double> >
    {
      PiPi(void)
        : ffBase_t<std::complex<double> >(radmat::PiPi::genList())  
      {  } 

      PiPi& operator=(const PiPi &o)
      {

        if(this != &o)
          ffBase_t<std::complex<double> >::operator=(o);
        return *this;
      }

      // no slicing
      PiPi(const PiPi &o)
        : ffBase_t<std::complex<double> >(o)
      {  }

      private:
      // I'm not sure if these could inherit so we will hide them as well
      PiPi(const ffBase_t<std::complex<double> >::ff_list &);
      PiPi(const ffBase_t<std::complex<double> >::ff_list);
    };


  } // close PiPi

} // close radmat

#endif
