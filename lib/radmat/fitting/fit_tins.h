#ifndef FIT_TINS_H_H_GUARD
#define FIT_TINS_H_H_GUARD

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include "ensem/ensem.h"
#include "semble/semble_vector.h"
#include "fit_three_point_ff.h"
#include "radmat_fit_forms.h"
#include "adat/handle.h"
#include "radmat/llsq/llsq_q2_pack.h"

namespace radmat
{
  // find the flat region, figure out if we want the real or imag part
  struct TinsFitter
  {
    TinsFitter(void)
      : didFit(false) 
    { 
      Q2.resize(1);
      Q2 = ENSEM::toDouble(10000000.);
     }

    template<typename T>
      void fit(const std::string &filenameBase, const LLSQRet_ff_Q2Pack<T> &data);

    std::pair<ENSEM::EnsemReal, SEMBLE::SembleVector<double> > fetchFF(void) const;
    ENSEM::EnsemReal getFF(const int ffnum) const; 
    ADAT::Handle<Fit3PtFF_Base> getFit(const int ffnum) const;
    ENSEM::EnsemReal getQ2(void) const {return Q2;}
    int nFF(void) const {return ff.getN();}

    private:

    void doFit(const std::string &filenameBase, const ENSEM::EnsemVectorReal &data, const int ffnum);

    bool didFit;
    std::map<int,ADAT::Handle<Fit3PtFF_Base> > fitters; 
    ENSEM::EnsemReal Q2;
    SEMBLE::SembleVector<double> ff;
  };


  template<>
    void TinsFitter::fit<double>(const std::string &, const LLSQRet_ff_Q2Pack<double> &);

  template<>
    void TinsFitter::fit<std::complex<double> >(const std::string &,
        const LLSQRet_ff_Q2Pack<std::complex<double> > &);





} // namespace radmat




#endif
