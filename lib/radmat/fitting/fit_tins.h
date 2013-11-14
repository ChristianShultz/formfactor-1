#ifndef FIT_TINS_H_H_GUARD
#define FIT_TINS_H_H_GUARD

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include "ensem/ensem.h"
#include "semble/semble_vector.h"
#include "radmat/utils/handle.h"
#include "radmat/llsq/llsq_q2_pack.h"
#include "jackFitter/three_point_fit_forms.h"

namespace radmat
{


  // find the flat region, figure out if we want the real or imag part
  struct TinsFitter
  {

    // construct
    TinsFitter(void)
      : didFit(false) 
    { 
      Q2.resize(1);
      Q2 = ENSEM::toDouble(10000000.);
    }

    // template to do the fit with real or complex data -- needs work
    template<typename T>
      void fit(const std::string &filenameBase, 
          const LLSQRet_ff_Q2Pack<T> &data,
          const ThreePointComparatorProps_t &fitProps, 
          const int tsrc,
          const int tsnk);

    // get the form factors at this q2
    std::pair<ENSEM::EnsemReal, SEMBLE::SembleVector<double> > fetchFF(void) const;

    // get a single form factor at this q2
    ENSEM::EnsemReal getFF(const int ffnum) const; 

    // get the fit associated with ffnum
    rHandle<FitThreePoint> getFit(const int ffnum) const;

    // get the ensemble value for q2
    ENSEM::EnsemReal getQ2(void) const {return Q2;}

    // how many form factors are there
    int nFF(void) const {return ff.getN();}

    // write out the fit log files 
    void writeFitLogs(const std::string &path) const;

    // write the fits with the components also plotted
    void writeFitPlotsWithComponents(const std::string &path) const; 

    private:

    // plays nicely with the "fit" function
    void doFit(const std::string &filenameBase, 
        const ENSEM::EnsemVectorReal &data, 
        const int ffnum,
        const ThreePointComparatorProps_t &fitProps,
        const int tsrc, 
        const int tsnk);

    // data store
    bool didFit;
    std::map<int,rHandle<FitThreePoint> > fitters; 
    ENSEM::EnsemReal Q2;
    SEMBLE::SembleVector<double> ff;
  };


  // pre declare the template specializations so we can put them in the cc file
  // won't compile if its called with int or char or something
  template<>
    void TinsFitter::fit<double>(const std::string &, 
        const LLSQRet_ff_Q2Pack<double> &,
        const ThreePointComparatorProps_t &,
        const int tsrc,
        const int tsnk);

  template<>
    void TinsFitter::fit<std::complex<double> >(const std::string &,
        const LLSQRet_ff_Q2Pack<std::complex<double> > &,
        const ThreePointComparatorProps_t &fitProps,
        const int tsrc,
        const int tsnk);





} // namespace radmat




#endif
