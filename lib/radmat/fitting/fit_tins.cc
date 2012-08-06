/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : fit_tins.cc

 * Purpose :

 * Creation Date : 01-08-2012

 * Last Modified : Mon Aug  6 10:58:51 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "fit_tins.h"
#include "fit_three_point_ff.h"
#include "radmat_fit_forms.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "semble/semble_vector.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "adat/handle.h"
#include <string>
#include <vector>
#include <iostream>
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"




namespace radmat
{


  namespace 
  {
    inline double rad(const double deg)
    {
      return 3.14159 * deg/180.;
    }

    ENSEM::EnsemVectorReal convertToReal(const ENSEM::EnsemVectorComplex & in)
    {
      ENSEM::EnsemVectorReal real, imag;
      ENSEM::EnsemReal ephase;

      real = ENSEM::real(in);
      imag = ENSEM::imag(in);

      ephase = ENSEM::atan2(ENSEM::peekObs(imag,0),ENSEM::peekObs(real,0));

      double phase = SEMBLE::toScalar( ENSEM::mean(ephase) );

      if( (phase < 0.174528) && (phase > -0.174708) ) // +/- 10 degree about 0 in rad
        return real;
      else if( (phase > 1.39622) && (phase < 1.74528)) // 90deg
        return imag;
      else if( (phase > 2.96697) || (phase < -2.96715))  //180deg // this is atan2 specific, it returns (-pi,pi)
        return real;
      else if( (phase > -1.74546) && (phase < -1.3964)) // 270 de/doFit
        return imag;
      else
      {
        std::cout << "The calculated phase was " << phase*180./3.14159 << " (deg)" << std::endl; 
        SPLASH("An error occured while trying to pull the real/imag part of the solution vector,exiting.");
        exit(1);
      }
    }

  } // anonymous


  template<>
    void TinsFitter::fit<double>(const std::string &fname, const LLSQRet_ff_Q2Pack<double> &pack)
    {
      const unsigned int sz = pack.size();
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,sz);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<double>::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,it->second,it->first);

      didFit = true;
    }



  template<>
    void TinsFitter::fit<std::complex<double> >(const std::string &fname , 
        const LLSQRet_ff_Q2Pack<std::complex<double> > &pack)
    {
      const unsigned int sz = pack.size();
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,sz);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<std::complex<double> >::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,convertToReal(it->second),it->first);

      didFit = true;
    }


  void TinsFitter::doFit(const std::string &filenameBase, 
      const ENSEM::EnsemVectorReal &data, 
      const int ffnum)
  {
    std::stringstream ss,jack,ax;
    ss << filenameBase;
    jack << ss.str() << ".jack";
    ax << ss.str() << ".ax";

    std::vector<double> time; 
    const int Lt = data.numElem();

    for(int t = 0; t < Lt; t++)
      time.push_back(t);

    EnsemData corrData(time,data);
    ADAT::Handle<FitComparator> fitComp(new CompareFitsByChisqPerNDoF);
    ADAT::Handle<Fit3PtFF_Base> fitCorr(new Fit3PtFF_Constant(corrData,fitComp,0.1,5));

    // send to files
    fitCorr->saveFitPlot(ax.str());
    write(jack.str(),data);

    ff.loadEnsemElement(ffnum,fitCorr->getFF());

    fitters.insert(std::map<int,ADAT::Handle<Fit3PtFF_Base> >::value_type(ffnum,fitCorr));
  }


  std::pair<ENSEM::EnsemReal,SEMBLE::SembleVector<double> > TinsFitter::fetchFF(void) const
  {
    POW2_ASSERT(didFit);

    return std::pair<ENSEM::EnsemReal,SEMBLE::SembleVector<double> >(Q2,ff);
  }


  ENSEM::EnsemReal TinsFitter::getFF(const int index) const
  {
    POW2_ASSERT(index < ff.getN());
    return ff.getEnsemElement(index);
  }

  ADAT::Handle<Fit3PtFF_Base> TinsFitter::getFit(const int ffnum) const
  {
    POW2_ASSERT(fitters.find(ffnum) != fitters.end());
    return ADAT::Handle<Fit3PtFF_Base>(fitters.find(ffnum)->second);
  }


} // namespace radmat
