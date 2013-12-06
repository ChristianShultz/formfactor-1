/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : fit_tins.cc

 * Purpose :

 * Creation Date : 01-08-2012

 * Last Modified : Fri 06 Dec 2013 12:32:32 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "fit_tins.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "semble/semble_vector.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/plot.h"
#include <string>
#include <math.h>
#include <vector>
#include <iostream>
#include <complex>
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"
#include "adat/handle.h"




namespace radmat
{


  namespace 
  {
    inline double rad(const double deg)
    {
      return 3.14159 * deg/180.;
    }

    // run an avg fit on the real and imag bits
    //      to try to find a phase
    ENSEM::EnsemVectorReal 
      convertToReal(const TinsFitter &fitter,
          const ENSEM::EnsemVectorComplex & in, 
          const int tlow,
          const int thigh)
    {
      ENSEM::EnsemVectorReal real, imag;

      double const_real, const_imag; 

      real = ENSEM::real(in);
      imag = ENSEM::imag(in);

      std::vector<double> t;
      for(int i = 0; i < in.numElem(); ++i)
        t.push_back(double(i)); 

      EnsemData ereal(t,real),eimag(t,imag); 
      ereal.hideDataAboveX(thigh - 0.1); 
      eimag.hideDataBelowX(tlow -0.1); 

      ADAT::Handle<FitFunction> freal(new ThreePointConstant), fimag(new ThreePointConstant);  
      JackFit fit_real(ereal,freal), fit_imag(eimag,fimag); 

      fit_real.runAvgFit(); 
      fit_imag.runAvgFit(); 

      const_real = fit_real.getAvgFitParValue(0);
      const_imag = fit_imag.getAvgFitParValue(0);

      // what if it is a longitudial factor or crossing 
      //    we are only getting about 2% precision, 
      //    assume a FF of O(1)
      if( (const_real < 2e-2) && (const_imag < 2e-2) ) 
      {
        // doesn't matter, both are zero to precision  
        return real;
      }

      if ( isnan(const_real) && isnan(const_imag))
      {
        return real; 
      }


      double fit_phase = std::arg(std::complex<double>(const_real,const_imag)); 

      if(fit_phase < -3.*3.14159/4.)
        fit_phase = - fit_phase; 


      double phase = fit_phase ; 

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
        std::cout << "rl = " << const_real << " im = " << const_imag << std::endl;
        std::cout << "for Q2 = " << SEMBLE::toScalar(ENSEM::mean(fitter.getQ2())) << std::endl;  
        std::cout << "used tlow = " << tlow << " thigh = " << thigh << std::endl; 
        SPLASH("check bad_corr.jack, bad_corr.ax for the correlator"); 

        AxisPlot plot; 
        plot.addEnsemData(ENSEM::real(in),"\\sq",1);
        plot.addEnsemData(ENSEM::imag(in),"\\sq",2);
        plot.sendToFile("bad_corr.ax");

        ENSEM::write("bad_corr.jack",in); 

        SPLASH("An error occured while trying to pull the real/imag part of the solution vector,exiting.");
        exit(1);
      }
    }


  } // anonymous

  template<>
    void TinsFitter::single_fit<double>(const std::string &fname, 
        const LLSQRet_ff_Q2Pack<double> &pack,
        const int ff_max, 
        const ThreePointComparatorProps_t &fitProps,
        const int tsrc, 
        const int tsnk)
    {
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,ff_max);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<double>::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,it->second,it->first,fitProps,tsrc,tsnk);

      didFit = true;
    }



  template<>
    void TinsFitter::single_fit<std::complex<double> >(const std::string &fname , 
        const LLSQRet_ff_Q2Pack<std::complex<double> > &pack,
        const int ff_max,
        const ThreePointComparatorProps_t &fitProps,
        const int tsrc,
        const int tsnk)
    {
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,ff_max);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<std::complex<double> >::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,convertToReal(*this,it->second,tsrc,tsnk),it->first,fitProps,tsrc,tsnk);

      didFit = true;
    }


  template<>
    void TinsFitter::fit<double>(const std::string &fname, 
        const LLSQRet_ff_Q2Pack<double> &pack,
        const ThreePointComparatorProps_t &fitProps,
        const int tsrc, 
        const int tsnk)
    {
      const unsigned int sz = pack.size();
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,sz);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<double>::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,it->second,it->first,fitProps,tsrc,tsnk);

      didFit = true;
    }



  template<>
    void TinsFitter::fit<std::complex<double> >(const std::string &fname , 
        const LLSQRet_ff_Q2Pack<std::complex<double> > &pack,
        const ThreePointComparatorProps_t &fitProps,
        const int tsrc,
        const int tsnk)
    {
      const unsigned int sz = pack.size();
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,sz);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<std::complex<double> >::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,convertToReal(*this,it->second,tsrc,tsnk),it->first,fitProps,tsrc,tsnk);

      didFit = true;
    }


  void TinsFitter::doFit(const std::string &filenameBase, 
      const ENSEM::EnsemVectorReal &data, 
      const int ffnum,
      const ThreePointComparatorProps_t &fitProps,
      const int tsrc,
      const int tsnk)
  {
    std::stringstream ss,jack,ax,fit;
    ss << filenameBase;
    fit << ss.str() << "FF_" << ffnum << "_fit.jack";
    jack << ss.str() << "FF_" << ffnum << ".jack";
    ax << ss.str() << "FF_" << ffnum << ".ax";

    std::vector<double> time; 
    const int Lt = data.numElem();

    for(int t = 0; t < Lt; t++)
      time.push_back(t);

    EnsemData corrData(time,data);


    ADAT::Handle<FitComparator> fitComp = constructThreePointFitComparator(fitProps);
    // NB: I Have assumed that no chopping has gone on in the data
    rHandle<FitThreePoint> fitCorr (new FitThreePoint(corrData,tsnk,tsrc,
          fitProps.thigh,fitProps.tlow,fitComp,fitProps.minTSlice,fitProps.fit_type));

    // send to files
    fitCorr->saveFitPlot(ax.str());
    write(jack.str(),data);

    ff.loadEnsemElement(ffnum,fitCorr->getFF());
    write(fit.str(),fitCorr->getFF()); 

    fitters.insert(std::map<int,rHandle<FitThreePoint> >::value_type(ffnum,fitCorr));
  }


  void TinsFitter::writeFitLogs(const std::string &path) const
  {
    std::map<int,rHandle<FitThreePoint> >::const_iterator it;
    for(it = fitters.begin(); it != fitters.end(); ++it)
    {
      std::stringstream ss;
      ss << path << "F_" << it->first << ".fitlog";
      std::ofstream out(ss.str().c_str());
      out << it->second->getFitSummary();
      out.close();
    }

  }

  void TinsFitter::writeFitPlotsWithComponents(const std::string & path) const
  {
    std::map<int,rHandle<FitThreePoint> >::const_iterator it;
    for(it = fitters.begin(); it != fitters.end(); ++it)
    {
      std::stringstream ss; 
      ss << path <<  "FF_" << it->first << ".ax";
      std::ofstream out(ss.str().c_str());
      out << it->second->getFitPlotStringWithComponents(); 
      out.close();
    }

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

  rHandle<FitThreePoint> TinsFitter::getFit(const int ffnum) const
  {
    POW2_ASSERT(fitters.find(ffnum) != fitters.end());
    return rHandle<FitThreePoint>(fitters.find(ffnum)->second);
  }


} // namespace radmat
