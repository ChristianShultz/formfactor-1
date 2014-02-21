/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_formfactor_data.cc

 * Purpose :

 * Creation Date : 21-02-2014

 * Last Modified : Fri 21 Feb 2014 04:26:30 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "llsq_formfactor_data.h"
#include "ensem/ensem.h"
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

    enum PHASE
    {
      RP,
      RM,
      IP,
      IM,
      ZERO,
      ERROR
    };

    inline double rad(const double deg)
    {
      return 3.14159 * deg/180.;
    }

    // carpal tunnel prevention 
    typedef std::pair<PHASE,ENSEM::EnsemVectorReal> phase_pair;

    // run an avg fit on the real and imag bits
    //      to try to find a phase
    phase_pair
      convert_to_real(const ENSEM::EnsemVectorComplex & in, 
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

        // no on has time for this 
        //   fit_real.runJackFit(); 
        //   fit_imag.runJackFit(); 
        //   const_real = fit_real.getAvgFitParValue(0);
        //   const_imag = fit_imag.getAvgFitParValue(0);

        double const_real_var = fit_real.getAvgFitParError(0);
        double const_imag_var = fit_imag.getAvgFitParError(0);

        // is the constant consistent with zero?
        if( fabs(const_real) - 3.*fabs(const_real_var) < 0.)
        {
          if (const_imag > 0. )
            return phase_pair(IP,imag); 
          return phase_pair(IM,imag);
        }

        if( fabs(const_imag) - 3.*fabs(const_imag_var) < 0.)
        {
          if( const_real > 0. )
            return phase_pair(RP,real); 
          return phase_pair(RM,real);
        }

        // what if it is a longitudial factor or crossing 
        //    we are only getting about 2% precision, 
        //    assume a FF of O(1)
        if( (const_real < 2e-2) && (const_imag < 2e-2) ) 
        {
          // doesn't matter, both are zero to precision  
          return phase_pair(ZERO,real);
        }

        // something unexpected happened
        if ( isnan(const_real) || isnan(const_imag))
        {
          return phase_pair(ERROR,real); 
        }


        double fit_phase = std::arg(std::complex<double>(const_real,const_imag)); 

        if(fit_phase < -3.*3.14159/4.)
          fit_phase = - fit_phase; 


        double phase = fit_phase ; 

        if( (phase < 0.174528) && (phase > -0.174708) ) // +/- 10 degree about 0 in rad
          return phase_pair(RP,real);
        else if( (phase > 1.39622) && (phase < 1.74528)) // 90deg
          return phase_pair(IP,imag);
        else if( (phase > 2.96697) || (phase < -2.96715))  //180deg // this is atan2 specific, it returns (-pi,pi)
          return phase_pair(RM,real);
        else if( (phase > -1.74546) && (phase < -1.3964)) // 270 de/doFit
          return phase_pair(IM,imag);
        else
        {
          std::cout << "The calculated phase was " << phase*180./3.14159 << " (deg)" << std::endl;
          std::cout << "rl = " << const_real << " im = " << const_imag << std::endl;
          std::cout << "used tlow = " << tlow << " thigh = " << thigh << std::endl; 
          SPLASH("check bad_corr.jack, bad_corr.ax for the correlator, returning zero"); 

          AxisPlot plot; 
          plot.addEnsemData(ENSEM::real(in),"\\sq",1);
          plot.addEnsemData(ENSEM::imag(in),"\\sq",2);
          plot.sendToFile("bad_corr.ax");

          ENSEM::write("bad_corr.jack",in); 

        }

        return phase_pair( ERROR, SEMBLE::toScalar( double(0.) ) * real ); 
      }


    // find the phase for a single vector
    phase_pair
      find_phase( const SEMBLE::SembleVector<std::complex<double> > &d)
      {
        ENSEM::EnsemVectorComplex e; 
        int sz = d.getN(); 
        e.resizeObs(sz); 
        for(int i = 0; i < sz; ++i)
          ENSEM::pokeObs( e, d.getEnsemElement(i) , i);

        return convert_to_real( e, int(sz*0.15) , int(sz*0.85) ); 
      }

    // fine the phase for all of the vectors
    std::map<std::string, phase_pair> 
      find_phases( const LLSQComplexFormFactorData_t &d )
      {
        std::map<std::string,phase_pair> ret; 
        LLSQComplexFormFactorData_t::const_iterator it; 
        for(it = d.begin(); it != d.end(); ++it)
          ret.insert( std::pair<std::string,phase_pair>(it->first, find_phase(it->second) ) ); 

        return ret; 
      }

    // check that all phases are consistent
    void check_phases( const std::map<std::string, phase_pair> &mappy)
    {
      std::map<std::string,phase_pair>::const_iterator it; 
      for(it = mappy.begin(); it != mappy.end(); ++it)
        if( it->second.first == ERROR )
        {
          std::cout << __PRETTY_FUNCTION__ << ": error encountered, exiting" << std::endl;
          exit(1); 
        }

      // above guards ERROR
      PHASE expectedA,expectedB;  
      for(it = mappy.begin(); it != mappy.end(); ++it)
        if( it->second.first != ZERO )
          expectedA = it->second.first; 

      // so it was something
      if( expectedA == RP )
        expectedB == RM; 
      if( expectedA == RM )
        expectedB == RP; 
      if( expectedA == IP )
        expectedB == IM; 
      if( expectedA == IM )
        expectedB == IP; 

      // check that they are all either real , imag , or zero
      for(it = mappy.begin(); it != mappy.end(); ++it)
        if( (it->second.first != expectedA )
            &&(it->second.first != expectedB) 
            &&(it->second.first != ZERO) )
        {
          std::cout << __PRETTY_FUNCTION__ 
            << ": error encountered, unexpected phases, exiting" << std::endl;
          exit(1); 
        }
    }

    // run checks then push the result into the return data
    LLSQRealFormFactorData_t
      do_work_local(const LLSQComplexFormFactorData_t &d)
      {
        std::map<std::string, phase_pair> mappy = find_phases(d); 
        check_phases(mappy); 

        LLSQRealFormFactorData_t ret; 
        std::map<std::string,phase_pair>::const_iterator it; 
        for(it = mappy.begin(); it != mappy.end(); ++it)
          ret.mappy[it->first] = it->second.first; 

        return ret; 
      }


  } // anonymous



  // callback 
  LLSQRealFormFactorData_t 
    rephase_formfactor_data( const LLSQComplexFormFactorData_t &d)
    {
      return do_work_local(d); 
    }

} // radmat
