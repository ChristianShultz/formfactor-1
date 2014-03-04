/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_formfactor_data.cc

 * Purpose :

 * Creation Date : 21-02-2014

 * Last Modified : Tue 04 Mar 2014 02:51:57 PM EST

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
#include "radmat/utils/printer.h"
#include "adat/handle.h"



namespace radmat
{


  namespace 
  {
    // swapping in empty functions should make the compiler 
    // optimize these print statements away.. plus its fun
    struct dimension_printer
    {
      static void print(const std::string &msg)
      {}
      //      { std::cout << msg << std::endl; }
    };

    struct inp_dimension_printer
    {
      static void print(const std::string &msg)
      {}
      //      { std::cout << msg << std::endl; }
    };

    struct ret_dimension_printer
    {
      static void print(const std::string &msg)
      {}
      //      { std::cout << msg << std::endl; }
    };

    struct mean_printer
    {
      static void print(const std::string &s)
      {}
      // { std::cout << s << std::endl;}
    };

    struct case_printer
    {
      static void print(const std::string &s)
      {}
     // { std::cout << s << std::endl;}
    };

    struct ensem_printer
    {
      static void print(const std::string &s)
      {}
      // {std::cout << s << std::endl;}
    };

    template<typename T>
      std::string to_string( T t )
      {
        std::stringstream ss; 
        ss << t ;
        return ss.str(); 
      }



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

        // check code
        std::stringstream ssr,ssi;
        ssr << "tl " << tlow << " th " << thigh;
        ssi << "tl " << tlow << " th " << thigh;

        for(int i = 0; i < in.numElem(); ++i)
        {
          ssr << " " << SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(real,i)));
          ssi << " " << SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(imag,i)));
        }

        printer_function<ensem_printer>("real - " + ssr.str());
        printer_function<ensem_printer>("imag - " + ssi.str());
        printer_function<ensem_printer>("in.numElem() " + to_string(in.numElem()));
        printer_function<ensem_printer>("in.size() " + to_string(in.size()));

        // time variable, make ensem data 
        std::vector<double> t;
        for(int i = 0; i < in.numElem(); ++i)
          t.push_back(double(i)); 

        EnsemData ereal(t,real),eimag(t,imag); 
        ereal.hideDataAboveX(thigh + 0.1); 
        eimag.hideDataBelowX(tlow -0.1); 

 
//        // the fitting is acting up, check that we get out what 
//        // we put into the ensem data -- passes
//        ENSEM::EnsemVectorReal insanityi, insanityr; 
//        insanityr = ereal.getAllYData(); 
//        insanityi = eimag.getAllYData(); 
//
//        ssr.str("");
//        ssr.clear();
//        ssi.str("");
//        ssi.clear(); 
//
//        for(int i = 0; i < in.numElem(); ++i)
//        {
//          ssr << " " << SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(real,i)));
//          ssi << " " << SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(imag,i)));
//        }
//
//        printer_function<ensem_printer>("insanity real - " + ssr.str());
//        printer_function<ensem_printer>("insanity imag - " + ssi.str());


        // do fits
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


//        // the fitting is acting up, check that we get out what 
//        // we put into the ensem data -- also passes 
//        insanityr = fit_real.data.getAllYData(); 
//        insanityi = fit_real.data.getAllYData(); 
//
//        ssr.str("");
//        ssr.clear();
//        ssi.str("");
//        ssi.clear(); 
//
//        for(int i = 0; i < in.numElem(); ++i)
//        {
//          ssr << " " << SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(real,i)));
//          ssi << " " << SEMBLE::toScalar( ENSEM::mean( ENSEM::peekObs(imag,i)));
//        }
//
//        printer_function<ensem_printer>("insanity2 real - " + ssr.str());
//        printer_function<ensem_printer>("insanity2 imag - " + ssi.str());
//


        // is the constant consistent with zero?
        if( fabs(const_real) - 3.*sqrt(fabs(const_real_var)) < 0.)
        {
          printer_function<case_printer>("real is consistent with zero");
          printer_function<case_printer>("real = " + to_string(const_real) 
              + "+/-" + to_string(sqrt(fabs(const_real_var)))); 
          printer_function<case_printer>("imag = " + to_string(const_imag) 
              + "+/-" + to_string(sqrt(fabs(const_imag_var)))); 
          if (const_imag > 0. )
            return phase_pair(IP,imag); 
          return phase_pair(IM,imag);
        }

        if( fabs(const_imag) - 3.*sqrt(fabs(const_imag_var)) < 0.)
        {
          printer_function<case_printer>("imag is consistent with zero");
          printer_function<case_printer>("imag = " + to_string(const_imag) 
              + "+/-" + to_string(sqrt(fabs(const_imag_var)))); 
          if( const_real > 0. )
            return phase_pair(RP,real); 
          return phase_pair(RM,real);
        }

        // what if it is a longitudial factor or crossing 
        //    we are only getting about 2% precision, 
        //    assume a FF of O(1)
        if( (const_real < 2e-2) && (const_imag < 2e-2) ) 
        {
          printer_function<case_printer>("ff is consistent with zero");
          // doesn't matter, both are zero to precision  
          return phase_pair(ZERO,real);
        }



        // something unexpected happened
        if ( isnan(const_real) || isnan(const_imag))
        {
          // the fitter seems to fail on ensembles with zero variance??
          // 
          // with svd resetting it is very easy to get an ensemble 
          // with mean zero and variance zero thus we have to 
          // work harder which makes christian cranky since its sunday 

          // guard zero variance here.. 
          // this is completely nuts, someone needs to update the stupid fitter

          // need to use assignment operator, no builtin ensem constructors
          SEMBLE::SembleVector<double> foor; foor = real; 
          SEMBLE::SembleVector<double> fooi; fooi = imag; 

          // pull down the copies, grab the variance using semble
          // then compare it to zero, if the test passes check for 
          // either the real part or the imag part being explicitly 
          // zero, then return the other with the correct phase
          itpp::Vec<double> bar( foor.getN() ), bazr,bazi; 
          bar.zeros(); 
          bazr = foor.variance(); 
          bazi = fooi.variance(); 

          if( bazr == bar ) 
            if( bazi == bar )
            {
              printer_function<case_printer>("zero variance encountered"); 
              bazr = foor.mean(); 
              bazi = fooi.mean(); 
              if ( ( bazr == bar ) && ( bazi == bar ) )
                return phase_pair(ZERO,real); 
              if ( bazr == bar )
                return bazi(0) > 0 ? phase_pair(IP,imag) : phase_pair(IM,imag);
              if( bazi == bar )
                return bazr(0) > 0 ? phase_pair(RP,real) : phase_pair(RM,real); 
            }

          // otherwise die one function up since this is stupid 

          SPLASH("encountered nan: check bad_corr.nan.jack, bad_corr.nan.ax for the correlator"); 

          AxisPlot plot; 
          plot.addEnsemData(ENSEM::real(in),"\\sq",1);
          plot.addEnsemData(ENSEM::imag(in),"\\sq",2);
          plot.sendToFile("bad_corr.nan.ax");

          ENSEM::write("bad_corr.nan.jack",in); 

          return phase_pair(ERROR,real); 
        }


        double fit_phase = std::arg(std::complex<double>(const_real,const_imag)); 

        if(fit_phase < -3.*3.14159/4.)
          fit_phase = - fit_phase; 


        double phase = fit_phase ; 

        printer_function<case_printer>("general cases encountered"); 

        if( (phase < 0.174528) && (phase > -0.174708) ) // +/- 10 degree about 0 in rad
        {
          printer_function<case_printer>("RP encountered"); 
          return phase_pair(RP,real);
        }
        else if( (phase > 1.39622) && (phase < 1.74528)) // 90deg
        {
          printer_function<case_printer>("IP encountered"); 
          return phase_pair(IP,imag);
        }
        else if( (phase > 2.96697) || (phase < -2.96715))  //180deg // this is atan2 specific, it returns (-pi,pi)
        {
          printer_function<case_printer>("RM encountered"); 
          return phase_pair(RM,real);
        }
        else if( (phase > -1.74546) && (phase < -1.3964)) // 270 deg
        {
          printer_function<case_printer>("IM encountered"); 
          return phase_pair(IM,imag);
        }
        else
        {
          printer_function<case_printer>("unknown phase encountered"); 
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
        int bns = d.getB(); 
        int sz = d.getN(); 
        e.resize(bns); 
        e.resizeObs(sz); 
        for(int i = 0; i < sz; ++i)
          ENSEM::pokeObs( e, d.getEnsemElement(i) , i);

        // run the fits on the central 70% of the correlator, this 
        // is a completely arbitrary choice 
        return convert_to_real( e, int(sz*0.15) , int(sz*0.85) ); 
      }

    // fine the phase for all of the vectors
    std::map<std::string, phase_pair> 
      find_phases( const LLSQComplexFormFactorData_t &d )
      {
        std::map<std::string,phase_pair> ret; 
        LLSQComplexFormFactorData_t::const_iterator it; 
        for(it = d.begin(); it != d.end(); ++it)
        {
          printer_function<inp_dimension_printer>(
              "input " +  it->first 
              + " N " + to_string( it->second.getN() )
              + " B " + to_string( it->second.getB() ) ); 

          printer_function<mean_printer>(
              "input_corr " +  it->first 
              + to_string( it->second.mean() ) ); 
          ret.insert( std::make_pair(it->first, find_phase(it->second) ) ); 
        }

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
          std::cout << it->first << " was bad" << std::endl;
          exit(1); 
        }

      // above guards ERROR
      PHASE expectedA,expectedB;  
      for(it = mappy.begin(); it != mappy.end(); ++it)
        if( it->second.first != ZERO )
          expectedA = it->second.first; 

      // so it was something
      if( expectedA == RP )
        expectedB = RM; 
      if( expectedA == RM )
        expectedB = RP; 
      if( expectedA == IP )
        expectedB = IM; 
      if( expectedA == IM )
        expectedB = IP; 

      // check that they are all either real , imag , or zero
      for(it = mappy.begin(); it != mappy.end(); ++it)
      {
        printer_function<dimension_printer>( "check_phases" + it->first 
            + " numElem = " + to_string(it->second.second.numElem()) 
            + " size = " + to_string(it->second.second.size()) ); 
        if( (it->second.first != expectedA )
            &&(it->second.first != expectedB) 
            &&(it->second.first != ZERO) )
        {
          std::cout << __PRETTY_FUNCTION__ 
            << "\nerror: encountered, unexpected phases, exiting" << std::endl;
          std::cout << "expected " << expectedA << " " << expectedB 
            << " or " << ZERO << " got " << it->second.first << std::endl;
          std::cout << "keys \nRP->" << RP 
            << "\nRM->" << RM
            << "\nIP->" << IP
            << "\nIM->" << IM
            << "\nZZ->" << ZERO 
            << std::endl;
          exit(1); 
        }
      }
    }

    // run checks then push the result into the return data
    LLSQRealFormFactorData_t
      do_work_local(const LLSQComplexFormFactorData_t &d)
      {
        printer_function<inp_dimension_printer>("entering rephase"); 

        std::map<std::string, phase_pair> mappy = find_phases(d); 
        check_phases(mappy); 

        printer_function<ret_dimension_printer>(
            " map size " + to_string(mappy.size()) ); 

        LLSQRealFormFactorData_t ret; 

        // update momentum transfer        
        ret.Qsq = d.Q2(); 

        // load elements for fit
        std::map<std::string,phase_pair>::const_iterator it; 
        for(it = mappy.begin(); it != mappy.end(); ++it)
        {
          printer_function<ret_dimension_printer>( 
              "output " + it->first 
              + " numElem " + to_string(it->second.second.numElem()) 
              + " size " + to_string(it->second.second.size()));  
          ret.mappy[it->first] = it->second.second; 
        }

        // possibly print output result
        LLSQRealFormFactorData_t::const_iterator pit; 
        for(pit = ret.begin(); pit != ret.end(); ++pit)
          printer_function<mean_printer>(
              "output_corr " +  pit->first 
              + to_string( pit->second.mean() ) ); 

        printer_function<ret_dimension_printer>("exiting rephase");

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
