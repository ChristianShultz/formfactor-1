/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : axis_plotter.cc

 * Purpose :

 * Creation Date : 31-07-2012

 * Last Modified : Mon Nov 26 14:34:44 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "axis_plotter.h"
#include "ensem/ensem.h"
#include "jackFitter/plot.h"
#include "semble/semble_file_management.h"
#include "radmat/load_data/build_q2_packs.h"
#include "radmat/utils/splash.h"
#include <sstream>



namespace radmat
{




  template<>
    void MakeAxisPlots::plot(const ThreePointCorrelator<double> &C3pt,
        const std::string &path, 
        const std::string &fname) const

    {

      SEMBLE::SEMBLEIO::makeDirectoryPath(path);

      ENSEM::EnsemVectorReal cor = C3pt.C3pt;
      std::stringstream axis, jack;

      axis << fname << ".ax";
      jack << fname << ".jack";


      ENSEM::write(jack.str(),cor);

      AxisPlot plot;
      plot.addEnsemData(cor,"\\sq",1);

      plot.sendToFile(axis.str());

    }


  template<>
    void MakeAxisPlots::plot(const ThreePointCorrelator<std::complex<double> > &C3pt,
        const std::string &path,
        const std::string &fname) const
    {
      SEMBLE::SEMBLEIO::makeDirectoryPath(path);

      ENSEM::EnsemVectorComplex tmp;
      tmp = C3pt.C3pt;


      std::stringstream ss_rawj , ss_rawx, ss_real, ss_cplx;

      ss_rawj << fname << "_raw.jack";
      ss_rawx << fname << "_raw.ax" ;
      ss_real << fname << "_real.ax";
      ss_cplx << fname << "_cplx.ax";

      write(ss_rawj.str(),tmp);

      AxisPlot plot,Rplot,Cplot;
      ENSEM::EnsemVectorReal Creal,CComplex;

      Creal = ENSEM::real(tmp);
      CComplex = ENSEM::imag(tmp);

      Rplot.addEnsemData(Creal,"\\sq",1);
      Cplot.addEnsemData(CComplex,"\\sq",2);
      plot.addEnsemData(Creal,"\\sq",1);
      plot.addEnsemData(CComplex,"\\sq",2);

      Rplot.sendToFile(std::string(ss_real.str()));
      Cplot.sendToFile(std::string(ss_cplx.str()));
      plot.sendToFile(std::string(ss_rawx.str()));


    }

  namespace
  {
    template<typename T> 
      typename SEMBLE::PromoteEnsem<T>::Type getNorm(const Fake3ptKey &k, const int t);

    template<>
      ENSEM::EnsemComplex getNorm<std::complex<double> >(const Fake3ptKey &k, const int t)
      {
        ThreePtPropagationFactor<std::complex<double> > p;
        return p(k.E_sink,k.Z_sink,k.t_sink,t,k.E_source,k.Z_source,k.t_source);  
      }

    template<>
      ENSEM::EnsemReal getNorm<double>(const Fake3ptKey &k, const int t)
      {
        ThreePtPropagationFactor<double> p;
        return p(k.E_sink,ENSEM::real(k.Z_sink),k.t_sink,t,k.E_source,ENSEM::real(k.Z_source),k.t_source);  
      }
  } // namespace anonomyous 


  template<> 
    void MakeAxisPlots::plot(const Fake3ptCorr<double> &C3pt)const
    {
      std::string  path =  SEMBLE::SEMBLEIO::getPath();

      // convert Q2 to a string format
      double Q2 = SEMBLE::toScalar(ENSEM::mean(C3pt.getKey().Q2));
      std::stringstream ssQ2;
      ssQ2 << Q2;
      std::string sQ2 = ssQ2.str();
      std::replace(sQ2.begin(),sQ2.end(),'.','p');
      std::replace(sQ2.begin(),sQ2.end(),'-','m');

      // set up the output paths
      path += std::string("FakeDataPlots/");
      SEMBLE::SEMBLEIO::makeDirectoryPath(path);

      std::string corr_path = path + std::string("Corrs/");
      std::string plot_path = path + std::string("Plots/"); 

      SEMBLE::SEMBLEIO::makeDirectoryPath(corr_path);
      SEMBLE::SEMBLEIO::makeDirectoryPath(plot_path); 

      // set up the base filename

      std::stringstream fname_base;
      fname_base << "C3pt_l" << C3pt.getKey().lorentz  << "_hf" << C3pt.getKey().hel_sink 
        << "_hi" << C3pt.getKey().hel_source << "_ID_" << C3pt.getKey().elemIDBase << "_pf" << C3pt.getKey().mom.momSink[0] << "_"
        << C3pt.getKey().mom.momSink[1] << "_" << C3pt.getKey().mom.momSink[2] << "_pi" << C3pt.getKey().mom.momSource[0] << "_"
        << C3pt.getKey().mom.momSource[1] << "_" << C3pt.getKey().mom.momSource[2] << "_Q2_" << sQ2;

      // get the loop variable bounds

      int nsource = C3pt.getNSource(); 
      int nsink = C3pt.getNSink();
      int Lt = C3pt.getLt();
      int ncfg = C3pt.getNCfg(); 

      ENSEM::EnsemVectorReal norm, normed_corr; 
      norm.resize(ncfg);
      norm.resizeObs(Lt);
      norm = SEMBLE::toScalar(double(0));
      normed_corr = C3pt.get3pt();

      std::vector<double> time, mean_norm_corr,mme,mpe; 

      for(int t = 0; t < Lt; ++t)
      {
        ENSEM::pokeObs(norm,getNorm<double>(C3pt.getKey(),t), t);
        time.push_back(double(t));
        ENSEM::pokeObs(normed_corr, ENSEM::peekObs(normed_corr,t)/ENSEM::peekObs(norm,t),t); 
        double mean = SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(normed_corr,t)));
        double err =  sqrt( SEMBLE::toScalar(ENSEM::variance(ENSEM::peekObs(normed_corr,t))));
        mme.push_back(mean - err);
        mpe.push_back(mean + err); 
      }

      AxisPlot plot; 

      double high = * std::max_element(mpe.begin(),mpe.end());
      double low = * std::min_element(mme.begin(), mme.end()); 

      // do work!
      for(int source = 0; source < nsource; ++source)
        for(int sink = 0; sink < nsink; ++sink)
        {
          // record the ith_jth component of the sum that contributed to the correlator
          std::stringstream jack;
          jack << corr_path << fname_base << "_source" << source << "_sink" << sink << ".jack";
          ENSEM::EnsemVectorReal work = C3pt.get3ptComponent(source,sink); 
          ENSEM::write(jack.str(),work); 

          std::vector<double> line;

          for(int t = 0; t < Lt; ++t)
          {
            ENSEM::pokeObs(work,ENSEM::peekObs(work,t)/ENSEM::peekObs(norm,t),t); 
            line.push_back(SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(work,t))));
          }

          plot.addLineData(time,line, (sink*nsource + source) % 5 + 1); 


          double tmp = *std::min(line.begin(),line.end());
          low = std::min(low,tmp);
          tmp = *std::max(line.begin(),line.end());
          high = std::max(high,tmp); 
        }


      double height = high - low; 
      double unit = height / double(nsource * nsink);


      for(int source = 0; source < nsource; ++source)
        for(int sink = 0; sink < nsink; ++sink)
        {        
          std::stringstream ss; 
          ss << sink << "_" << source;
          plot.addLabel( Lt + 1, unit * (sink*nsource + source) + low, ss.str(), (sink*nsource + source) % 5 + 1, 0.8);
        }

      plot.addEnsemData(normed_corr, "\\sq" , 0); 

      plot.setXRange(0,Lt+2);
      plot.setYRange(low,high);

      std::stringstream plot_file;
      plot_file << plot_path << fname_base << ".ax"; 
      plot.sendToFile(plot_file.str());

    }

  // cant really reuse doPlot -- sigh
  template<> 
    void MakeAxisPlots::plot(const Fake3ptCorr<std::complex<double> >&C3pt) const
    {
      std::string  path =  SEMBLE::SEMBLEIO::getPath();

      // convert Q2 to a string format
      double Q2 = SEMBLE::toScalar(ENSEM::mean(C3pt.getKey().Q2));
      std::stringstream ssQ2;
      ssQ2 << Q2;
      std::string sQ2 = ssQ2.str();
      std::replace(sQ2.begin(),sQ2.end(),'.','p');
      std::replace(sQ2.begin(),sQ2.end(),'-','m');

      // set up the output paths
      path += std::string("FakeDataPlots/");
      SEMBLE::SEMBLEIO::makeDirectoryPath(path);

      std::string corr_path = path + std::string("Corrs/");
      std::string plot_path = path + std::string("Plots/"); 

      SEMBLE::SEMBLEIO::makeDirectoryPath(corr_path);
      SEMBLE::SEMBLEIO::makeDirectoryPath(plot_path); 

      // set up the base filename

      std::stringstream fname_base;
      fname_base << "C3pt_l" << C3pt.getKey().lorentz  << "_hf" << C3pt.getKey().hel_sink 
        << "_hi" << C3pt.getKey().hel_source << "_ID_" << C3pt.getKey().elemIDBase << "_pf" << C3pt.getKey().mom.momSink[0] << "_"
        << C3pt.getKey().mom.momSink[1] << "_" << C3pt.getKey().mom.momSink[2] << "_pi" << C3pt.getKey().mom.momSource[0] << "_"
        << C3pt.getKey().mom.momSource[1] << "_" << C3pt.getKey().mom.momSource[2] << "_Q2_" << sQ2;

      // get the loop variable bounds

      int nsource = C3pt.getNSource(); 
      int nsink = C3pt.getNSink();
      int Lt = C3pt.getLt();
      int ncfg = C3pt.getNCfg(); 

      // get some variables
      ENSEM::EnsemVectorComplex norm, normed_corr;
      ENSEM::EnsemVectorReal real,imag;
      norm.resize(ncfg);
      norm.resizeObs(Lt);
      norm = SEMBLE::toScalar(std::complex<double>(0,0));
      normed_corr = C3pt.get3pt();

      std::vector<double> time;

      // fill in the normalization
      for(int t = 0; t < Lt; ++t)
      {
        ENSEM::pokeObs(norm,getNorm<double>(C3pt.getKey(),t), t);
        time.push_back(double(t));
        ENSEM::pokeObs(normed_corr, ENSEM::peekObs(normed_corr,t)/ENSEM::peekObs(norm,t),t); 
      }

      // write the normalization factor as a function of time to a log
      std::stringstream normcorr; 
      normcorr << corr_path << fname_base.str() << "_normalization_by_time.jack";
      ENSEM::write(normcorr.str(),norm);


      real = ENSEM::real(normed_corr);
      imag = ENSEM::imag(normed_corr);

      std::vector<double> rmpe,rmme,impe,imme; 

      // get the initial values from the normalized correlator
      for(int t = 0; t < Lt; ++t)
      {
        double rmean = SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(real,t)));
        double rerr =  sqrt( SEMBLE::toScalar(ENSEM::variance(ENSEM::peekObs(real,t))));
        double imean =  SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(imag,t)));
        double ierr =  sqrt( SEMBLE::toScalar(ENSEM::variance(ENSEM::peekObs(imag,t))));

        rmme.push_back(rmean - rerr);
        rmpe.push_back(rmean + rerr); 
        imme.push_back(imean - ierr);
        impe.push_back(imean + ierr); 

      }

      AxisPlot rplot,iplot; 

      double rhigh = * std::max_element(rmpe.begin(),rmpe.end());
      double rlow = * std::min_element(rmme.begin(), rmme.end()); 

      double ihigh = * std::max_element(impe.begin(),impe.end());
      double ilow = * std::min_element(imme.begin(), imme.end()); 



      // do work!
      for(int source = 0; source < nsource; ++source)
        for(int sink = 0; sink < nsink; ++sink)
        {
          // record the ith_jth component of the sum that contributed to the correlator
          std::stringstream jack, jacknorm;
          jack << corr_path << fname_base.str() << "_sink" << sink << "_source" << source <<   ".jack";
          ENSEM::EnsemVectorComplex work = C3pt.get3ptComponent(source,sink); 
          ENSEM::write(jack.str(),work);

          std::vector<double> rline, iline;

          // normalize this piece of the correlator
          for(int t = 0; t < Lt; ++t)
          {
            ENSEM::pokeObs(work,ENSEM::peekObs(work,t)/ENSEM::peekObs(norm,t),t); 
            rline.push_back(SEMBLE::toScalar(ENSEM::mean(ENSEM::real(ENSEM::peekObs(work,t)))));
            iline.push_back(SEMBLE::toScalar(ENSEM::mean(ENSEM::imag(ENSEM::peekObs(work,t)))));
          }

          jacknorm << corr_path << fname_base.str() << "_sink" << sink << "_source" << source <<   "_normalized.jack";
          ENSEM::write(jacknorm.str(), work); 



          rplot.addLineData(time,rline, (sink*nsource + source) % 5 + 1); 
          iplot.addLineData(time,iline, (sink*nsource + source) % 5 + 1); 


          double tmpl = *std::min_element(rline.begin(),rline.end());
          rlow = std::min(rlow,tmpl);
          double tmph = *std::max_element(rline.begin(),rline.end());
          rhigh = std::max(rhigh,tmph); 

      //    std::cout << jack.str() << " " << sink << " " << source << " " << tmpl << " " << tmph << std::endl;

          tmpl = *std::min_element(iline.begin(),iline.end());
          ilow = std::min(ilow,tmpl);
          tmph = *std::max_element(iline.begin(),iline.end());
          ihigh = std::max(ihigh,tmph); 

        }


      double rheight = rhigh - rlow; 
      double runit = rheight / double(nsource * nsink);

      double iheight = ihigh - ilow; 
      double iunit = iheight / double(nsource * nsink);


      for(int source = 0; source < nsource; ++source)
        for(int sink = 0; sink < nsink; ++sink)
        {

          std::stringstream ss; 
          ss << sink << "_" << source;
          rplot.addLabel( Lt + 1, rlow + runit * (sink*nsource + source + 0.5), ss.str(), (sink*nsource + source) % 5 + 1, 0.8);
          iplot.addLabel( Lt + 1, rlow + iunit * (sink*nsource + source + 0.5), ss.str(), (sink*nsource + source) % 5 + 1, 0.8);

        }

      std::stringstream lab;
      lab << "sink_source";
      rplot.addLabel(Lt + 1, rlow + runit*0.1, lab.str() , 5, 0.8);
      iplot.addLabel(Lt + 1, ilow + iunit*0.1, lab.str() , 5, 0.8);
      
      std::stringstream gnu_info;
      gnu_info << "///dim=" << nsource*nsink;
      rplot.addLabel(Lt+1, rlow - runit*10 , gnu_info.str(), 1, 0.8); 
      iplot.addLabel(Lt+1, rlow - runit*10 , gnu_info.str(), 1, 0.8); 


      rplot.setXRange(0,Lt+2);
      rplot.setYRange(rlow,rhigh);

      iplot.setXRange(0,Lt+2);
      iplot.setYRange(ilow,ihigh);



      rplot.addEnsemData(real, "\\sq" , 0); 
      iplot.addEnsemData(imag,"\\sq", 0);

      std::stringstream rplot_file, iplot_file;
      rplot_file << plot_path << fname_base.str() <<  "_real.ax"; 
      iplot_file << plot_path << fname_base.str() <<  "_imag.ax"; 


      rplot.sendToFile(rplot_file.str());
      iplot.sendToFile(iplot_file.str());

    } 



} // namespace radmat
