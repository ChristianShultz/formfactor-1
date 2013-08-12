/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat_driver.cc

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Tue Jun 25 17:01:40 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/driver/radmat_driver.h"
#include "radmat/driver/radmat_driver_aux.h"

#include "radmat/utils/splash.h"
#include "radmat/load_data/load_fake_data.h"
#include "radmat/load_data/build_q2_packs.h"
#include "radmat/load_data/three_point.h"
#include "radmat/load_data/build_correlators.h"
#include "radmat/load_data/disconnected_graph_nuker.h"
#include "radmat/llsq/llsq_driver.h"
#include "radmat/llsq/llsq_q2_pack.h"
#include "radmat/driver/radmat_driver_props.h"

#include "jackFitter/plot.h"

#include "semble/semble_semble.h"

#include "io/adat_xmlio.h"

#include "ensem/ensem.h"

#include "adat/adat_stopwatch.h"

#include "hadron/ensem_filenames.h"

#include <complex>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>
#include <omp.h>

using namespace radmat;
using namespace ENSEM;
using namespace ADAT;
using namespace ADATIO;



#define LOAD_LLSQ_PARALLEL 
#define FIT_LLSQ_PARALLEL
#define CHISQ_ANALYSIS_PARALLEL


namespace radmat
{


  void RadmatDriver::run_program(const std::string &inifile)
  {
    init_false(); 
    read_xmlini(inifile);
    if(m_ini.maxThread > 1)
      omp_set_num_threads(m_ini.maxThread);
    build_correlators();
    solve_llsq();
    fit_ffs();
    do_chisq_analysis();
    make_FF_of_Q2_plots();
  }


  void RadmatDriver::xml_handler(const std::string &ini, const std::string &mode)
  {
    std::map<std::string, (RadmatDriver::*)(const std::string &)> handler; 
    std::map<std::string, (RadmatDriver::*)(const std::string &)>::const_iterator it; 
    handler["all"] = &RadmatDriver::build_xml; 
    handler["split"] = &RadmatDriver::build_xml_split_p2;
    
    it = handler.find(mode); 

    if (it != handler.end())
      this->(*(it->second))(ini); // get second which is a pointer, deref it 
    else
    {
      std::cerr << __PRETTY_FUNCTION__ << ": error mode, " << mode
        << " not recognized, try one of the following" << std::endl;
      for (it = handler.begin(); it != handler.end(); ++it)
        std::cout << it->first << std::endl; 
    }
  }


  void RadmatDriver::build_xml(const std::string &inifile)
  {
    read_xmlini(inifile);
    if(m_ini.maxThread > 1)
      omp_set_num_threads(m_ini.maxThread);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    keys = m_correlators.build_correlator_xml(m_ini.threePointIni); 

    

    ADATXML::XMLBufferWriter corrs;
    ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;

    bc.resize(keys.size()); 
    for(unsigned int i = 0; i < keys.size(); ++i)
      bc[i] = keys[i];

    write(corrs,"NPointList",bc);

    std::ofstream out("npt.list.xml");
    corrs.print(out);
    out.close();

    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it; 

    out.open("npt.ensemFileNames.list"); 
    for(it = keys.begin(); it != keys.end(); ++it)
      out << Hadron::ensemFileName(*it) << "\n";
    out.close(); 
  }


  void RadmatDriver::nuke_graph(const std::string &inifile, 
      const std::string &graph_db,
      const std::string &nuke_xml_out)
  {

    read_xmlini(inifile);
    if(m_ini.maxThread > 1)
      omp_set_num_threads(m_ini.maxThread);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    keys = m_correlators.build_correlator_xml(m_ini.threePointIni); 

    DisconnectedGraphNuker n;
    n.find_nukes(keys,graph_db); 
    n.dump_nukes(nuke_xml_out); 
  }


  void RadmatDriver::build_stub_xml(const std::string &inifile)
  {
    read_xmlini(inifile);
    if(m_ini.maxThread > 1)
      omp_set_num_threads(m_ini.maxThread);

    std::vector<Hadron::KeyHadronNPartNPtCorr_t> keys;
    keys = m_correlators.build_correlator_xml(m_ini.threePointIni); 


    stubify(keys);

    ADATXML::XMLBufferWriter corrs;
    ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;

    bc.resize(keys.size()); 
    for(unsigned int i = 0; i < keys.size(); ++i)
      bc[i] = keys[i];

    write(corrs,"NPointList",bc);

    std::ofstream out("npt.list.xml");
    corrs.print(out);
    out.close();

    std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it; 

    out.open("npt.ensemFileNames.list"); 
    for(it = keys.begin(); it != keys.end(); ++it)
      out << Hadron::ensemFileName(*it) << "\n";
    out.close(); 
 
  }



  void RadmatDriver::init_false(void)
  {
    read_ini = false;
    built_correlators = false;
    init_llsq = false; 
    solved_llsq = false;
    fit_formfacs = false;
    chisq_analysis = false; 
  }

  void RadmatDriver::check_exit(const bool &b, const char *c) const 
  {
    if(!!!b)
    {
      std::cerr << __func__ << ": error: called by " << c << ", exiting." << std::endl;
      exit(1); 
    }
  }


  void RadmatDriver::read_xmlini(const std::string &xmlini)
  {


    try
    {
      XMLReader xml(xmlini);
      read(xml,"/DriverProps",m_ini);
    }
    catch(std::exception &e)
    {
      std::cout << "std exception: " << e.what();
    }
    catch(std::string &e)
    {
      std::cout << __func__ << ": ERROR: can't read xmlinifile ("
        << xmlini << "): " << e << std::endl;
      exit(1);
    }
    catch(...)
    {
      SPLASH("An error occured while loading the fake data inifile");
      exit(1);
    }

    read_ini = true; 
  }



  void RadmatDriver::build_correlators(void)
  {
    check_exit_ini(); 

    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 

    std::cout << "Loading correlators and inverting subduction.. " << std::endl; 

    multi_lattice_data = m_correlators.build_multi_correlators(m_ini.threePointIni);

    my_stopwatch.stop(); 
    std::cout << "Loading correlators and inverting subduction took " 
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    built_correlators = true; 
  }


  void RadmatDriver::solve_llsq(void)
  {
    check_exit_corrs(); 

    int idx, sz = multi_lattice_data.size(); 
    std::string soln_ID = std::string ("SVDNonSquare");

    if(sz == 0)
    {
      std::cerr << __func__ << ": error nothing to solve" << std::endl;
      exit(1); 
    }

    linear_systems_of_Q2.resize(multi_lattice_data.size());

    std::cout << "Solving LLSQ.. " << std::endl;

    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 


    good_qs.resize(sz,false); 

    std::cout << __func__ << ": sz = " << sz << std::endl;


#ifdef LOAD_LLSQ_PARALLEL 

#pragma omp parallel for shared(idx)  schedule(dynamic,1)

#endif    
    // POSSIBLE PARALLEL HERE
    for(idx =0; idx < sz; ++idx)
      good_qs[idx] =  linear_systems_of_Q2[idx].load_llsq(multi_lattice_data[idx],m_ini.poleMass);
    // END PARALLEL

#pragma omp barrier

    int ngood(0);
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        ++ngood; 

    init_llsq = true;

    std::cout << __func__ << ": " << ngood << " good Q^2 points out of " << sz << std::endl;

    // print the list here in case the solver flakes we can easily determine where it went wrong
    print_Q2_list(); 

#ifdef LOAD_LLSQ_PARALLEL 

#pragma omp parallel for shared(idx)  schedule(dynamic,1)

#endif    
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].solve_llsq(soln_ID); 

#pragma omp barrier

    my_stopwatch.stop();
    std::cout << "Solving LLSQ took "     
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    solved_llsq = true;

    // probably shouldn't parallel here since the file system will get pissed off
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])                        // can only print the successful guys
        linear_systems_of_Q2[idx].dump_llsq(); 
  }

  void RadmatDriver::fit_ffs(void)
  {
    check_exit_llsq(); 

    int idx, sz = multi_lattice_data.size(); 

    std::cout << "Fitting FF(t_ins) " << std::endl;

    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 


#ifdef FIT_LLSQ_PARALLEL

#pragma omp parallel for shared(idx)  schedule(dynamic,1)

#endif 
    // POSSIBLE PARALLEL HERE
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].fit_data(m_ini.threePointComparatorProps);

#pragma omp barrier

    my_stopwatch.stop();
    std::cout << "Fitting took " 
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    fit_formfacs = true; 

    // probably shouldn't parallel here since the file system will get pissed off
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])                        // can only print the successful guys
        linear_systems_of_Q2[idx].dump_fits(); 
  }


  void RadmatDriver::do_chisq_analysis(void)
  {
    check_exit_fit();

    int idx, sz = multi_lattice_data.size(); 

    std::cout << "chisq_analysis" << std::endl;

    Util::StopWatch my_stopwatch; 
    my_stopwatch.start(); 

#ifdef CHISQ_ANALYSIS_PARALLEL

#pragma omp parallel for shared(sz) schedule(dynamic,1)

#endif

    // POSSIBLE PARALLEL HERE
    for(idx = 0; idx < sz; ++idx)
      if(good_qs[idx])
        linear_systems_of_Q2[idx].chisq_analysis();

#pragma omp barrier

    my_stopwatch.stop();
    std::cout << "chisq_analysis took " 
      << my_stopwatch.getTimeInSeconds() << " seconds " << std::endl;

    chisq_analysis = true; 

  }


  void RadmatDriver::make_FF_of_Q2_plots(void)
  {
    std::string pth = SEMBLE::SEMBLEIO::getPath();
    std::stringstream path;
    path << pth << "FF_of_Q2/";
    SEMBLE::SEMBLEIO::makeDirectoryPath(path.str()); 


    int nff = 0; 
    int nQs = 0;
    int ncfg = 0; 


    for(unsigned int allQ = 0; allQ < linear_systems_of_Q2.size(); ++allQ)
    {
      if(linear_systems_of_Q2[allQ].check_fits()) 
      {
        std::cout << "allQ = " << allQ << std::endl; 
        nff = (linear_systems_of_Q2[allQ].fetchFF()).second.getN();
        nQs = linear_systems_of_Q2.size(); 
        ncfg =  linear_systems_of_Q2[allQ].Q2().size(); 
        break; // only do this once but don't assume any ordering
      }
      if(allQ == linear_systems_of_Q2.size() -1)
      {
        std::cerr <<__func__ << ": the linear system was empty. exiting" << std::endl;
        exit(1);
      }
    }

    std::vector<double> q2s;

    for(int Q = 0; Q < nQs; ++Q)
      if(good_qs[Q])
        q2s.push_back(SEMBLE::toScalar(ENSEM::mean(linear_systems_of_Q2[Q].Q2()))); 

    for(int ff; ff < nff; ++ff)
    {
      ENSEM::EnsemVectorReal FF;
      FF.resize(ncfg);
      FF.resizeObs(nQs);

      std::vector<ENSEM::EnsemReal> tmp;  

      // need to make sure to match up allQ vs good Q so we get the data
      // at the correct position
      for(int allQ = 0; allQ < nQs; ++allQ)
        if(good_qs[allQ])
          tmp.push_back( (linear_systems_of_Q2[allQ].fetchFF()).second.getEnsemElement(ff) ); 

      for(unsigned int goodQ = 0; goodQ < q2s.size(); ++goodQ)
        ENSEM::pokeObs(FF, tmp[goodQ],goodQ); 

      AxisPlot plt; 
      plt.addEnsemData(q2s,FF,"\\sq",1);

      std::stringstream ss; 
      ss << path.str() << "FF_" << ff << ".ax"; 
      plt.sendToFile(ss.str()); 
    }

  }


  void RadmatDriver::print_Q2_list(void) 
  {
    check_exit_init_llsq(); 
    std::stringstream ss;
    std::vector<RadmatSingleQ2Driver>::const_iterator it; 
    std::string pth = SEMBLE::SEMBLEIO::getPath() + std::string("Q2_to_mat_elems.txt"); 
    std::ofstream out(pth.c_str());
    std::string delim("--------------------------------\n"); 
    for(it = linear_systems_of_Q2.begin(); it != linear_systems_of_Q2.end(); ++it)
      if(it->check_linear_system())
        out << it->tags_at_this_Q2() << delim; 
    out.close();
  }











}


#undef LOAD_LLSQ_PARALLEL 
#undef FIT_LLSQ_PARALLEL
#undef CHISQ_ANALYSIS_PARALLEL



