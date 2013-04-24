/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Wed Apr 24 11:06:14 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/driver/radmat_single_q2_driver.h"
#include "semble/semble_semble.h"
#include "ensem/ensem.h"
#include <sstream>



namespace radmat
{



  namespace 
  {
    template<typename T>
      typename SEMBLE::PromoteEnsemVec<T>::Type
      get_ensem_row(const int row, const SEMBLE::SembleMatrix<T> &in)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type out;
        out.resize(in.getB());
        out.resizeObs(in.getM());
        for(int elem = 0; elem < in.getM(); ++elem)
          ENSEM::pokeObs(out,in.getEnsemElement(row,elem),elem);

        return out;
      }

  } // anonomyous 



  RadmatSingleQ2Driver::RadmatSingleQ2Driver(void)
  {
    init_false(); 
  }


  RadmatSingleQ2Driver::RadmatSingleQ2Driver(const RadmatSingleQ2Driver &o)
    : init_linear_system(o.init_linear_system) , init_fits(o.init_fits) ,
    linear_system(o.linear_system) , fit_across_time(o.fit_across_time)
  {  }


  RadmatSingleQ2Driver& RadmatSingleQ2Driver::operator=(const RadmatSingleQ2Driver &o)
  {
    if(this != &o)
    {
      init_linear_system = o.init_linear_system;
      init_fits = o.init_fits;
      linear_system = o.linear_system;
      fit_across_time = o.fit_across_time;
    }
    return *this; 
  }

  bool RadmatSingleQ2Driver::load_llsq(const ADAT::Handle<LLSQLatticeMultiData> &d,
      const std::string &soln_ID, const double pole_mass_squared)
  { 
    if(!!!linear_system.load_data(d))
      return false;

    if( (pole_mass_squared + SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) < 0.))
    {
      std::cout << "Killing Q2 = " << SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) 
        << "  q2tag = " << linear_system.qsq_label() << "  beacause pole_mass^2 + Q^2 < 0 "
        << "\nvalue is " << pole_mass_squared + SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) << std::endl;
      return false; 
    }


    std::cout << "Solving Q2 = " << linear_system.qsq_label() 
      << "  " << linear_system.peek_tags().begin()->mom_string() << std::endl;

    linear_system.solve_llsq(soln_ID);
    init_linear_system = true; 
    return true;
  }


  void RadmatSingleQ2Driver::fit_data(const ThreePointComparatorProps_t &fit_props)
  {
    check_exit_linear_system(); 
    SEMBLE::SembleMatrix<std::complex<double> > FF_of_t = linear_system.peek_FF(); 
    LLSQRet_ff_Q2Pack<std::complex<double> > tmp;

    for(int row = 0; row < FF_of_t.getN(); ++row)
      tmp.insert(row,get_ensem_row(row,FF_of_t));

    tmp.setQ2(linear_system.Q2()); 
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("t_ins_fits/"));
    fit_across_time.fit<std::complex<double> >(base_path() + std::string("t_ins_fits/"),
        tmp,
        fit_props);
    init_fits = true; 
  }


  void RadmatSingleQ2Driver::chisq_analysis(void)
  {
    check_exit_fits();
    std::cerr << __func__ 
      << ": warning, there is a stupid hardwire here that needs to be updated" << std::endl;

    int tlow = 10;
    int thigh = 24; 

    for(int fn = 0; fn < fit_across_time.nFF(); ++fn)
    {
      ADAT::Handle<FitThreePoint> some_fit = fit_across_time.getFit(fn); 
      std::string fit_name = some_fit->getFitName();
      std::cout << __func__ << " fit name " << fit_name << std::endl;
    }

    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("chisq/"));
    linear_system.chisq_analysis( fit_across_time.fetchFF().second,
        base_path() + std::string("chisq/") , tlow , thigh,1e-6);
  }


  ENSEM::EnsemReal RadmatSingleQ2Driver::Q2(void) const
  {
    check_exit_linear_system();
    return linear_system.Q2(); 
  }


  std::pair<ENSEM::EnsemReal, SEMBLE::SembleVector<double> > 
    RadmatSingleQ2Driver::fetchFF(void) const
    { 
      check_exit_fits(); 
      return fit_across_time.fetchFF(); 
    }


  std::string RadmatSingleQ2Driver::tags_at_this_Q2(void) const
  {
    check_exit_linear_system(); 
    double qq = SEMBLE::toScalar(ENSEM::mean(Q2())); 
    std::vector<LatticeMultiDataTag> tt = linear_system.peek_tags(); 
    std::vector<LatticeMultiDataTag>::const_iterator it;
    std::stringstream ss;

    for(it = tt.begin(); it != tt.end(); ++it)
      ss << qq << " (tag val = " <<  linear_system.qsq_label() << ") " << it->file_id << std::endl;

    return ss.str(); 
  }


  void RadmatSingleQ2Driver::dump_fits(void) 
  {
    check_exit_fits();
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("fit_logs/"));
    fit_across_time.writeFitLogs(base_path() + std::string("fit_logs/"));
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("component_fits/")); 
    fit_across_time.writeFitPlotsWithComponents(base_path() + std::string("component_fits/"));
  }


  void RadmatSingleQ2Driver::dump_llsq(void)
  {
    check_exit_linear_system();
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("llsq"));
    linear_system.dump_llsq(base_path() + std::string("llsq/")); 
  }



  void  RadmatSingleQ2Driver::check_exit(const bool b, const char *c) const 
  {
    if(!!!b)
    {
      std::cerr << __func__ << ": error: called by " << c << ", exiting." << std::endl;
      exit(1); 
    }
  }


  void RadmatSingleQ2Driver::init_false(void)
  { 
    init_linear_system = false;
    init_fits = false;
  }


  std::string RadmatSingleQ2Driver::base_path(void) const 
  {

    std::stringstream ss; 
    ss << SEMBLE::SEMBLEIO::getPath() << "Q2_" << linear_system.qsq_label() << "/";
    SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
    return ss.str(); 
  }


}



