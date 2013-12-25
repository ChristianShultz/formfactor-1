/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Mon 23 Dec 2013 09:29:37 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/driver/radmat_single_q2_driver.h"
#include "radmat/llsq/llsq_multi_data_serialize.h"
#include "radmat/ff/lorentzff_canonical_rotations.h"
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
  
    void write_jackfile_fit_report ( const TinsFitter &fits, const std::string &bpth)
    {
      for (int ff = 0; ff < fits.nFF(); ++ff)
      {
        std::stringstream ss;
        ss << bpth << "FF_" << ff << ".jack";
        ENSEM::write(ss.str(),fits.getFF(ff));  
      }

      std::stringstream ss; 
      ss << bpth << "Q2.jack";
      ENSEM::write(ss.str(),fits.getQ2()); 
    }

    void write_single_jackfile_fit_report ( const TinsFitter &fits, const std::string &bpth, const int ff)
    {
      std::stringstream ss;
      ss << bpth << "FF_" << ff << ".jack";
      ENSEM::write(ss.str(),fits.getFF(ff));  

      std::stringstream s; 
      s << bpth << "Q2.jack";
      ENSEM::write(s.str(),fits.getQ2()); 
    }

  } // anonomyous 



  RadmatSingleQ2Driver::RadmatSingleQ2Driver(void)
  {
    init_false(); 
  }


  RadmatSingleQ2Driver::RadmatSingleQ2Driver(const RadmatSingleQ2Driver &o)
    : init_linear_system(o.init_linear_system) , init_fits(o.init_fits) ,
    linear_system(o.linear_system) , fit_across_time(o.fit_across_time)
  { rot_id = std::string(""); }


  RadmatSingleQ2Driver& RadmatSingleQ2Driver::operator=(const RadmatSingleQ2Driver &o)
  {
    if(this != &o)
    {
      init_linear_system = o.init_linear_system;
      init_fits = o.init_fits;
      linear_system = o.linear_system;
      fit_across_time = o.fit_across_time;
      rot_id = o.rot_id; 
    }
    return *this; 
  }

  bool RadmatSingleQ2Driver::load_llsq(const rHandle<LLSQLatticeMultiData> &d, 
      const double pole_mass_squared,
      const double tolerance)
  {

    if(!!!linear_system.load_data(d,tolerance))
      return false;

    if(linear_system.peek_tags().empty())
    {
      std::cerr << __func__ << ": warning, no tags" << std::endl; 
      return false;
    }

    if( (pole_mass_squared + SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) < 0.))
    {
      std::cout << "Killing Q2 = " << SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) 
        << "  q2tag = " << linear_system.qsq_label() << "  beacause pole_mass^2 + Q^2 < 0 "
        << "\nvalue is " << pole_mass_squared + SEMBLE::toScalar(ENSEM::mean(linear_system.Q2())) << std::endl;
      return false; 
    }

    // append the label
    append_rotation_group_label(rotation_group_label()); 
    

    init_linear_system = true; 

    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("llsq"));
    linear_system.dump_llsq_lattice(base_path() + std::string("llsq/"));
    linear_system.save_llsq_state( base_path() + std::string("llsq/") ); 
    return true;
  }

  bool RadmatSingleQ2Driver::load_llsq(const rHandle<LLSQLatticeMultiData> &d, 
      const double tolerance)
  {

    if(!!!linear_system.load_data(d,tolerance))
      return false;

    if(linear_system.peek_tags().empty())
    {
      std::cerr << __func__ << ": warning, no tags" << std::endl; 
      return false;
    }

    init_linear_system = true; 

    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("llsq"));
    linear_system.dump_llsq_lattice(base_path() + std::string("llsq/"));
    return true;
  }

  void RadmatSingleQ2Driver::solve_llsq(const std::string &soln_ID)
  {
    check_exit_linear_system();
    std::cout << "Solving Q2 = " << linear_system.qsq_label() 
      << "  " << linear_system.peek_tags().begin()->mom_string() << std::endl;

    linear_system.solve_llsq(soln_ID);
    init_solved_llsq = true; 
  }

  void RadmatSingleQ2Driver::save_llsq_state(void) const
  {
    check_exit_solved_llsq(); 
    linear_system.save_ff_state( base_path() ); 
  }

  void RadmatSingleQ2Driver::save_ff_of_t(void) const
  {
    check_exit_solved_llsq(); 
    SEMBLE::SembleMatrix<std::complex<double> > FF_of_t = linear_system.peek_FF(); 
    std::string path;
    path = base_path() + std::string("ff_of_t/");
    SEMBLE::SEMBLEIO::makeDirectoryPath( path ); 

    for(int row = 0; row < FF_of_t.getN(); ++row)
    {
      std::stringstream id;
      id << path << "unphasedFF_" << row << ".jack"; 
      ENSEM::write( id.str() , get_ensem_row(row,FF_of_t)); 
    }

  }

  void RadmatSingleQ2Driver::fit_data(const ThreePointComparatorProps_t &fit_props, 
      const int tsrc,
      const int tsnk)
  {
    check_exit_linear_system();
    check_exit_solved_llsq(); 
    SEMBLE::SembleMatrix<std::complex<double> > FF_of_t = linear_system.peek_FF(); 
    LLSQRet_ff_Q2Pack<std::complex<double> > tmp;

    for(int row = 0; row < FF_of_t.getN(); ++row)
      tmp.insert(row,get_ensem_row(row,FF_of_t));

    tmp.setQ2(linear_system.Q2()); 
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("t_ins_fits/"));
    fit_across_time.fit<std::complex<double> >(base_path() + std::string("t_ins_fits/"),
        tmp,
        fit_props,
        tsrc,
        tsnk);
    init_fits = true; 
  }

  void 
    RadmatSingleQ2Driver::fit_and_dump_single_ffs(const ThreePointComparatorProps_t &fit_props, 
        const SEMBLE::SembleMatrix<std::complex<double> > &FF_of_t,
        const ENSEM::EnsemReal &Q2,
        const int tsrc,
        const int tsnk,
        const int ff,
        const int ff_max) const
    {
      // snarf (RGEism)  ffs
      LLSQRet_ff_Q2Pack<std::complex<double> > tmp;

      // build a "fake" pack
      tmp.insert(ff,get_ensem_row(ff,FF_of_t));
      tmp.setQ2(Q2); 

      // i/o junk 
      std::string base_p = std::string("REFIT/");
      SEMBLE::SEMBLEIO::makeDirectoryPath(base_p); 
      SEMBLE::SEMBLEIO::makeDirectoryPath(base_p + std::string("t_ins_fits/"));

      // run a fit on this ff 
      TinsFitter local_fit_across_time; 
      local_fit_across_time.single_fit<std::complex<double> >(base_p + std::string("t_ins_fits/"),
          tmp,
          ff_max, 
          fit_props,
          tsrc,
          tsnk);


      // dump results
      SEMBLE::SEMBLEIO::makeDirectoryPath(base_p + std::string("fit_logs/"));
      local_fit_across_time.writeFitLogs(base_p + std::string("fit_logs/"));
      SEMBLE::SEMBLEIO::makeDirectoryPath(base_p + std::string("component_fits/")); 
      local_fit_across_time.writeFitPlotsWithComponents(base_p + std::string("component_fits/"));
      SEMBLE::SEMBLEIO::makeDirectoryPath(base_p + std::string("jack_files/")); 
      write_single_jackfile_fit_report ( local_fit_across_time, base_p + std::string("jack_files/"), ff);
    }

  void RadmatSingleQ2Driver::chisq_analysis(const int tlows, const int thighs)
  {
    check_exit_fits();

    int tlow = tlows; 
    int thigh = thighs; 

    for(int fn = 0; fn < fit_across_time.nFF(); ++fn)
    {
      rHandle<FitThreePoint> some_fit = fit_across_time.getFit(fn); 
      if ( some_fit->tlow() > tlow) 
        tlow = some_fit->tlow(); 

      if ( some_fit->thigh() < thigh)
        thigh = some_fit->thigh();

    }

    if ( thigh < tlow ) 
    {
      std::cout << __func__ << ": range [" << tlow << "," << thigh << "]" << std::endl; 
      std::cout << "thigh < tlow, unable to continue " << std::endl; 
      return ; 
    }

    std::cout << __func__ << ": using range [" << tlow << "," << thigh << "]" << std::endl; 

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
      ss << qq << " (tag val = " <<  linear_system.qsq_sort() << ") " << it->file_id << std::endl;

    return ss.str(); 
  }

  std::string RadmatSingleQ2Driver::rotation_group_label(void) const
  {
    // check_exit_linear_system(); 
    std::vector<LatticeMultiDataTag> tt = linear_system.peek_tags(); 
    std::vector<LatticeMultiDataTag>::const_iterator it;
  
    if( tt.empty() )
      return std::string(); 

    it = tt.begin(); 
    std::string chk = radmat::LatticeRotationEnv::rotation_group_label(it->p_f,it->p_i); 
    for(it = tt.begin(); it != tt.end(); ++it)
      if( chk != radmat::LatticeRotationEnv::rotation_group_label(it->p_f,it->p_i) )
      {
        std::cout << __func__ << ": returning empty string since more "
          << " than one rotation group is present, unless one of the "
          << "particles is spin zero this is rather problematic for you "
          << std::endl;
       return std::string();  
      }
    return chk; 
  }

  void RadmatSingleQ2Driver::dump_fits(void) 
  {
    check_exit_fits();
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("fit_logs/"));
    fit_across_time.writeFitLogs(base_path() + std::string("fit_logs/"));
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("component_fits/")); 
    fit_across_time.writeFitPlotsWithComponents(base_path() + std::string("component_fits/"));
    SEMBLE::SEMBLEIO::makeDirectoryPath(base_path() + std::string("jack_files/")); 
    write_jackfile_fit_report ( fit_across_time, base_path() + std::string("jack_files/"));
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
    rot_id = std::string(""); 
  }


  std::string RadmatSingleQ2Driver::base_path(void) const
  {
    std::stringstream ss; 
    ss << SEMBLE::SEMBLEIO::getPath() << "Q2_" << linear_system.qsq_sort();

    // this is an ugly naming convention 
    if( rot_id != std::string("") )
    {
      SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
      ss << "/" << rot_id; 
    }

    ss << "/";
    SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
    return ss.str(); 
  }


}



