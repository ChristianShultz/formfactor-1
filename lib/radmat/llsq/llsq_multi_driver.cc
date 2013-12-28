/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_multi_driver.cc

 * Purpose :

 * Creation Date : 22-02-2013

 * Last Modified : Fri Dec 27 19:31:06 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "llsq_multi_driver.h"
#include "llsq_solvers.h"
#include "llsq_chisq.h"
#include "llsq_multi_data_serialize.h"
#include "llsq_solution.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/ff/ff_gen_llsq_row.h"
#include "ensem/ensem.h"
#include <complex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>

#include "adat/handle.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/plot.h"
// #define DEBUG_AT_MAKE_MOM_INV_TAGS



namespace radmat
{

  namespace 
  {
    bool is_singular(const SEMBLE::SembleMatrix<std::complex<double> > &M)
    {
      itpp::Mat<std::complex<double> > mm( M.mean() ); 
      int bound = mm.cols();
      itpp::Vec<double> s = itpp::svd(mm * itpp::hermitian_transpose(mm));
      if ( s(bound - 1) / s(0) < 1e-4 )
      {
        std::cout << __func__ << bound << std::endl;
        std::cout << __func__ << ": M is singular, sings " << s << std::endl;
        std::cout << __func__ << ": M \n" << mm << std::endl;
        return true; 
      }
      return false; 
    }



    void init_dim(SEMBLE::SembleMatrix<std::complex<double> > &to_init, 
        const SEMBLE::SembleVector<std::complex<double> > &first_row)
    {
      SEMBLE::SembleMatrix<std::complex<double> > foo(first_row.getB(),1,first_row.getN());
      foo.zeros(); 
      for(int elem = 0; elem < first_row.getN(); ++elem)
        foo.loadEnsemElement(0,elem,first_row.getEnsemElement(elem)); 

      to_init = foo; 
    }

    int remap_jmu_4_to_0(const int &i)
    {
      if ( i == 4 ) 
        return 0; 
      else if ( i < 4 && i > 0 ) 
        return i;
      else
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ 
          << "error index " << i << " is out of range" 
          << std::endl; 
        exit(1337);
      }
    }

    SemblePInv makeMomInvariants(const LatticeMultiDataTag &t)
    {
#ifdef DEBUG_AT_MAKE_MOM_INV_TAGS 

      std::cout << __func__ << ": debuggin on" << std::endl;
      t.print_me();
      std::cout << "mom_string = " << t.mom_string() << std::endl;
      std::cout << "E_string = " << t.E_string() << std::endl;
#endif

      // need to scope for this to compile
      return radmat::makeMomInvariants(t.E_f,
          t.E_i,
          t.p_f,
          t.p_i,
          t.mom_fac);
    }


    std::string sort_string(const LatticeMultiDataTag &t)
    {
      std::stringstream ss; 
      ss << t.p_f[0] <<  t.p_f[1] <<  t.p_f[2] <<  t.p_i[0] <<  t.p_i[1] <<  t.p_i[2] << t.jmu; 
      return ss.str(); 
    }

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


    template<typename T>
      void my_writer_mean(const std::string &fname, const SEMBLE::SembleMatrix<T> &in)
      {
        itpp::Mat<T> mean = itpp::round_to_zero( in.mean() , 1e-5 );
        std::ofstream out(fname.c_str());
        out << mean ;
        out.close(); 
      }

    template <typename T>
      std::string number_to_string ( T Number )
      {
        std::stringstream ss;
        ss << Number;
        return ss.str();
      }

    template<typename T>
      void my_writer_rows(const std::string &fname, const SEMBLE::SembleMatrix<T> &in)
      {
        const int nr = in.getN();

        for( int row = 0; row < nr; ++ row)
        {
          ENSEM::write(fname + number_to_string(row) + std::string(".jack"), 
              get_ensem_row(row,in)); 
        }
      }

    void my_writer_cont_expr(const std::string &fname, const std::vector<LatticeMultiDataTag> &tt)
    {
      std::ofstream out(fname.c_str());
      for(unsigned int idx = 0; idx < tt.size(); ++idx)
        out << idx << " " << tt[idx].file_id << std::endl;
      out.close(); 
    }

    template<typename T>
      void my_writer(const std::string &fname, const T &t)
      {
        std::ofstream out(fname.c_str());
        out << t;
        out.close(); 
      }


    // take the fit form factors and compute the chisq of the fit across time 
    std::string chisq_per_t(const SEMBLE::SembleVector<std::complex<double> > &FF,
        const rHandle<LLSQLatticeMultiData> &lattice_data,
        const SEMBLE::SembleMatrix<std::complex<double> > &K,
        const double tol)
    {
      SEMBLE::SembleMatrix<std::complex<double> > data = lattice_data->data(); 
      SEMBLE::SembleVector<std::complex<double> > thy = K * FF; 
      std::stringstream ss; 
      ss << "#t chisq/DoF nDoF \n";  // think # makes gnuplot ignore the line..

      const int Lt = data.getM();

      for(int t = 0; t < Lt; ++t)
      {
        std::pair<double,int> chisq = ChisqAndDoF(data.getCol(t),thy,tol);
        ss << t << " " << chisq.first / double(chisq.second) << " " << chisq.second << std::endl;
      }

      return ss.str(); 
    }


    std::string chisq_per_data_of_fit_range(const SEMBLE::SembleVector<std::complex<double> > &FF,
        const rHandle<LLSQLatticeMultiData> &lattice_data,
        const SEMBLE::SembleMatrix<std::complex<double> > &K,
        const double tol,
        const int tlow,
        const int thigh)
    {
      const int Lt = thigh - tlow; 
      SEMBLE::SembleMatrix<std::complex<double> > data = lattice_data->data();
      std::vector<LatticeMultiDataTag> tags = lattice_data->tags(); 
      SEMBLE::SembleVector<std::complex<double> > thy = K * FF; 
      const int sz = tags.size(); 
      std::stringstream ss;
      ss << "#id chisq/DoF nDoF \n";

      for(int corr = 0; corr < sz; ++corr)
      {
        SEMBLE::SembleVector<std::complex<double> > work_thy,work_data(data.getB(),Lt); 
        work_thy = work_data; 
        for(int t = 0; t < Lt; ++t)
        {
          work_thy.loadEnsemElement(t,thy.getEnsemElement(corr));
          work_data.loadEnsemElement(t,data.getEnsemElement(corr,t+tlow));
        }
        std::pair<double,int> chisq = ChisqAndDoF(work_data,work_thy,tol);
        ss << tags[corr].file_id << " " << chisq.first / double(chisq.second) 
          << " " << chisq.second << std::endl;
      }

      return ss.str(); 
    }


    std::string chisq_of_system_of_fit_range(const SEMBLE::SembleVector<std::complex<double> > &FF,
        const rHandle<LLSQLatticeMultiData> &lattice_data,
        const SEMBLE::SembleMatrix<std::complex<double> > &K,
        const double tol,
        const int tlow,
        const int thigh)
    {
      const int Lt = thigh - tlow; 
      SEMBLE::SembleMatrix<std::complex<double> > data = lattice_data->data();
      SEMBLE::SembleVector<std::complex<double> > thy = K * FF; 
      std::stringstream ss;
      const int ncorr = data.getN(); 
      const int rankV = ncorr * Lt; 
      SEMBLE::SembleVector<std::complex<double> > wdata(data.getB(),rankV),wthy(data.getB(),rankV); 

      // construct vectors for chisq
      for(int corr = 0; corr < ncorr; ++corr)      
        for(int t = 0; t < Lt; ++t)
        {
          wdata.loadEnsemElement(corr*Lt + t, data.getEnsemElement(corr,t+tlow));
          wthy.loadEnsemElement(corr*Lt + t, thy.getEnsemElement(corr)); 
        }

      std::pair<double,int> chisq = ChisqAndDoF(wdata,wthy,tol); 
      ss<< "tlow = " << tlow << " thigh = " << thigh << " tol = " << tol << std::endl;
      ss << "chisq = " << chisq.first / double(chisq.second) << " nDoF = " << chisq.second;
      return ss.str(); 
    }



  } // anonomyous 





  LLSQMultiDriver_t::LLSQMultiDriver_t(void) 
  {
    init_false(); 
  }


  bool LLSQMultiDriver_t::load_data(const rHandle<LLSQLatticeMultiData> &d,
      const double tolerance)
  {
    init_false();
    lattice_data = d;
    POW2_ASSERT(&*d);
    init_lat = true; 

    sort_data(); 
    bool success = run_zero_filter(tolerance);
    return success; 
  }


  void LLSQMultiDriver_t::splash_tags(void) const
  {
    lattice_data->splash_tags(); 
  }


  void LLSQMultiDriver_t::sort_data(void)
  {
    rHandle<LLSQLatticeMultiData> sorted_data(new LLSQLatticeMultiData); 
    std::vector<LatticeMultiDataTag> tags = lattice_data->tags(); 
    std::map<std::string,std::vector<int> > elems;
    std::map<std::string,std::vector<int> >::iterator it; 

    unsigned int sz = tags.size(); 
    for(unsigned int elem = 0; elem < sz; ++elem)
    {
      std::string tt = sort_string(tags[elem]);
      it = elems.find(tt);
      if(it == elems.end())
      {
        elems[tt] = std::vector<int>(1,elem); 
      }
      else
      {
        it->second.push_back(elem); 
      }
    }


    for(it = elems.begin(); it != elems.end(); ++it)
    {
      std::vector<int>::const_iterator sorted; 
      for(sorted = it->second.begin(); sorted != it->second.end(); ++sorted)
        sorted_data->append_row_semble(lattice_data->get_row_semble(*sorted),tags[*sorted]);
    }

    lattice_data = sorted_data; 
  }



  bool LLSQMultiDriver_t::run_zero_filter(const double tolerance)
  {
    check_exit_lat(); 

    // splash_tags();

    rHandle<LLSQLatticeMultiData> non_zero_data(new LLSQLatticeMultiData);
    rHandle<LLSQLatticeMultiData> check_singular( new LLSQLatticeMultiData); 
    std::vector<LatticeMultiDataTag> old_tags;
    SEMBLE::SembleMatrix<std::complex<double> > Junk; 
    SEMBLE::SembleVector<std::complex<double> > Zero; 
    old_tags = lattice_data->tags(); 

    const unsigned int sz = old_tags.size(); 

    if ( sz == 0 ) 
    {
      std::cerr << __func__ << "no tags?? " << std::endl; 
      return false; 
    }


    ffKinematicFactors_t<std::complex<double> >
      KJunk(FormFactorDecompositionFactoryEnv::callFactory(
            old_tags.begin()->mat_elem_id)); 

    Junk = KJunk.genFactors(makeMomInvariants(*old_tags.begin())); 
    Zero = Junk.getRow(0);
    Zero.zeros(); 

    for(unsigned int elem = 0; elem < sz; ++elem)
    {
      ffKinematicFactors_t<std::complex<double> >
        KK(FormFactorDecompositionFactoryEnv::callFactory(old_tags[elem].mat_elem_id)); 
      SEMBLE::SembleMatrix<std::complex<double> > workM;
      SEMBLE::SembleVector<std::complex<double> > workV;

      workM = KK.genFactors(makeMomInvariants(old_tags[elem]));
     // std::cout << __FILE__ << __func__ << workM.mean() << std::endl;
      workV = SEMBLE::round_to_zero(
          workM.getRow(
            remap_jmu_4_to_0(old_tags[elem].jmu)), tolerance);

      if(workV == Zero)
      {
        zeroed_data.append_row_semble(
            lattice_data->get_row_semble(elem),
            old_tags[elem]);    
      }  
      else
      {
        non_zero_data->append_row_semble(
            lattice_data->get_row_semble(elem),
            old_tags[elem]);
        check_singular->append_row_semble(workV, old_tags[elem]); 
      }
    }

    lattice_data = non_zero_data;

    // warn that we are killing this data point 
    if(non_zero_data->nrows() < KJunk.nFacs())
    {
      std::cout << __func__ << ": not enough data points to solve the llsq" << std::endl;
      std::cout << "passed in " << sz << " elements of which " << zeroed_data.nrows() 
        << " failed the zero test, needed " << KJunk.nFacs() << " elems, had " 
        << non_zero_data->nrows() << "elements " << std::endl;
    }

    // check if the matrix is singular 
    if( lattice_data->nrows() >= KJunk.nFacs() )
      if ( is_singular( check_singular->data() ) )
        return false; 

    return (non_zero_data->nrows() >= KJunk.nFacs());
  }




  void LLSQMultiDriver_t::solve_llsq(const std::string &soln_ID)
  {
    check_exit_lat();

    rHandle<LLSQBaseSolver_t<std::complex<double> > >
      foo = LLSQSolverFactoryEnv::callFactory(soln_ID);

    bool success(true); 

    if(foo->invertable())
      success &= solve_fast(soln_ID); 
    else
      success &= solve_slow(soln_ID);

    check_exit(success,__func__);
  }


  void LLSQMultiDriver_t::chisq_analysis(const SEMBLE::SembleVector<double> &ff,
      const std::string &path,
      const int tlow,
      const int thigh,
      const double tol)
  {
    check_exit_lat();
    check_exit_K();
    check_exit_Kinv();
    check_exit_FF(); 



    // itpp/semble dont like mixed operations
    SEMBLE::SembleVector<std::complex<double> > ff_cmplx(ff.getB(),ff.getN());
    const int B = ff.getB();
    const int N = ff.getN();

    for(int bin = 0; bin < B; ++bin)
      for(int elem = 0; elem < N; ++elem)
        ff_cmplx.setElement(bin,elem,std::complex<double>(ff[bin][elem],0.));

    std::string chisq_per_tt = chisq_per_t(ff_cmplx,lattice_data,K,tol);
    std::string chisq_per_data = chisq_per_data_of_fit_range(ff_cmplx,
        lattice_data,K,tol,tlow,thigh);
    std::string chisq_of_sys = chisq_of_system_of_fit_range(ff_cmplx,
        lattice_data,K,tol,tlow,thigh); 

    std::stringstream ss,per_t,per_data,sys;
    ss << path; 
    per_t << ss.str() << "chisq_per_t.txt";
    per_data << ss.str() << "chisq_per_data.txt";
    sys << ss.str() << "chisq_of_sys_of_fit_range.txt"; 

    my_writer(per_t.str(),chisq_per_tt);
    my_writer(per_data.str(),chisq_per_data);
    my_writer(sys.str(),chisq_of_sys); 
  }


  void LLSQMultiDriver_t::dump_llsq_lattice(const std::string &path)
  {
    check_exit_lat();
    std::stringstream b,cont,zero,db;
    cont << path << "row_index_to_continuum_elem.txt";
    zero << path << "zeroed_matrix_elems";
    b << path << "lattice_mat_elem_";

    my_writer_rows(b.str(), lattice_data->data());
    my_writer_cont_expr(cont.str(), lattice_data->tags());

    SEMBLE::SEMBLEIO::makeDirectoryPath(zero.str()); 
    zero << "/zeroed_mat_elem_"; 
    my_writer_rows(zero.str(),zeroed_data.data()); 
    zero.str(std::string()); 
    zero << path << "zeroed_matrix_elems/row_index_to_continuum_elem.txt"; 
    my_writer_cont_expr(zero.str(),zeroed_data.tags()); 
  }



  void LLSQMultiDriver_t::dump_llsq(const std::string &path)
  {
    check_exit_K();
    check_exit_Kinv();
    check_exit_FF(); 


    std::stringstream ss, A,Ainv, x; 
    std::stringstream unity; 
    ss << path;
    A << ss.str() << "K.mean";
    Ainv << ss.str() << "Kinv.mean"; 
    unity << ss.str() <<  "Kinv_x_K.mean"; 
    x << ss.str() << "ff_";

    my_writer_mean(A.str(),K);
    my_writer_mean(Ainv.str(),Kinv);
    my_writer_mean(unity.str(),Kinv*K); 
    my_writer_rows(x.str(), FF_t); 

    std::string pth = path + std::string("solver.log"); 
    std::ofstream out( pth.c_str() ); 
    out << solver_log; 
    out.close(); 
  }



  void LLSQMultiDriver_t::save_llsq_state(const std::string &path) const
  {
    //    std::cout << __func__ << ": entering" << std::endl; 
    check_exit_lat();
    std::stringstream ss; 
    ss << path << "state_database.rad"; 
    ADATIO::BinaryFileWriter bin(ss.str());
    write(bin,*lattice_data); 
    bin.close();
  }


  void LLSQMultiDriver_t::save_ff_state(const std::string &path) const
  {
    //    std::cout << __func__ << ": entering" << std::endl; 
    check_exit_lat();
    std::stringstream ss; 
    ss << path << "ff_database.rad"; 
    ADATIO::BinaryFileWriter bin(ss.str());

    FormFacSolutions<std::complex<double> >
      ff(peek_FF(),peek_tags()); 

    write(bin,ff); 
    bin.close();
  }

  SEMBLE::SembleMatrix<std::complex<double> > 
    LLSQMultiDriver_t::get_rephase_FF(void) const
    {
      typedef std::complex<double> T; 
      SEMBLE::SembleMatrix<T> un = peek_FF(); 
      ENSEM::EnsemVectorComplex in; 
      in = get_ensem_row(0,un); 

      ENSEM::EnsemVectorReal real, imag;

      double const_real, const_imag; 

      real = ENSEM::real(in);
      imag = ENSEM::imag(in);

      std::vector<double> t;
      for(int i = 0; i < in.numElem(); ++i)
        t.push_back(double(i)); 

      // take the middle 70 %
      int thigh = int( double(t.size()) * .85); 
      int tlow = int( double(t.size()) *.15);  
      EnsemData ereal(t,real),eimag(t,imag); 
      ereal.hideDataAboveX(thigh - 0.1); 
      eimag.hideDataBelowX(tlow -0.1); 

      ADAT::Handle<FitFunction> freal(new ThreePointConstant), fimag(new ThreePointConstant);  
      JackFit fit_real(ereal,freal), fit_imag(eimag,fimag); 

      fit_real.runAvgFit(); 
      fit_imag.runAvgFit(); 
      const_real = fit_real.getAvgFitParValue(0);
      const_imag = fit_imag.getAvgFitParValue(0);

      // apparently this guy is zero, just dump the whole thing 
      if( (const_real < 2e-2) && (const_imag < 2e-2) ) 
        return un;

      // something went wrong, deal with it higher up 
      if ( isnan(const_real) && isnan(const_imag))
        return un; 

      double fit_phase = std::arg(std::complex<double>(const_real,const_imag)); 
      double apply_phase = -fit_phase; 

      // apply this phase to everyone!
      std::complex<double> phase(cos(apply_phase),sin(apply_phase));  
      un *= phase; 

      return un; 
    }


  bool LLSQMultiDriver_t::solve_fast(const std::string &soln_ID)
  {
    check_exit_lat();    
    rHandle<LLSQBaseSolver_t<std::complex<double> > > 
      my_solver = LLSQSolverFactoryEnv::callFactory(soln_ID);
    POW2_ASSERT(&*my_solver); 

    bool success(true); 

    generate_kinematic_factors();

    Kinv = my_solver->inv(K); 
    solver_log = my_solver->solution_log(); 
    init_Kinv = true; 

    FF_t = Kinv * lattice_data->data(); 
    init_FF = true; 

    /*
       std::cout << __func__ << std::endl;
       std::cout << "K = \n" << K.mean() << std::endl;
       std::cout << "Kinv = \n" << Kinv.mean() << std::endl;
       std::cout << "lat = \n" << SEMBLE::mean(lattice_data->data()) << std::endl;
       std::cout << "FF = \n" << FF_t.mean() << std::endl;
       */

    return success; 
  }


  bool LLSQMultiDriver_t::solve_slow(const std::string &soln_ID)
  {
    bool success(true); 
    check_exit(false,__func__); 
    return success;
  }


  void LLSQMultiDriver_t::generate_kinematic_factors(void)
  {
    check_exit_lat(); 
    std::vector<LatticeMultiDataTag> tags = lattice_data->tags(); 
    std::vector<LatticeMultiDataTag>::const_iterator it; 

    for(it = tags.begin(); it != tags.end(); ++it)
    {
      ffKinematicFactors_t<std::complex<double> >
        KK(FormFactorDecompositionFactoryEnv::callFactory(it->mat_elem_id)); 
      SEMBLE::SembleMatrix<std::complex<double> > work;
      work = KK.genFactors(makeMomInvariants(*it));

      /*  std::cout << __FILE__ << __func__ 
          << "\n" <<  it->file_id << std::endl;
          std::cout << SEMBLE::mean(work) << std::endl;
          */
      if(it == tags.begin())
        init_dim(K, work.getRow(
              remap_jmu_4_to_0(it->jmu)));
      else
        K.append_row(work.getRow(
              remap_jmu_4_to_0(it->jmu)));
    }

    init_K = true; 
  }

  void LLSQMultiDriver_t::init_false(void)
  { 
    init_lat = false;
    init_K = false;
    init_Kinv = false;
    init_FF = false;
  }


  void LLSQMultiDriver_t::check_exit(const bool b, const char *c) const
  {
    if(!!!b)
    {
      std::cerr << __func__ << ": error: called by " << c << ", exiting." << std::endl;
      exit(1); 
    }
  }







}


#undef DEBUG_AT_MAKE_MOM_INV_TAGS
