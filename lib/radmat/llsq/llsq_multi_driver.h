#ifndef LLSQ_MULTI_DRIVER_H
#define LLSQ_MULTI_DRIVER_H


#include "llsq_driver.h"
#include "llsq_solvers.h"
#include "llsq_gen_system.h"
#include "radmat/load_data/build_correlators.h"
#include <complex>
#include <string>

// roughly this is an optimization class, if its possible to just invert the matrix once 
// then we should be intelligent enough to do so, Jo's fancy explicit chisq extremization 
// method is an example where we actually have to do the inversion on every t_ins
namespace radmat
{



  struct LLSQMultiDriver_t
  {
    typedef std::complex<double> T;

    LLSQMultiDriver_t(void);

    bool load_data(const ADAT::Handle<LLSQLatticeMultiData> &dat); 
    bool run_zero_filter(void); 
    void solve_llsq(const std::string &soln_ID);
    void chisq_analysis(const SEMBLE::SembleVector<double> &ff, 
        const std::string &path,
        const int tlow,
        const int thigh,
        double tol=1e-6); 
    void dump_llsq(const std::string &path); 

    SEMBLE::SembleMatrix<T> peek_K(void) const {check_exit_K();return K;}
    SEMBLE::SembleMatrix<T> peek_Kinv(void) const {check_exit_Kinv();return Kinv;}
    SEMBLE::SembleMatrix<T> peek_FF(void) const {check_exit_FF(); return FF_t;}
    SEMBLE::SembleMatrix<T> peek_data(void) const {check_exit_lat(); return lattice_data->data();}
    std::vector<LatticeMultiDataTag> peek_tags(void) const {check_exit_lat(); return lattice_data->tags();}
    ENSEM::EnsemReal Q2(void) const {check_exit_lat(); return lattice_data->tags().begin()->Q2();}
    double qsq_label(void) const {check_exit_lat(); return lattice_data->tags().begin()->get_qsq_label();}

    private:

    bool solve_slow(const std::string &soln_ID); 
    bool solve_fast(const std::string &soln_ID); 

    void generate_kinematic_factors(void); 

    void init_false(void);
    void check_exit_lat(void) const {check_exit(init_lat,__func__);} 
    void check_exit_K(void) const {check_exit(init_K,__func__);}
    void check_exit_Kinv(void) const {check_exit(init_Kinv,__func__);}
    void check_exit_FF(void) const {check_exit(init_FF,__func__);}
    void check_exit(const bool, const char *) const ;     

    bool init_lat,init_K,init_Kinv,init_FF; 
    ADAT::Handle<LLSQLatticeMultiData> lattice_data;
    SEMBLE::SembleMatrix<T> K,Kinv,FF_t;
    std::vector<LatticeMultiDataTag> zeroed_elems; 
  };



}






#endif /* LLSQ_MULTI_DRIVER_H */
