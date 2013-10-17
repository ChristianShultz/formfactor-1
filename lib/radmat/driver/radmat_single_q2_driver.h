#ifndef RADMAT_SINGLE_Q2_DRIVER_H
#define RADMAT_SINGLE_Q2_DRIVER_H 

#include "radmat/load_data/build_correlators.h"
#include "radmat/llsq/llsq_multi_driver.h"
#include "radmat/fitting/fit_tins.h"
#include "ensem/ensem.h"

#include <string>
#include <utility>

namespace radmat
{


  struct RadmatSingleQ2Driver
  { 
    RadmatSingleQ2Driver(void);
    RadmatSingleQ2Driver(const RadmatSingleQ2Driver &o);
    RadmatSingleQ2Driver& operator=(const RadmatSingleQ2Driver &o);

    bool load_llsq(const ADAT::Handle<LLSQLatticeMultiData> &lattice_data,
        const double pole_mass_squared);

    // single q2 backdoor
    bool load_llsq(const ADAT::Handle<LLSQLatticeMultiData> &lattice_data);

    void solve_llsq(const std::string &soln_ID);
    void fit_data(const ThreePointComparatorProps_t &fitProps, 
        const int tsrc,
        const int tsnk);
    void chisq_analysis(const int tlow, const int thigh);

    ENSEM::EnsemReal Q2(void) const; 
    std::pair<ENSEM::EnsemReal, SEMBLE::SembleVector<double> > fetchFF(void) const;

    std::string tags_at_this_Q2(void) const; 

    void dump_fits(void);
    void dump_llsq(void);

    void check_exit_linear_system(void) const {check_exit(init_linear_system,__func__);}
    void check_exit_solved_llsq(void) const {check_exit(init_solved_llsq,__func__);}
    void check_exit_fits(void) const {check_exit(init_fits,__func__);}

    bool check_linear_system(void) const {return init_linear_system;}
    bool check_solved_llsq(void) const {return init_solved_llsq;}
    bool check_fits(void) const {return init_fits;}


    private:
    std::string base_path(void) const; 
    void init_false(); 
    void check_exit(const bool, const char *) const;    

    bool init_linear_system, init_solved_llsq, init_fits; 
    LLSQMultiDriver_t linear_system;
    TinsFitter fit_across_time;
  };


}












#endif /* RADMAT_SINGLE_Q2_DRIVER_H */
