#ifndef RADMAT_SINGLE_Q2_DRIVER_H
#define RADMAT_SINGLE_Q2_DRIVER_H 

#include "radmat/utils/handle.h"
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

    bool load_llsq(const rHandle<LLSQLatticeMultiData> &lattice_data,
        const double pole_mass_squared, 
        const double tolerance,
        const bool mix_irreps);

    // single q2 backdoor
    bool load_llsq(const rHandle<LLSQLatticeMultiData> &lattice_data,
        const double tolerance);

    void solve_llsq(const std::string &soln_ID);

    void save_llsq_state(void) const;
    void save_ff_of_t(void) const; 

    void fit_data(const ThreePointComparatorProps_t &fitProps, 
        const int tsrc,
        const int tsnk);

    // single q2 fit individual ffs
    void fit_and_dump_single_ffs(const ThreePointComparatorProps_t &props,
        const SEMBLE::SembleMatrix<std::complex<double> > &F_t,
        const ENSEM::EnsemReal &Q2,
        const int tsrc, 
        const int tsnk,
        const int ff,
        const int ff_max) const;

    void fit_and_dump_single_ffs(const ThreePointComparatorProps_t &props,
        const int tsrc, 
        const int tsnk,
        const int ff,
        const int ff_max) const
    {
      // do we have a solved llsq system
      check_exit_linear_system();
      check_exit_solved_llsq(); 
      fit_and_dump_single_ffs(props,
          linear_system.peek_FF(),
          linear_system.Q2(),
          tsrc,
          tsnk,
          ff,
          ff_max);
    }

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
    
    void append_rotation_group_label(const std::string &s) {rot_id = s;}

    private:
    std::string rotation_group_label(const bool mix_irreps) const; 
    std::string rot_id; 
    std::string base_path(void) const; 
    void init_false(); 
    void check_exit(const bool, const char *) const;    

    bool init_linear_system, init_solved_llsq, init_fits; 
    LLSQMultiDriver_t linear_system;
    TinsFitter fit_across_time;
  };


}












#endif /* RADMAT_SINGLE_Q2_DRIVER_H */
