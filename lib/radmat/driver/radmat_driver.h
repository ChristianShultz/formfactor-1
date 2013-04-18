#ifndef RADMAT_DRIVER_H
#define RADMAT_DRIVER_H


#include <string>
#include "adat/handle.h"
#include "radmat/load_data/build_correlators.h"
#include "radmat/llsq/llsq_multi_data.h"
#include "radmat_driver_props.h"
#include "radmat_single_q2_driver.h"

namespace radmat
{

  struct RadmatDriver
  {
    RadmatDriver(void) {}
    
    void run_program(const std::string &inifile);

    private:

    bool read_ini, built_correlators, solved_llsq, fit_formfacs;
    bool chisq_analysis; 

    void init_false(void); 
    void check_exit_ini(void) const {check_exit(read_ini,__func__);}
    void check_exit_corrs(void) const {check_exit(built_correlators,__func__);}
    void check_exit_llsq(void) const {check_exit(solved_llsq,__func__);}
    void check_exit_fit(void) const {check_exit(fit_formfacs,__func__);}
    void check_exit_chisq(void) const {check_exit(chisq_analysis,__func__);}
    void check_exit(const bool &, const char *) const;

    void read_xmlini(const std::string &ini); 
    void build_correlators(void);
    void solve_llsq(void); 
    void fit_ffs(void); 
    void do_chisq_analysis(void);  
    void make_FF_of_Q2_plots(void);
    void print_Q2_list(void);  

    RDriverProps_t m_ini;  
    BuildCorrelators m_correlators; 
    std::vector<bool> good_qs;
    std::vector<ADAT::Handle<LLSQLatticeMultiData> > multi_lattice_data; 
    std::vector<RadmatSingleQ2Driver> linear_systems_of_Q2; 
  };


}


#endif /* RADMAT_DRIVER_H */
