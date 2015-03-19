#ifndef RADMAT_DRIVER_H
#define RADMAT_DRIVER_H


#include <string>
#include "radmat/utils/handle.h"
#include "radmat/construct_data/construct_correlators.h"
#include "radmat/llsq/llsq_multi_data.h"
#include "radmat_driver_props.h"
#include "radmat_single_q2_driver.h"

namespace radmat
{

  struct RadmatDriver
  {
    typedef std::map<std::string, std::vector<Hadron::KeyHadronNPartNPtCorr_t> > xml_map_t; 

    RadmatDriver(void) {}

    //! do a three point analysis and extract form factors    
    void run_program(const std::string &inifile);

    //! run various incarnations of xml generation (check_op=true will only return xml of known ops)
    void xml_handler(const std::string &inifile, const std::string &mode, bool check_op=false); 

    //! just run up to the point of building xml 
    xml_map_t build_xml(void); 

    //! split up xml on p^2
    xml_map_t build_xml_split_p2(void);

    //! split up xml on p^2 then split into N sets 
    //    -- N is set in the .cc file 
    xml_map_t build_xml_split_p2_N(void);

    //! split up xml on octant 
    xml_map_t build_xml_split(void); 

    //! two point xml hack
    xml_map_t build_xml_twopoint(void);

    //! build 2 tsrc xml -- mod 128 
    xml_map_t build_xml_2tsrc(void); 

    //! build 4 tsrc xml -- mod 128
    xml_map_t build_xml_4tsrc(void); 


    private:

    bool read_ini, built_correlators, init_llsq, solved_llsq, fit_formfacs;
    bool chisq_analysis; 

    void init_false(void); 
    void check_exit_ini(void) const {check_exit(read_ini,__func__);}
    void check_exit_corrs(void) const {check_exit(built_correlators,__func__);}
    void check_exit_init_llsq(void) const {check_exit(init_llsq,__func__);}
    void check_exit_llsq(void) const {check_exit(solved_llsq,__func__);}
    void check_exit_fit(void) const {check_exit(fit_formfacs,__func__);}
    void check_exit_chisq(void) const {check_exit(chisq_analysis,__func__);}
    void check_exit(const bool &, const char *) const;

    // for running the full analysis 
    bool read_xmlini(const std::string &ini); 
    bool build_correlators(void);
    bool solve_llsq(void);
    bool fit_ffs(void); 
    bool do_chisq_analysis(void);  
    bool make_FF_of_Q2_plots(void);
    bool print_Q2_list(void);  


    RDriverProps_t m_ini;  
    ConstructCorrelators m_correlators; 
    std::vector<bool> good_qs;
    std::vector<rHandle<LLSQLatticeMultiData> > multi_lattice_data; 
    std::vector<RadmatSingleQ2Driver> linear_systems_of_Q2; 
  };


}


#endif /* RADMAT_DRIVER_H */
