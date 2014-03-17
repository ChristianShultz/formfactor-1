#ifndef LATTICE_MULIT_DATA_TAG_H
#define LATTICE_MULIT_DATA_TAG_H 


#include "radmat/ff_interface/formfactor_invariants.h"
#include "radmat/utils/handle.h"
#include "ensem/ensem.h"
#include "io/adat_xmlio.h"
#include <complex>
#include <string>
#include <iostream>
#include <iomanip>



namespace radmat
{

  // this is a tag for a row in the linear system 
  struct LatticeMultiDataTag
  {
    //! Constructor
    LatticeMultiDataTag(void);

    LatticeMultiDataTag& operator=(const LatticeMultiDataTag &o); 

    //! the lorentz index
    int mu(void) const {return jmu;}

    //! ensemble qsquared
    ENSEM::EnsemReal Q2(void) const;

    //! splash
    void print_me(void) const;

    //! got sick of typing this
    std::string mom_string(void) const;

    //! energies
    std::string E_string(void) const;

    //! the value of Q2 we use when sorting and labeling
    void set_qsq_label(const double &q2) {qsq_label = q2;}

    //! the value of Q2 we use when sorting and labeling 
    double get_qsq_label(void) const {return qsq_label;}

    //! splash it to the screen
    std::string unique_tag(void) const;

    // tags

    std::string file_id; // some unique string telling us what this is

    // for llsq system
    double qsq_label;        
    int jmu,hf,hi;                 
    std::string mat_elem_id; 
    ADATXML::Array<int> p_f;
    ADATXML::Array<int> p_i;
    ADATXML::Array<int> q; 
    ENSEM::EnsemReal E_f;
    ENSEM::EnsemReal E_i;
    double mom_fac; 

    bool have_reps;  
    std::string lefty,gamma,righty; 

  };


  // binary read write functions for serialization
  void write(ADATIO::BinaryWriter &bin, const LatticeMultiDataTag &t); 

  void read(ADATIO::BinaryReader &bin, LatticeMultiDataTag &t); 


} // radmat




#endif /* LATTICE_MULIT_DATA_TAG_H */
