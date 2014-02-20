/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : lattice_mulit_data_tag.cc

* Purpose :

* Creation Date : 01-10-2013

* Last Modified : Thu 20 Feb 2014 01:53:50 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "lattice_multi_data_tag.h"
#include "radmat_database_interface.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/perThreadStorage.h"
#include "radmat/llsq/llsq_gen_system.h"
#include "semble/semble_semble.h"
#include "radmat_overlap_key_val_db.h"
#include "hadron/ensem_filenames.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"
#include "jackFitter/plot.h"
#include "ensem/ensem.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <math.h>
#include <omp.h>


namespace radmat
{
  LatticeMultiDataTag::LatticeMultiDataTag(void)
  {
    qsq_label = 1000.;
    E_f.resize(1); 
    E_f = SEMBLE::toScalar(double(0.));
    E_i = E_f; 
  }

  LatticeMultiDataTag& LatticeMultiDataTag::operator=(const LatticeMultiDataTag &o)
  {
    if(this != &o)
    {
      qsq_label = o.qsq_label; 
      jmu = o.jmu;
      hf = o.hf;
      hi = o.hi; 
      mat_elem_id = o.mat_elem_id;
      p_f = o.p_f;
      p_i = o.p_i;
      E_f = o.E_f;
      E_i = o.E_i;
      mom_fac = o.mom_fac;
      file_id = o.file_id;  
    }
    return *this;
  }

  ENSEM::EnsemReal LatticeMultiDataTag::Q2(void) const
  {
    double pp(0); 
    pp = mom_fac*mom_fac*((p_f[0] - p_i[0])*(p_f[0] - p_i[0])
        + (p_f[1] - p_i[1])*(p_f[1] - p_i[1])
        + (p_f[2] - p_i[2])*(p_f[2] - p_i[2]));


    return ( - (E_f-E_i)*(E_f-E_i) + SEMBLE::toScalar(pp));
  }

  void LatticeMultiDataTag::print_me(void) const
  {
    std::cout << file_id << " " << jmu << " " << mat_elem_id << std::endl;  
  }

  std::string LatticeMultiDataTag::splash_tag(void) const
  {
    std::stringstream ss; 
    ss << "f = " << file_id << " mu = " << jmu << " id = " << mat_elem_id; 
    ss << " Q2 = " << SEMBLE::toScalar(ENSEM::mean(Q2())) << " "; 
    ss << mom_string() << " " << E_string() <<  std::endl;
    return ss.str();
  }


  std::string LatticeMultiDataTag::mom_string(void) const
  {

    if(p_f.size() != 3 || p_i.size() != 3)
    {
      std::cerr <<__func__ << ": error, momenta don't have correct size" << std::endl;
      exit(1);
    }

    std::stringstream ss;
    ss << "pf = " << p_f[0] << "," << p_f[1] << ","
      << p_f[2] << "  pi = "  << p_i[0] << "," 
      << p_i[1] << "," << p_i[2] ;
    return ss.str();
  }


  std::string LatticeMultiDataTag::E_string(void) const
  {
    std::stringstream ss;
    ss << "E_f = " << std::setw(3) << SEMBLE::toScalar(ENSEM::mean(E_f)) 
      << " E_i = " << std::setw(3) << SEMBLE::toScalar(ENSEM::mean(E_i));
    return ss.str(); 
  }



  void write(ADATIO::BinaryWriter &bin, const LatticeMultiDataTag &t) 
  {

    ADATIO::write(bin,t.qsq_label);
    ADATIO::write(bin,t.jmu); 
    ADATIO::write(bin,t.hf); 
    ADATIO::write(bin,t.hi); 
    ADATIO::writeDesc(bin,t.mat_elem_id);
    ADATIO::write(bin,t.p_f);
    ADATIO::write(bin,t.p_i);
    ENSEM::write(bin,t.E_f);
    ENSEM::write(bin,t.E_i);
    ADATIO::write(bin,t.mom_fac); 
  }


  void read(ADATIO::BinaryReader &bin, LatticeMultiDataTag &t) 
  {
    ADATIO::read(bin,t.qsq_label);
    ADATIO::read(bin,t.jmu); 
    ADATIO::read(bin,t.hf); 
    ADATIO::read(bin,t.hi); 
    ADATIO::readDesc(bin,t.mat_elem_id);
    ADATIO::read(bin,t.p_f);
    ADATIO::read(bin,t.p_i);
    ENSEM::read(bin,t.E_f);
    ENSEM::read(bin,t.E_i);
    ADATIO::read(bin,t.mom_fac); 
  }


} // radmat
