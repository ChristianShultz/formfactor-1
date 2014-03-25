/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_tag_three_point.cc

* Purpose :

* Creation Date : 24-03-2014

* Last Modified : Mon 24 Mar 2014 11:10:09 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "data_tag_three_point.h"
#include "semble/semble_meta.h"


namespace radmat
{

  ThreePointDataTag::ThreePointDataTag()
  {
    qsq_label = -10000.;
    left_E.resize(1); 
    left_E = SEMBLE::toScalar(double(0.)); 
    right_E = left_E;  
  }


  void write(ADATIO::BinaryWriter &bin, const ThreePointDataTag &d)
  {
    write(bin,d.origin_rep);
    write(bin,d.data_rep);
    ADATIO::writeDesc(bin,d.file_id);
    ADATIO::writeDesc(bin,d.mat_elem_id); 
    ADATIO::write(bin,d.left_row);
    ADATIO::write(bin,d.gamma_row);
    ADATIO::write(bin,d.right_row);
    ADATIO::write(bin,d.left_mom);
    ADATIO::write(bin,d.right_mom);
    ENSEM::write(bin,d.left_E);
    ENSEM::write(bin,d.right_E); 
    ADATIO::write(bin,d.qsq_label); 
    ADATIO::write(bin,d.mom_fac); 
  }


  void read(ADATIO::BinaryReader &bin, ThreePointDataTag &d)
  {
    read(bin,d.origin_rep);
    read(bin,d.data_rep);
    ADATIO::readDesc(bin,d.file_id);
    ADATIO::readDesc(bin,d.mat_elem_id); 
    ADATIO::read(bin,d.left_row);
    ADATIO::read(bin,d.gamma_row);
    ADATIO::read(bin,d.right_row);
    ADATIO::read(bin,d.left_mom);
    ADATIO::read(bin,d.right_mom);
    ENSEM::read(bin,d.left_E);
    ENSEM::read(bin,d.right_E); 
    ADATIO::read(bin,d.qsq_label); 
    ADATIO::read(bin,d.mom_fac); 
  }



} 
