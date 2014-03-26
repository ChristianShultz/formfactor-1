/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : data_tag_three_point.cc

* Purpose :

* Creation Date : 24-03-2014

* Last Modified : Wed 26 Mar 2014 12:53:01 PM EDT

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


  ENSEM::EnsemReal 
    ThreePointDataTag::Q2() const
    {
      ENSEM::EnsemReal Q0, Q1, Q2, Q3; 
      ENSEM::Real zero,q1,q2,q3; 
      zero = ENSEM::toDouble(0.); 
      Q0 = (left_E - right_E) * (left_E - right_E); 
      Q1 = zero*Q0;
      Q2 = zero*Q0;
      Q3 = zero*Q0;

      q1 = ENSEM::toDouble( double( left_mom[0] - right_mom[0] ) ); 
      q2 = ENSEM::toDouble( double( left_mom[1] - right_mom[1] ) ); 
      q3 = ENSEM::toDouble( double( left_mom[2] - right_mom[1] ) ); 

      Q1 = q1*q1; 
      Q2 = q2*q2; 
      Q3 = q3*q3; 

      return -Q0 + Q1 + Q2 + Q2;
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
