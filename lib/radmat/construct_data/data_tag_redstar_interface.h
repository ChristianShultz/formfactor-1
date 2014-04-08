#ifndef DATA_TAG_REDSTAR_INTERFACE_H
#define DATA_TAG_REDSTAR_INTERFACE_H 

#include "radmat/data_representation/data_representation.h"
#include "radmat/redstar_interface/redstar_interface.h"
#include "radmat/ff/ff.h"


namespace radmat
{

  struct TaggedEnsemRedstarNPtBlock
  {
    TaggedEnsemRedstarNPtBlock() {}

    TaggedEnsemRedstarNPtBlock( const EnsemRedstarNPtBlock &e, 
        const ThreePointDataTag &t)
      : coeff_lattice_xml(e) , data_tag(t) 
    { }

    double qsq_tag() const { return data_tag.get_qsq_label(); }
    std::string rot_qsq_tag(const bool mode) const
    {
      return  ( mode ) ? rot_qsq_tag( data_tag.origin_rep ) : rot_qsq_tag( data_tag.data_rep ); 
    } 
    std::string rot_qsq_tag( const DataRep3pt &data_rep ) const
    {
      std::stringstream ss; 
      ss << qsq_tag() << "_" 
        << LatticeRotationEnv::rotation_group_label(
            data_tag.left_mom,data_tag.right_mom); 
      ss << "__" << data_rep.l 
        << "," << data_rep.r; 
      return ss.str();  
    }

    EnsemRedstarNPtBlock coeff_lattice_xml; 
    ThreePointDataTag data_tag; 
  }; 


  std::vector<TaggedEnsemRedstarNPtBlock>
    tag_lattice_xml(const std::vector<ThreePointData> &, 
        const double mom_factor, 
        const double m_lefty, 
        const double m_righty, 
        const std::string &elem_id_base); 

}

#endif /* DATA_TAG_REDSTAR_INTERFACE_H */
