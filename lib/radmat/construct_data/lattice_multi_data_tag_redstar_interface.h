#ifndef LATTICE_MULTI_DATA_TAG_REDSTAR_INTERFACE_H
#define LATTICE_MULTI_DATA_TAG_REDSTAR_INTERFACE_H 

#include "radmat/redstar_interface/redstar_interface.h"
#include "lattice_multi_data_tag.h"
#include <vector>

namespace radmat
{

  // data storage structure
  struct TaggedEnsemRedstarNPtBlock
  {
    TaggedEnsemRedstarNPtBlock() {}

    TaggedEnsemRedstarNPtBlock(const EnsemRedstarNPtBlock &e, 
        const LatticeMultiDataTag &t)
      : coeff_lattice_xml(e) , continuum_tag(t) 
    {  }

    double qsq_tag(void) const {return continuum_tag.get_qsq_label();}

    EnsemRedstarNPtBlock coeff_lattice_xml; 
    LatticeMultiDataTag continuum_tag; 
  }; 

  // translate from redstar lattice language 
  // to continuum lorentz language, no more 
  // abstract base class pointers too!! 
  //
  // const pointer to const data to avoid any 
  // copies of the npoints more than once
  std::vector<TaggedEnsemRedstarNPtBlock>
    tag_lattice_xml( const AbsRedstarMergeNPtData_t * const,
        const double mom_factor, 
        const double m_snk,
        const double m_src,
        const std::string &elem_id); 

  LatticeMultiDataTag 
    retrieve_tag(const std::vector<AbsRedstarInput_t*> &elem,
        const double mom_factor, 
        const double m_snk,
        const double m_src,
        const std::string &elem_id); 

} // radmat



#endif /* LATTICE_MULTI_DATA_TAG_REDSTAR_INTERFACE_H */
