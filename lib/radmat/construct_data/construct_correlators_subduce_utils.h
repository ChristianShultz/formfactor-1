#ifndef CONSTRUCT_CORRELATORS_SUBDUCE_UTILS_H
#define CONSTRUCT_CORRELATORS_SUBDUCE_UTILS_H 


#include "construct_correlators_xml.h"
#include "lattice_multi_data_tag_redstar_interface.h"
#include "lattice_multi_data_object.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "radmat_database_interface.h"
#include "radmat_overlap_key_val_db.h"
#include "radmat/utils/mink_qsq.h"
#include "ensem/ensem.h"
#include "adat/map_obj.h"
#include <vector>
#include <string>
#include <map> 
#include <utility>

namespace radmat
{

  std::vector<TaggedEnsemRedstarNPtBlock> 
    retag_subduced_lattice_xml( 
        const std::vector<TaggedEnsemRedstarNPtBlock> &cont_variant); 


}



#endif /* CONSTRUCT_CORRELATORS_SUBDUCE_UTILS_H */
