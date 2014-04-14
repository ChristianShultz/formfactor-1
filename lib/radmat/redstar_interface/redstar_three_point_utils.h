#ifndef REDSTAR_THREE_POINT_UTILS_H
#define REDSTAR_THREE_POINT_UTILS_H 

#include "redstar_particle_handler_utils.h"
#include "radmat/data_representation/data_representation.h"

namespace radmat
{
  // the fundamental xml work unit 
  typedef ListObjExpr_t<ENSEM::Complex,
          Hadron::KeyHadronNPartNPtCorr_t> EnsemRedstarNPtBlock; 


  // the origin_rep is some lorentz info 
  // the data_rep is either lorentz or cubic
  // the row structure corresponds to data_rep 
  // the data is a list of coefficient times xml 
  struct ThreePointData
  {
    DataRep3pt origin_rep; 
    DataRep3pt data_rep; 
    int left_row, right_row, gamma_row; 
    EnsemRedstarNPtBlock data;
  }; 


  std::vector<ThreePointData>
    merge_blocks(const std::vector<BlockData> &lefty,
        const std::vector<BlockData> &gamma,
        const std::vector<BlockData> &righty,
        const std::string &ensemble); 


  std::vector<ThreePointData>
    merge_blocks(const std::vector<BlockData> &lefty,
        const std::vector<VectorCurrentImprovedBlockData> &gamma,
        const std::vector<BlockData> &righty,
        const std::string &ensemble); 

}

#endif /* REDSTAR_THREE_POINT_UTILS_H */
