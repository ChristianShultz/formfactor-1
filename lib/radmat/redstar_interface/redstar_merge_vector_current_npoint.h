#ifndef REDSTAR_MERGE_VECTOR_CURRENT_NPOINT_H
#define REDSTAR_MERGE_VECTOR_CURRENT_NPOINT_H 

#include "redstar_abstract_merge_npoint.h"

namespace radmat
{

  // what do we do with a vector current three point 
  struct RedstarMergeVectorCurrentThreePoint;
  REGISTER_STRINGIFY_TYPE(RedstarMergeVectorCurrentThreePoint); 

  struct RedstarMergeVectorCurrentThreePoint
    : public AbsRedstarMergeNPt
  {
    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarMergeVectorCurrentThreePoint>()); 
    }

    // multiply this thing out 
    virtual void do_work(void) ;

  };


} // radmat 


#endif /* REDSTAR_MERGE_VECTOR_CURRENT_NPOINT_H */
