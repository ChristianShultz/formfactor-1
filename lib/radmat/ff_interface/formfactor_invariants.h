#ifndef FORMFACTOR_INVARIANTS_H
#define FORMFACTOR_INVARIANTS_H 

#include <string>

namespace radmat
{

  struct FFRep_p
  {
    virtual ~FFRep_p() {}
    virtual std::string rep_id(void) const = 0; 
    virtual int rep_row(void) const = 0; 
  };

}; 


#endif /* FORMFACTOR_INVARIANTS_H */
