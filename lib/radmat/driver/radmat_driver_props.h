#ifndef RADMAT_DRIVER_PROPS_H
#define RADMAT_DRIVER_PROPS_H

#include <string>
#include "io/adat_xmlio.h"
#include "jackFitter/three_point_fit_forms.h"
#include "radmat/load_data/build_correlators.h"

namespace radmat
{
  struct RDriverProps_t
  {
    ThreePointComparatorProps_t threePointComparatorProps;
    ThreePointCorrIni_t threePointIni;    
  };


  std::string toString (const RDriverProps_t &);
  std::ostream& operator<<(std::ostream& , const RDriverProps_t &);
  void read(ADATXML::XMLReader &xml, const std::string &path, RDriverProps_t &); 

}// namespace radmat


#endif /* RADMAT_DRIVER_PROPS_H */
