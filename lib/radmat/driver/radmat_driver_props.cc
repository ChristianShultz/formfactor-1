/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : radmat_driver_props.cc

* Purpose :

* Creation Date : 29-11-2012

* Last Modified : Wed Apr 24 10:14:16 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/driver/radmat_driver_props.h"

#include <sstream>


namespace
{

  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
          << ", path was empty, exiting" << std::endl;
        exit(1);
      }
    }

} // namespace anonomyous 




namespace radmat
{

  std::string toString (const RDriverProps_t &prop)
  {
    std::stringstream ss;
      ss << "threePointComparatorProps = " << prop.threePointComparatorProps;
      ss << " \nthreePointIni = " << prop.threePointIni;
      ss << "\npoleMass^2 = " <<  prop.poleMass; 
      return ss.str(); 
  }

  std::ostream& operator<<(std::ostream& o, const RDriverProps_t &p)
  {
    o << toString(p);
    return o; 
  }
  void read(ADATXML::XMLReader &xml, const std::string &path, RDriverProps_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"threePointComparatorProps",prop.threePointComparatorProps,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"threePointIni",prop.threePointIni,__PRETTY_FUNCTION__);
    
    // NB: we square the mass here!
    double pole_mass; 
    doXMLRead(ptop,"poleMass",pole_mass,__PRETTY_FUNCTION__); 
    prop.poleMass = pole_mass*pole_mass; 
  } 


} // namespace radmat 
