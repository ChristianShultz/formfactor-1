/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : radmat_driver_props.cc

* Purpose :

* Creation Date : 29-11-2012

* Last Modified : Mon 30 Sep 2013 11:06:09 AM EDT

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
      ss << "version = " << prop.version; 
      ss << "\nthreePointComparatorProps = " << prop.threePointComparatorProps;
      ss << "\nthreePointIni = " << prop.threePointIni;
      ss << "\nmaxThread = " << prop.maxThread; 
      ss << "\npoleMass^2 = " <<  prop.poleMass; 
      return ss.str(); 
  }

  std::ostream& operator<<(std::ostream& o, const RDriverProps_t &p)
  {
    o << toString(p);
    return o; 
  }

  
  namespace
  {
    void check_version(const int v)
    {
      int my_version  = 1; 
      if ( my_version != 1 ) 
      {
        std::cerr << "version " << v << "is not supported, must be at " 
          << my_version << " try again later" << std::endl; 
        exit(1); 
      }
    }
  }


  void read(ADATXML::XMLReader &xml, const std::string &path, RDriverProps_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"version",prop.version,__PRETTY_FUNCTION__); 

    check_version(prop.version); 

    doXMLRead(ptop,"threePointComparatorProps",prop.threePointComparatorProps,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"threePointIni",prop.threePointIni,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"maxThread",prop.maxThread,__PRETTY_FUNCTION__);
    
    // NB: we square the mass here!
    double pole_mass; 
    doXMLRead(ptop,"poleMass",pole_mass,__PRETTY_FUNCTION__); 
    prop.poleMass = pole_mass*pole_mass; 
  } 


} // namespace radmat 
