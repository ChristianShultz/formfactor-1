/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat.cc

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Tue 03 Dec 2013 09:05:04 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include "radmat/driver/radmat_driver.h"
#include "radmat/construct_data/lattice_multi_data_object.h"
#include "radmat/llsq/llsq_multi_data_serialize.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include "radmat/utils/handle.h"
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "jackFitter/three_point_fit_forms.h"


namespace
{
  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
          << f << " trying to read path, " << path
          << ", path was empty, exiting" << std::endl;
        exit(1);
      }
    }

  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f, const T &val)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
          << f << " trying to read path, " << path
          << ", path was empty, reverting to default value " 
          << val << std::endl;
        place = val;
      }
    }

  struct SingleQ2Prop_t
  {
    int ff;                                                 // which form factor are we refitting
    ThreePointComparatorProps_t threePointComparatorProps;  // how are we fitting it
  };

  struct ArrSingleQ2Prop_t
  {
    ADATXML::Array<SingleQ2Prop_t>  ffs;                  // the list of ffs that we want to refit
    int tsrc;                                               // some duplicate info
    int tsnk;                                               
    std::string dbfile;                                     // where does it live
    std::string solnID;                                     // how are we inverting   
    double tolerance;                                       // tolerance
    ADATXML::Array<int> lat_elems;                          // which elements are we using in the refit
  };


  void echo(const std::string &s , const char * c)
  {
    std::cout << c << s << std::endl; 
  }

  void read(ADATXML::XMLReader &xml, const std::string &path, SingleQ2Prop_t &p) 
  {
    ADATXML::XMLReader ptop(xml,path); 

    echo(path,__PRETTY_FUNCTION__); 

    doXMLRead(ptop,"ff",p.ff,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"threePointComparatorProps",p.threePointComparatorProps,__PRETTY_FUNCTION__); 
  }


  void read(ADATXML::XMLReader &xml, const std::string &path, ArrSingleQ2Prop_t &p)
  {
    echo(path,__PRETTY_FUNCTION__); 
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"ffs",p.ffs,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"tsrc",p.tsrc,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"tsnk",p.tsnk,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"dbfile",p.dbfile,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"solnID",p.solnID,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"tolerance",p.tolerance,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"lat_elems",p.lat_elems,__PRETTY_FUNCTION__); 
  }


  void pull_elems(radmat::rHandle< radmat::LLSQLatticeMultiData > &inout , const ArrSingleQ2Prop_t &p)
  {
    radmat::LLSQLatticeMultiData trim;
    for (int i = 0; i < p.lat_elems.size(); ++i)
    {
#if 0 
      std::cout << "pulling row " << p.lat_elems[i] << " from the original llsq \n" << std::endl; 
      std::cout << "    tag -> " << (inout->get_tag(p.lat_elems[i])).splash_tag() << std::endl; 
#endif
      trim.append_row_semble(inout->get_row_semble(p.lat_elems[i]),inout->get_tag(p.lat_elems[i])); 
    }
    *inout = trim;
  }

} // namespace anonomyous 






int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cerr << "usage: radmat : <xmlinifile> " << std::endl;
    exit(1); 
  }


  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[1]);
  val >> xmlini;

  // read the xml array of ffs that we want to refit
  ArrSingleQ2Prop_t arr_ini; 

  try
  {
    std::cout << "reading " << xmlini << std::endl;
    ADATXML::XMLReader xml(xmlini); 
    read(xml,"/props",arr_ini); 
  }
  catch( std::string &str) 
  {
    std::cout << "Error: " << str << std::endl; 
    exit(1); 
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the inifile") ; 
    exit(1); 
  }

  radmat::rHandle< radmat::LLSQLatticeMultiData > foo( new radmat::LLSQLatticeMultiData() ); 
  POW2_ASSERT( &*foo ) ; 

  try
  {
    ADATIO::BinaryFileReader bread(arr_ini.dbfile); 
    radmat::read(bread,*foo); 
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the state database");
    exit(1);
  }

  // driver
  radmat::RadmatSingleQ2Driver my_driver;

  // get the lattice elems as a function of insertion time from 
  // the database that was saved in the orig run 
  pull_elems(foo,arr_ini); 

  // check that we can load the thing
  POW2_ASSERT( my_driver.load_llsq(foo,arr_ini.tolerance) ); 

  // solve the linear system 
  my_driver.solve_llsq(arr_ini.solnID); 

  // loop them 
  for (int elem = 0; elem < arr_ini.ffs.size(); ++elem)
  {
    std::cout << "\n\n** refiting ff_" << arr_ini.ffs[elem].ff << std::endl;
    std::cout << "*********************************" << std::endl;
    // fit out the insertion time dependence
    my_driver.fit_and_dump_single_ffs(
        arr_ini.ffs[elem].threePointComparatorProps,
        arr_ini.tsrc,
        arr_ini.tsnk,
        arr_ini.ffs[elem].ff);
  } // ff loop


  return 0;
}
