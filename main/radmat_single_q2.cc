/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat.cc

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Tue 01 Oct 2013 04:48:30 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include "radmat/driver/radmat_driver.h"
#include "radmat/load_data/lattice_multi_data_object.h"
#include "radmat/llsq/llsq_multi_data_serialize.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include "adat/handle.h"
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

  struct SingleQ2Prop_t
  {
    ThreePointComparatorProps_t threePointComparatorProps;
    std::string dbfile; 
    std::string solnID; 
    ADATXML::Array<int> lat_elems; 
  };

  void read(ADATXML::XMLReader &xml, const std::string &pth, SingleQ2Prop_t &p) 
  {
    doXMLRead(xml,"threePointComparatorProps",p.threePointComparatorProps,__PRETTY_FUNCTION__); 
    doXMLRead(xml,"dbfile",p.dbfile,__PRETTY_FUNCTION__); 
    doXMLRead(xml,"solnID",p.solnID,__PRETTY_FUNCTION__); 
    doXMLRead(xml,"lat_elems",p.lat_elems,__PRETTY_FUNCTION__); 
  }


  void pull_elems(ADAT::Handle< radmat::LLSQLatticeMultiData > &inout , const SingleQ2Prop_t &p)
  {
    radmat::LLSQLatticeMultiData trim;
    for (int i = 0; i < p.lat_elems.size(); ++i)
    {
      std::cout << "pulling row " << p.lat_elems[i] << " from the original llsq \n" << std::endl; 
      std::cout << "    tag -> " << (inout->get_tag(p.lat_elems[i])).splash_tag() << std::endl; 
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


  std::string xmlini;
  std::istringstream val(argv[1]);
  val >> xmlini;

  ADAT::Handle< radmat::LLSQLatticeMultiData > foo( new radmat::LLSQLatticeMultiData() ); 
  POW2_ASSERT( &*foo ) ; 
  SingleQ2Prop_t ini; 

  try
  {
    ADATXML::XMLReader xml(xmlini); 
    read(xml,"/Props",ini); 
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the inifile") ; 
    exit(1); 
  }


  try
  {
    ADATIO::BinaryFileReader bread(ini.dbfile); 
    radmat::read(bread,*foo); 
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the state database");
    exit(1);
  }

  radmat::RadmatSingleQ2Driver my_driver;

  pull_elems(foo,ini); 

  POW2_ASSERT( my_driver.load_llsq(foo) ); 

  my_driver.solve_llsq(ini.solnID); 

  my_driver.fit_data(ini.threePointComparatorProps);

  my_driver.chisq_analysis(); 
  

  return 0;
}
