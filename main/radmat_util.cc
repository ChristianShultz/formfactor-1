/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat.cc

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Tue 24 Dec 2013 11:07:28 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <exception>
#include <iostream>

#include "radmat/ff/lorentzff_canonical_rotations_checker.h"
#include "radmat/register_all/register_all.h"
#include "radmat/driver/radmat_driver.h"
#include "radmat/construct_data/lattice_multi_data_object.h"
#include "radmat/llsq/llsq_multi_data_serialize.h"
#include "radmat/llsq/llsq_solution.h"
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



  void read(ADATXML::XMLReader &xml, const std::string &path, SingleQ2Prop_t &p) 
  {
    ADATXML::XMLReader ptop(xml,path); 

    doXMLRead(ptop,"ff",p.ff,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"threePointComparatorProps",p.threePointComparatorProps,__PRETTY_FUNCTION__); 
  }


  void read(ADATXML::XMLReader &xml, const std::string &path, ArrSingleQ2Prop_t &p)
  {
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
      trim.append_row_semble(inout->get_row_semble(p.lat_elems[i]),inout->get_tag(p.lat_elems[i])); 

    *inout = trim;
  }

} // anonomyous 



// UTILITY FUNCTIONS
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

//    generate redstar xml from the inifile 
//
/////////////////////////////////////////////////////
void gen_xml(int argc , char *argv[] )
{
  if(argc != 4)
  {
    std::cerr << "error: usage: radmat_util: gen_xml <xmlinifile> <mode> " << std::endl;
    exit(1); 
  }

  std::istringstream val(argv[2]); 
  std::string ini; 
  val >> ini; 

  std::istringstream val2(argv[3]); 
  std::string mode; 
  val2 >> mode; 

  radmat::RadmatDriver d; 
  d.xml_handler(ini,mode); 
}


//    use a statedatabase to check the rotation 
//    properties of the data 
//
/////////////////////////////////////////////////////
void rot_llsq(int argc, char *argv[])
{
  if(argc != 5)
  {
    std::cerr << "usage: radmat_util: rot_llsq <xmlinifile> <Jl> <Jr> " << std::endl;
    exit(1); 
  }


  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[2]);
  val >> xmlini;

  int Jl; 
  std::istringstream val1(argv[3]); 
  val1 >> Jl; 

  int Jr; 
  std::istringstream val2(argv[4]); 
  val2 >> Jr; 

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


  // get the lattice elems as a function of insertion time from 
  // the database that was saved in the orig run 
  pull_elems(foo,arr_ini); 

  radmat::LatticeRotationRelationChecker bar;

  try
  {
  bar.check(foo,Jl,Jr); 
  }
  catch (std::string &s)
  {
    std::cout << __func__ << ": caught " << s << std::endl;
  }

}




//    use a statedatabase to resolve a llsq 
//
/////////////////////////////////////////////////////
void Q2_llsq(int argc, char *argv[])
{
  if(argc != 4)
  {
    std::cerr << "usage: radmat_util: Q2_llsq <xmlinifile> <ffmax> " << std::endl;
    exit(1); 
  }


  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[2]);
  val >> xmlini;

  int ffmax; 
  std::istringstream val1(argv[3]); 
  val1 >> ffmax; 

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
        arr_ini.ffs[elem].ff,
        ffmax);
  } // ff loop

}




//    go back and refit form factors 
//
/////////////////////////////////////////////////////
void refit_ffs(int argc, char *argv[])
{
  if(argc != 4)
  {
    std::cerr << "usage: radmat_util: refit_ffs <xmlinifile> <ffmax> " << std::endl;
    exit(1); 
  }


  // read the name of the ini file 
  std::string xmlini;
  std::istringstream val(argv[2]);
  val >> xmlini;

  int ffmax; 
  std::istringstream val1(argv[3]); 
  val1 >> ffmax; 

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
  catch( std::exception &e)
  {
    std::cout << "excep: " << e.what() << std::endl; 
    exit(1);
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the inifile") ; 
    exit(1); 
  }

  radmat::FormFacSolutions<std::complex<double> > FF_of_t;

  try
  {
    ADATIO::BinaryFileReader bread(arr_ini.dbfile); 
    radmat::read(bread,FF_of_t); 
  }
  catch ( ... ) 
  {
    SPLASH("An error occurred reading the state database");
    exit(1);
  }

  // driver
  radmat::RadmatSingleQ2Driver my_driver;

  // loop them 
  for (int elem = 0; elem < arr_ini.ffs.size(); ++elem)
  {
    std::cout << "\n\n** refiting ff_" << arr_ini.ffs[elem].ff << std::endl;
    std::cout << "*********************************" << std::endl;
    // fit out the insertion time dependence
    my_driver.fit_and_dump_single_ffs(
        arr_ini.ffs[elem].threePointComparatorProps,
        FF_of_t.FF_t,
        FF_of_t.Ingredients.begin()->Q2(),
        arr_ini.tsrc,
        arr_ini.tsnk,
        arr_ini.ffs[elem].ff,
        ffmax);
  } // ff loop

}


//
//
//    WORK HANDLER STUFF
//
//

// typedef 
typedef void (*fptr)(int argc , char *argv[]) ; 

// a map of operation names and function pointers
std::map<std::string , fptr> options; 

// init the map 
void init_options(void)
{
  options.insert(std::pair<std::string,fptr>("gen_xml",&gen_xml)); 
  options.insert(std::pair<std::string,fptr>("rot_llsq",&rot_llsq)); 
  options.insert(std::pair<std::string,fptr>("Q2_llsq",&Q2_llsq)); 
  options.insert(std::pair<std::string,fptr>("refit_ffs",&refit_ffs)); 
}

// pick appropriate function and pass on command line inputs 
void do_work(std::string &op, int argc,char *argv[])
{
  init_options(); 

  if(options.find(op) == options.end())
  {
    std::cerr << " unrecognized op " << op 
      << " more intelligent choices are " << std::endl; 
    std::map<std::string , fptr>::const_iterator it; 
    for(it = options.begin(); it != options.end(); ++it)
      std::cerr << it->first << std::endl; 
    exit(1); 
  }

  fptr foo = options[op];

  foo(argc,argv); 
}


// main program wrapper
int main(int argc, char *argv[])
{
  radmat::AllFactoryEnv::registerAll(); 

  // we will always have at least 2 , radmat_util operation_with_no_inputs
  if(argc < 2)
  {
    std::cerr << "usage: radmat_util : <operation> <op inputs ...> " << std::endl;
    exit(1); 
  }

  std::string op;
  std::istringstream opi(argv[1]); 
  opi >> op; 

  do_work(op,argc,argv); 

  return 0;
}
