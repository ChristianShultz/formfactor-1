/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_npoint_function.cc

 * Purpose :

 * Creation Date : 11-11-2013

 * Last Modified : Mon 11 Nov 2013 09:21:41 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_npoint_function.h"
#include "redstar_abstract_block_base.h"
#include "redstar_abstract_block_factory.h"
#include <exception>
#include "io/adat_xmlio.h"

namespace radmat
{
  namespace
  {
    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, 
          const std::string &path, 
          T &place, const char * f)
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

  } // anonomyous 



  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      AbstractBlockNamedObject &obj)
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"object_name",obj.object_name,__PRETTY_FUNCTION__); 

    try
    {
      obj.param = TheRedstarAbstractBlockFactoryEnv::callFactory(obj.object_name); 
      obj.param->read(ptop,std::string("param"));
    }
    catch(std::exception &e) 
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": error, e.what() = " << e.what() << std::endl; 
      throw e; 
    }
    catch(...)
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": some non standard error" << std::endl; 
      throw std::string("penguin"); 
    }
  }



  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      NPointXML &npt)
  {
    radmat::TheRedstarAbstractBlockFactoryEnv::registerAll(); 
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"version",npt.version,__PRETTY_FUNCTION__); 

    switch ( npt.version )
    {
      case 1:
        ADATXML::read(xml,"NPoint",npt.npoint); 
        break; 

      default:
        std::cout << "unknown version " << npt.version << std::endl; 
        exit(1); 
    }

  }

} // radmat


