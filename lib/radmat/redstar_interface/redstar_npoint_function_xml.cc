/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_npoint_function_xml.cc

 * Purpose :

 * Creation Date : 12-11-2013

 * Last Modified : Mon 17 Mar 2014 11:43:25 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_npoint_function_xml.h"
#include "redstar_abstract_block_base.h"
#include "redstar_abstract_xml_factory.h"
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
      obj.param = TheRedstarAbstractXMLFactoryEnv::callFactory(obj.object_name); 
      std::cout << "constructed a " << obj.object_name << std::endl; 
      obj.param->read(ptop,std::string("param"));
    }
    catch(std::exception &e) 
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": error, e.what() = " << e.what() << std::endl; 
      throw e; 
    }
    catch(std::string &s)
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": error, " << s << std::endl; 
    }
    catch(...)
    {
      std::cout << __PRETTY_FUNCTION__ 
        << ": some non standard error" << std::endl; 
      throw std::string("in") + std::string(__PRETTY_FUNCTION__); 
    }
  }



  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      NPointXML &npt)
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"version",npt.version,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"ensemble",npt.ensemble,__PRETTY_FUNCTION__); 

    switch ( npt.version )
    {
      case 0:
        read(ptop,"NPoint",npt.npoint); 
        break; 

      default:
        std::cout << "unknown version " << npt.version << std::endl; 
        exit(1); 
    }
    // auto fill size info 
    npt.N = npt.npoint.size(); 

   //     std::cout << __PRETTY_FUNCTION__ << "\n    : read \n" << std::endl; 
   //     for(int i =0; i < npt.N; ++i)
   //       std::cout << "\nnpt[" << i << "]: " << npt.npoint[i].object_name 
   //         << " -> " << npt.npoint[i].param->type() 
   //         << "\n      " << npt.npoint[i].param->write() 
   //         << std::endl;

  }

} // radmat

