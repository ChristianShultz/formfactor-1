/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_merge_function_xml.cc

 * Purpose :

 * Creation Date : 12-11-2013

 * Last Modified : Thu 14 Nov 2013 08:24:17 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_merge_function_xml.h"
#include "redstar_abstract_merge_npoint_factory.h"


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

  }


  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      AbstractMergeNamedObject &p)
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"object_name",p.object_name,__PRETTY_FUNCTION__); 

    try
    {
      p.param = TheRedstarAbstractMergeNPtFactoryEnv::callFactory(p.object_name); 
      std::cout << "constructed a " << p.object_name << std::endl;
      p.param->read(ptop,std::string("param"));
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


} // radmat 



