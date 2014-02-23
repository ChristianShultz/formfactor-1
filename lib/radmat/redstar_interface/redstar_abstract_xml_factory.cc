/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_abstract_xml_factory.cc

* Purpose :

* Creation Date : 12-11-2013

* Last Modified : Sun 23 Feb 2014 10:38:05 AM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_abstract_xml_factory.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/printer.h"

#include "redstar_single_particle_meson_block.h"
#include "redstar_unimproved_vector_current.h"
#include "redstar_improved_vector_current.h"

#include <exception>

namespace FacEnv = radmat::TheRedstarAbstractXMLFactoryEnv;
typedef radmat::TheRedstarAbstractXMLFactory Factory; 

namespace radmat
{

  namespace TheRedstarAbstractXMLFactoryEnv
  {
    bool local_registration = false; 

    namespace
    {
      template<typename T> 
        AbsRedstarXMLInterface_t* 
        callback(void)
        {
          AbsRedstarXMLInterface_t *foo = new T; 
          return foo; 
        }
    }



    template<typename T> 
      bool do_reg(void)
      {
        bool foo = Factory::Instance().registerObject(Stringify<T>(),callback<T>); 
        if (!!!foo)
          std::cerr << __PRETTY_FUNCTION__ 
            << ":error reging " << Stringify<T>() << std::endl;

        return foo; 
      }


    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 
      if ( !!! local_registration ) 
      {
        success &= do_reg<RedstarSingleParticleMesonXML>();
        success &= do_reg<RedstarUnimprovedVectorCurrentXML>(); 
        success &= do_reg<RedstarImprovedVectorCurrentXML>(); 
        local_registration = true; 
      }

      if(!!! success)
        throw std::string("reg error in TheRedstarAbstractBlockFactoryEnv"); 

      return success; 
    }

    rHandle<AbsRedstarXMLInterface_t>
      callFactory(const std::string &id)
      {
        AbsRedstarXMLInterface_t *foo; 
        try
        {
          foo = Factory::Instance().createObject(id);  
        }
        catch(...)
        {
          POW2_ASSERT(false); 
        }

        POW2_ASSERT(foo); 

        return rHandle<AbsRedstarXMLInterface_t>(foo); 
      }

  } // TheRedstarAbstractBlockFactoryEnv

} // radmat


