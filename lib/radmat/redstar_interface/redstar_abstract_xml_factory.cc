/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_abstract_xml_factory.cc

* Purpose :

* Creation Date : 12-11-2013

* Last Modified : Thu 14 Nov 2013 09:24:11 AM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_abstract_xml_factory.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"

#include "redstar_single_particle_meson_block.h"
#include "redstar_unimproved_vector_current.h"

namespace FacEnv = radmat::TheRedstarAbstractXMLFactoryEnv;
typedef radmat::TheRedstarAbstractXMLFactory Factory; 

namespace radmat
{

  namespace TheRedstarAbstractXMLFactoryEnv
  {
    volatile bool local_registration = false; 

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
      bool success = true; 
      if ( !!! local_registration ) 
      {
        success &= do_reg<RedstarSingleParticleMesonXML>();
        success &= do_reg<RedstarUnimprovedVectorCurrentXML>(); 
        local_registration = true; 
      }
      return success; 
    }

    ADAT::Handle<AbsRedstarXMLInterface_t>
      callFactory(const std::string &id)
      {
        POW2_ASSERT( FacEnv::registerAll() ); 
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

        return ADAT::Handle<AbsRedstarXMLInterface_t>(foo); 
      }

  } // TheRedstarAbstractBlockFactoryEnv

} // radmat


