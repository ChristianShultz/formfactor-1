/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_abstract_block_factory.cc

 * Purpose :

 * Creation Date : 11-11-2013

 * Last Modified : Mon 11 Nov 2013 08:58:21 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_abstract_block_factory.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"

#include "redstar_single_particle_meson_block.h"
#include "redstar_unimproved_vector_current.h"

namespace FacEnv = radmat::TheRedstarAbstractBlockFactoryEnv;
typedef radmat::TheRedstarAbstractBlockFactory Factory; 

namespace radmat
{

  namespace TheRedstarAbstractBlockFactoryEnv
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



    template<typename T1,typename T2> 
      bool do_reg(void)
      {
        std::string s(Stringify<T1>()); 
        return Factory::Instance().registerObject(s,callback<T2>); 
      }


    bool registerAll(void)
    {
      bool success = true; 
      if ( !!! local_registration ) 
      {
        success &= do_reg<RedstarSingleParticleMesonXML_string,
                RedstarSingleParticleMesonXML>();
        success &= do_reg<RedstarUnimprovedVectorCurrentXML_string,
                RedstarUnimprovedVectorCurrentXML>(); 
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


