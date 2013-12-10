/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_abstract_merge_npoint_factory.cc

* Purpose :

* Creation Date : 12-11-2013

* Last Modified : Tue 10 Dec 2013 11:21:49 AM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_abstract_merge_npoint_factory.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"

#include "redstar_merge_vector_current_npoint.h"


namespace FacEnv = radmat::TheRedstarAbstractMergeNPtFactoryEnv;
typedef radmat::TheRedstarAbstractMergeNPtFactory Factory; 

namespace radmat
{

  namespace TheRedstarAbstractMergeNPtFactoryEnv
  {
    volatile bool local_registration = false; 

    namespace
    {
      template<typename T> 
        AbsRedstarMergeNPt* 
        callback(void)
        {
          T *foo = new T; 
          return foo; 
        }
    }



    template<typename T> 
      bool do_reg(void)
      {
        return Factory::Instance().registerObject(Stringify<T>(),callback<T>); 
      }


    bool registerAll(void)
    {
      bool success = true; 
      if ( !!! local_registration ) 
      {
        success &= do_reg<RedstarMergeVectorCurrentThreePoint>();
      }

      if ( !!! success ) 
        throw std::string("reg error in TheRedstarAbstractMergeNPtFactoryEnv"); 

      return success; 
    }

    rHandle<AbsRedstarMergeNPt>
      callFactory(const std::string &id)
      {
        POW2_ASSERT( FacEnv::registerAll() ); 
        AbsRedstarMergeNPt *foo; 
        try
        {
          foo = Factory::Instance().createObject(id);  
        }
        catch(...)
        {
          POW2_ASSERT(false); 
        }

        POW2_ASSERT(foo); 

        return rHandle<AbsRedstarMergeNPt>(foo); 
      }

  } // TheRedstarAbstractMergeNPtFactoryEnv

} // radmat


