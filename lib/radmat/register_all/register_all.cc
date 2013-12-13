/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 10-12-2013

 * Last Modified : Thu 12 Dec 2013 03:49:15 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat/ff/formfactor_factory.h"
#include "radmat/ff/lorentzff_canonical_rotations.h"
#include "radmat/llsq/llsq_solvers.h"
#include "radmat/redstar_interface/redstar_abstract_merge_npoint_factory.h"
#include "radmat/redstar_interface/redstar_abstract_xml_factory.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_invert_subduction.h"
#include <exception>

namespace radmat
{
  namespace AllFactoryEnv
  {

    namespace 
    {
      bool local_registration = false; 
    }


    bool registerAll(void)
    {
      bool success = true; 
      if(!!! local_registration ) 
      {

        try
        {
          success &= radmat::InvertSubductionEnv::registerAll(); 
          success &= radmat::FormFactorDecompositionFactoryEnv::registerAll(); 
          success &= radmat::LLSQSolverFactoryEnv::registerAll(); 
          success &= radmat::TheRedstarAbstractMergeNPtFactoryEnv::registerAll(); 
          success &= radmat::CanonicalRotationEnv::registerAll(); 
          success &= radmat::CanonicalLatticeRotationEnv::registerAll(); 
          success &= radmat::LatticeRotationEnv::registerAll(); 
          success &= radmat::TheRedstarAbstractXMLFactoryEnv::registerAll(); 
        }
        catch(std::exception &e)
        {
          std::cout << "std::exception: " << e.what() << std::endl; 
          std::cout << "dying in " << __PRETTY_FUNCTION__ << ":" 
            <<__FILE__ << ":" << __LINE__ << std::endl;
          exit(1); 
        }
        catch(std::string &s)
        {
          std::cout << "string exception: " << s << std::endl; 
          std::cout << "dying in " << __PRETTY_FUNCTION__ << ":" 
            <<__FILE__ << ":" << __LINE__ << std::endl;
          exit(1); 
        }
        catch(...)
        {
          std::cout << "some non standard non string exception" << std::endl; 
          std::cout << "dying in " << __PRETTY_FUNCTION__ << ":" 
            <<__FILE__ << ":" << __LINE__ << std::endl;
          exit(1); 
        }

        local_registration = true;  
      }

      return success; 
    }


  } 


}

