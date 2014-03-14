/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 10-12-2013

 * Last Modified : Fri 14 Mar 2014 11:58:33 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat/ff/lorentzff_canonical_rotations.h"
#include "radmat/ff/lorentzff_Wigner_D_matrix_factory.h"
#include "radmat/ff/lorentzff_formfactor_factory.h"
#include "radmat/ff_interface/formfactor_invariants.h"
#include "radmat/ff_interface/formfactor_factory.h"
#include "radmat/llsq/llsq_solvers.h"
#include "radmat/redstar_interface/redstar_abstract_merge_npoint_factory.h"
#include "radmat/redstar_interface/redstar_abstract_xml_factory.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_invert_subduction.h"
#include "radmat/utils/printer.h"
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
      printer_function<console_print>("Registering All Factories"); 
      bool success = true; 
      if(!!! local_registration ) 
      {

        try
        {
          // build subduction tables
          success &= radmat::InvertSubductionEnv::registerAll(); 
          // build lorentz helicity decompositions
          success &= radmat::LorentzffFormFactorDecompositionFactoryEnv::registerAll(); 
          // build  reps 
          success &= radmat::FormFactorInvariantsFactoryEnv::registerAll(); 
          // build radmat's internal decompositions
          success &= radmat::FormFactorDecompositionFactoryEnv::registerAll(); 
          // build llsq solution classes
          success &= radmat::LLSQSolverFactoryEnv::registerAll(); 
          // redstar xml 
          success &= radmat::TheRedstarAbstractMergeNPtFactoryEnv::registerAll(); 
          // redstar xml
          success &= radmat::TheRedstarAbstractXMLFactoryEnv::registerAll(); 
          // redstar canonical rotations 
          success &= radmat::CanonicalRotationEnv::registerAll(); 
          // redstar all lattice rotations
          success &= radmat::CanonicalLatticeRotationEnv::registerAll(); 
          // some more rotations
          success &= radmat::LatticeRotationEnv::registerAll(); 
          // d matricies 
          success &= radmat::WignerDMatrixEnv::registerAll(4); // up to J = 2 
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

