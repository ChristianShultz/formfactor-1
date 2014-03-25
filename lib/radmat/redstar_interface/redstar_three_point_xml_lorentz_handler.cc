/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_three_point_xml_lorentz_handler.cc

* Purpose :

* Creation Date : 20-03-2014

* Last Modified : Mon 24 Mar 2014 02:41:45 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_three_point_xml_lorentz_handler.h"
#include "redstar_three_point_utils.h"
#include "radmat/utils/utils.h"

namespace radmat
{



  std::vector<ThreePointData>
    RedstarThreePointXMLLorentzHandler::handle_work(const std::string &left_id,
        const std::string &right_id) const
    {
      ADATXML::Array<AbstractNamedObject<AbsRedstarXMLInterface_t> > xml_int; 
      std::string ensemble; 

      xml_int = this->get_npt(); 
      ensemble = this->get_ensemble(); 

      const RedstarSingleParticleMesonXML * lefty;
      const RedstarSingleParticleMesonXML * righty;
      const RedstarVectorCurrentXML *gammay; 

      POW2_ASSERT( xml_int[0].param->type() == Stringify<RedstarSingleParticleMesonXML>() );      
      POW2_ASSERT( xml_int[1].param->type() == Stringify<RedstarVectorCurrentXML>() );      
      POW2_ASSERT( xml_int[2].param->type() == Stringify<RedstarSingleParticleMesonXML>() );      

      lefty = dynamic_cast< const RedstarSingleParticleMesonXML * >( xml_int[0].param.get_ptr() ); 
      gammay = dynamic_cast< const RedstarVectorCurrentXML * >( xml_int[1].param.get_ptr() ); 
      righty = dynamic_cast< const RedstarSingleParticleMesonXML * >( xml_int[2].param.get_ptr() ); 

      std::vector<BlockData> left = generate_lorentz_block(lefty);
      std::vector<BlockData> gamma = generate_lorentz_block(gammay);
      std::vector<BlockData> right = generate_lorentz_block(right); 

      ADATXML::Array<int> timesliz(3); 
      timesliz[0] = left->t_slice; 
      timesliz[1] = gamma->t_slice; 
      timesliz[2] = right->t_slice; 

      this->set_timeslice(timesliz); 

      return merge_blocks( left, gamma, right, ensemble ); 
    }


} // radmat 

