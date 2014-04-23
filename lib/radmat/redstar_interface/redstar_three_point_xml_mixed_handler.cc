/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_three_point_xml_mixed_handler.cc

* Purpose :

* Creation Date : 23-04-2014

* Last Modified : Wed 23 Apr 2014 11:23:06 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_three_point_xml_mixed_handler.h"
#include "redstar_three_point_utils.h"
#include "radmat/utils/utils.h"

namespace radmat
{



  std::vector<ThreePointData>
    RedstarThreePointXMLMixedHandler::handle_work() 
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

      std::vector<BlockData> left = generate_cubic_block(lefty);
      std::vector<BlockData> gamma = generate_lorentz_block(gammay);
      std::vector<BlockData> right = generate_cubic_block(righty); 

      ADATXML::Array<int> timesliz(3); 
      timesliz[0] = lefty->t_slice; 
      timesliz[1] = gammay->t_slice; 
      timesliz[2] = righty->t_slice; 

      this->set_timeslice(timesliz); 

      return merge_blocks( left, gamma, right, ensemble ); 
    }


} // radmat 


