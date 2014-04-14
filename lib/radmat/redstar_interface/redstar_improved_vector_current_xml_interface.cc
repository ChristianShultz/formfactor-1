/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : redstar_improved_vector_current_xml_interface.cc

* Purpose :

* Creation Date : 14-04-2014

* Last Modified : Mon 14 Apr 2014 12:15:02 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_improved_vector_current_xml_interface.h"
#include "redstar_photon_props.h"


namespace radmat
{

  namespace
  {

    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, 
          const std::string &path, 
          T &place, 
          const char * f)
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

    std::string toString(const RedstarImprovedVectorCurrentXML::insertion &i)
    {
      std::stringstream ss; 
      ss << "active= " << i.active << " create= " << PHOTON_CREATE
        << " smear= " << i.smearedP << " photons: ";
      for(int j = 0; j < i.photons.size(); ++j)
        ss << "(" << i.photons[j].coeff_r 
          << " + " << i.photons[j].coeff_i 
          << "i) x" << i.photons[j].name << "   ";
      return ss.str(); 
    }


  } // anonomyous  




  void 
    RedstarImprovedVectorCurrentXML::read(ADATXML::XMLReader &xml, const std::string &path)
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"pmin",pmin,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"pmax",pmax,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"t_slice",t_slice,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"time",time,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"space",space,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"improvement",imp,__PRETTY_FUNCTION__); 


      if ( time.active && space.active ) 
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__
          << ": Warning, you decided to mix time and" 
          << " space and I don't think you should " << std::endl; 
      }
    }

  std::string 
    RedstarImprovedVectorCurrentXML::write(void) const
    {
      std::stringstream ss;
      ss << "pmin= " << pmin << " pmax=" << pmax << " t_slice= " << t_slice; 
      ss << "\ntime:\n" << toString(time) << std::endl;
      ss << "\nspace:\n" << toString(space) << std::endl;
      return ss.str(); 
    }

    void write(ADATXML::XMLWriter &xml,
        const std::string &path, 
        const RedstarImprovedVectorCurrentXML::improvement &i)
    {
      ADATXML::push(xml,path);
      ADATXML::write(xml,"Omega_s",i.Omega_s);
      ADATXML::write(xml,"Omega_t",i.Omega_t); 
      ADATXML::write(xml,"name",i.name); 
      ADATXML::pop(xml);
    }

    void read(ADATXML::XMLReader &xml, 
        const std::string &path,
        RedstarImprovedVectorCurrentXML::improvement &i)  
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"Omega_s",i.Omega_s,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"Omega_t",i.Omega_t,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"name",i.name,__PRETTY_FUNCTION__); 
    }

} // radmat
