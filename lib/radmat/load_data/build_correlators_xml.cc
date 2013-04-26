/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : build_correlators_xml.cc

* Purpose :

* Creation Date : 25-04-2013

* Last Modified : Fri Apr 26 17:30:55 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "build_correlators_xml.h"


namespace radmat
{

  namespace
  {
    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
      {
        if(ptop.count(path) > 0)
          read(ptop,path,place);
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
            << ", path was empty, exiting" << std::endl;
          exit(1);
        }
      }

  } // namespace anonomyous 

  //! write it to a string
  std::string toString(const ThreePointCorrXMLIni_t &o)
  { 
    std::stringstream ss;
    ss << "continuumMatElemXML = " <<  o.continuumMatElemXML << "\nsource_id " 
      << o.source_id << " sink_id " << o.sink_id  << " maSource = " << o.maSource
      << " maSink = " << o.maSink; 
    return ss.str();
  }

  //! stream it
  std::ostream& operator<<(std::ostream &o, const ThreePointCorrXMLIni_t &p)
  {
    o << toString(p);
    return o;
  }

  //! xml reader
  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrXMLIni_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"continuumMatElemXML",prop.continuumMatElemXML,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"source_id",prop.source_id,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"sink_id",prop.sink_id,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"maSource",prop.maSource,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"maSink",prop.maSink,__PRETTY_FUNCTION__); 
  }

  //! xml writer
  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrXMLIni_t &prop)
  {
    ADATXML::push(xml,path);
    write(xml,"continuumMatElemXML",prop.continuumMatElemXML);
    write(xml,"source_id",prop.source_id);
    write(xml,"sink_id",prop.sink_id);
    write(xml,"maSource",prop.maSource);
    write(xml,"maSink",prop.maSink); 
    ADATXML::pop(xml);
  }

  std::string toString(const ThreePointCorrIni_t &prop)
  {
    std::stringstream ss;
    ss << "threePointCorrXMLIni = " << prop.threePointCorrXMLIni
      << "\nradmatDBProp = "  << prop.radmatDBProp
      << "\nmatElemID = " << prop.matElemID
      << " xi = " << prop.xi << " L_s = " << prop.L_s;  
    return ss.str();
  }

  std::ostream& operator<<(std::ostream &o, const ThreePointCorrIni_t &prop)
  {
    o << toString(prop);
    return o;
  }

  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrIni_t &prop)
  { 
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"threePointCorrXMLIni",prop.threePointCorrXMLIni,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"radmatDBProp",prop.radmatDBProp,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"matElemID",prop.matElemID,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"xi",prop.xi,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"L_s",prop.L_s,__PRETTY_FUNCTION__);
  }

  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrIni_t &prop)
  {
    ADATXML::push(xml,path);
    write(xml,"threePointCorrXMLIni",prop.threePointCorrXMLIni);
    write(xml,"radmatDBProp",prop.radmatDBProp);
    write(xml,"matElemID",prop.matElemID); 
    write(xml,"xi",prop.xi);
    write(xml,"L_s",prop.L_s); 
    ADATXML::pop(xml);
  }

} // namespace radmat

