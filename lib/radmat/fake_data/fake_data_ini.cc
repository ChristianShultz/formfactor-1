// fake_data_ini.cc -
//
// Wednesday, June 20 2012
//

#include"fake_data_ini.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"
#include <string>
#include <iostream>
#include <sstream>
#include <complex>

using namespace ADATXML;
using namespace std;

namespace radmat
{

  void read(XMLReader &xml, const std::string &path, FakeDataIni_t &prop)
  {
    XMLReader ptop(xml,path);
    const int version = 1;

    //check version
    if(ptop.count("version") > 0)
      read(ptop,"version",prop.version);
    else
    {
      SPLASH("No version present, trying to proceed");
      prop.version = version;
    }

    if(version != prop.version)
    {
      std::cout << "warning: version " << prop.version << " is not supported. The supported version is "
        << version << std::endl;
      POW2_ASSERT(false);
    }
    if(ptop.count("stateProps") > 0)
      read(ptop,"stateProps",prop.stateProps);
    else
    {
      SPLASH("No stateProps present, exiting");
      exit(1);
    }

    if(ptop.count("matElemProps") > 0)
      read(ptop,"matElemProps",prop.matElemProps);
    else
    {
      SPLASH("No matElemProps present, exiting");
      exit(1);
    }

    if(ptop.count("timeProps") > 0)
      read(ptop,"timeProps", prop.timeProps);
    else
    {
      SPLASH("No timeProps present, exiting");
      exit(1);
    }

    if(ptop.count("dataProps") > 0)
      read(ptop,"dataProps", prop.dataProps);
    else
    {
      SPLASH("No dataProps present, exiting");
      exit(1);
    }

    if(ptop.count("dispersionProps") > 0)
      read(ptop,"dispersionProps",prop.dispersionProps);
    else
    {
      SPLASH("No dispersionProps present, exiting");
      exit(1);
    }
  }

  void read(XMLReader &xml, const std::string &path, MProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("sourceMasses") > 0)
      read(ptop,"sourceMasses",prop.sourceMasses);
    else
    {
      SPLASH("No source masses present");
      exit(1);
    }

    if(ptop.count("sourceVarO") > 0)
      read(ptop,"sourceVarO",prop.sourceVarO);
    else
      prop.sourceVarO = 0.1;

    if(ptop.count("sourceUpdateCovariance") > 0)
      read(ptop,"sourceUpdateCovariance",prop.sourceUpdateCovariance);
    else
      prop.sourceUpdateCovariance = true;

    if(ptop.count("sourceUpdateVariance") > 0)
      read(ptop,"sourceUpdateVariance",prop.sourceUpdateVariance);
    else
      prop.sourceUpdateVariance = true;

    if(ptop.count("sinkMasses") > 0) 
      read(ptop,"sinkMasses",prop.sinkMasses);
    else
    {
      SPLASH("No sink masses present");
      exit(1);
    }

    if(ptop.count("sinkVarO") > 0)
      read(ptop,"sinkVarO",prop.sinkVarO);
    else
      prop.sinkVarO = 0.1;

    if(ptop.count("sinkUpdateCovariance") > 0)
      read(ptop,"sinkUpdateCovariance",prop.sinkUpdateCovariance);
    else
      prop.sinkUpdateCovariance = true;

    if(ptop.count("sinkUpdateVariance") > 0)
      read(ptop,"sinkUpdateVariance",prop.sinkUpdateVariance);
    else
      prop.sinkUpdateVariance = true;
  } 


  void read(XMLReader &xml, const std::string &path, ZProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("overlapGenerator") > 0)
      read(ptop,"overlapGenerator",prop.overlapGenerator);
    else
      prop.overlapGenerator = std::string("unitary");

    if(ptop.count("suppress") > 0)
      read(ptop,"suppress",prop.suppress);
    else
      prop.suppress = true;

    if(ptop.count("targetZ_at_order1") > 0)
      read(ptop,"targetZ_at_order1",prop.targetZ_at_order1);
    else
      prop.targetZ_at_order1 = true;

    if(ptop.count("suppressionOrder") > 0)
      read(ptop,"suppressionOrder",prop.suppressionOrder);
    else
      prop.suppressionOrder = 1e-2;

    if(ptop.count("updateCovariance") > 0)
      read(ptop,"updateCovariance",prop.updateCovariance);
    else
      prop.updateCovariance = true;

    if(ptop.count("updateVariance") > 0)
      read(ptop,"updateVariance",prop.updateVariance);
    else
      prop.updateVariance = true;

    if(ptop.count("varianceOrder") > 0)
      read(ptop,"varianceOrder", prop.varianceOrder);
    else
      prop.varianceOrder = 0.1;
  }

  void read(XMLReader &xml, const std::string &path, ReadZProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("zsource_r") > 0)
      read(ptop,"zsource_r",prop.zsource_r);

    if(ptop.count("zsource_i") > 0)
      read(ptop,"zsource_i",prop.zsource_i);

    if(ptop.count("zsink_r") > 0)
      read(ptop,"zsink_r",prop.zsink_r);

    if(ptop.count("zsink_i") > 0)
      read(ptop,"zsink_i",prop.zsink_i);

  }

  void read(XMLReader &xml, const std::string &path, StateProps_t &prop)
  { 
    XMLReader ptop(xml,path);

    if(ptop.count("sameOp") > 0)
      read(ptop,"sameOp",prop.sameOp);
    else
      prop.sameOp = true;

    if(ptop.count("readZ") > 0)
      read(ptop,"readZ",prop.readZ);
    else
      prop.readZ = false;

    if(ptop.count("mProps") > 0)
      read(ptop,"mProps",prop.mProps);
    else
    {
      SPLASH("no mProps present, exiting");
      exit(1);
    } 

    if(ptop.count("zProps") > 0)
      read(ptop,"zProps",prop.zProps);
    else
    {
      SPLASH("no zProps present, exiting");
      exit(1);
    }

    if(prop.readZ)
    {
      if(ptop.count("readZProps") > 0)
        read(ptop,"readZProps",prop.readZProps);
      else
      {
        SPLASH("no zProps present, exiting");
        exit(1);
      }

      // sanity
      POW2_ASSERT(prop.readZProps.zsource_r.size() == prop.readZProps.zsource_i.size());
      POW2_ASSERT(prop.readZProps.zsink_r.size() == prop.readZProps.zsink_i.size());
      POW2_ASSERT(prop.readZProps.zsource_r.size() == prop.mProps.sourceMasses.size());
      POW2_ASSERT(prop.readZProps.zsink_r.size() == prop.mProps.sinkMasses.size());
    }

  }



  void read(XMLReader &xml, const std::string &path, MatElemProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("off") > 0 )
      read(ptop,"off",prop.off);
    else
    {
      SPLASH("need to provide decompositions");
      exit(1);
    }

    if(ptop.count("diag") > 0)
      read(ptop,"diag",prop.diag);
    else
    {
      SPLASH("need to provide decomposition");
      exit(1);
    }

    if(ptop.count("left_target") > 0)
      read(ptop,"left_target",prop.left_target);
    else
      prop.left_target = 0;

    if(ptop.count("right_target") > 0)
      read(ptop,"right_target",prop.right_target);
    else
      prop.right_target = 0;

  }

  void read(XMLReader &xml, const std::string &path, TimeProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("tsource") > 0)
      read(ptop,"tsource",prop.tsource);
    else
      prop.tsource = 0;

    if(ptop.count("tsink") > 0)
      read(ptop,"tsink",prop.tsink);
    else
      prop.tsink = 25;
  }

  void read(XMLReader &xml, const std::string &path, pProp_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("momSource")  > 0)
      read(ptop,"momSource",prop.momSource);
    else
    {
      prop.momSource.resize(3);
      prop.momSource[0] = prop.momSource[1] = prop.momSource[2] = 0;
    }

    if(ptop.count("momSink") > 0)
      read(ptop,"momSink",prop.momSink);
    else
    {
      prop.momSink.resize(3);
      prop.momSink[0] = prop.momSink[1] = prop.momSink[2] = 0;
    }
  }

  void read(XMLReader &xml, const std::string &path, DataProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("ncfg") > 0)
      read(ptop,"ncfg",prop.ncfg);
    else
      prop.ncfg = 500;

    if(ptop.count("hel_source") > 0)
      read(ptop,"hel_source",prop.hel_source);
    else
    {
      SPLASH("no source helicity provided, barfing");
      exit(1);
    }

    if(ptop.count("hel_sink") > 0)
      read(ptop,"hel_sink",prop.hel_sink);
    else
    {
      SPLASH("no sink helicity provided, barfing");
      exit(1);
    }

    if(ptop.count("momenta") > 0)
      read(ptop,"momenta",prop.momenta);
    else
    {
      SPLASH("no momenta provided, barfing");
      exit(1);
    }
  }

  void read(XMLReader &xml, const std::string &path, DispersionProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("a_t_inverse") > 0)
      read(ptop,"a_t_inverse",prop.a_t_inverse);
    else
      prop.a_t_inverse = 1./5.6;

    if(ptop.count("xi") > 0)
      read(ptop,"xi",prop.xi);
    else
      prop.xi = 3.5;

    if(ptop.count("L_s") > 0)
      read(ptop,"L_s",prop.L_s);
    else
      prop.L_s = 2.;
  }

  namespace
  {
    std::string n = std::string("\n");
  }

  std::string write_params(const FakeDataIni_t &ini)
  {
    std::stringstream ss;
    ss << "version: " << ini.version << n;
    ss << "stateProps: " << n << write_params(ini.stateProps) << n;
    ss << "matElemProps: " << n << write_params(ini.matElemProps) << n;
    ss << "timeProps: " << n << write_params(ini.timeProps) << n;
    ss << "dataProps: " << n << write_params(ini.dataProps) << n;
    ss << "dispersionProps: " << n << write_params(ini.dispersionProps) << n;
    return ss.str();
  }


  std::string write_params(const MProps_t &ini)
  {

    std::stringstream ss;

    ss << "source spectrum:";    
    for(int i = 0; i < ini.sourceMasses.size(); i++)
      ss << " " << ini.sourceMasses[i];
    ss << n;

    ss << "sourceVarO: " << ini.sourceVarO << n;
    ss << "sourceUpdateCovariance: " << ini.sourceUpdateCovariance << n;
    ss << "sourceUpdateVariance: " << ini.sourceUpdateVariance << n;

    ss << "sink spectrum:";    
    for(int i =0 ; i < ini.sinkMasses.size(); i++)
      ss << " " << ini.sinkMasses[i];
    ss << n;

    ss << "sinkVarO: " << ini.sinkVarO << n;
    ss << "sinkUpdateCovariance: " << ini.sinkUpdateCovariance << n;
    ss << "sinkUpdateCovariance: " << ini.sinkUpdateCovariance << n;

    return ss.str();
  }

  std::string write_params(const ZProps_t &ini)
  {
    std::stringstream ss;

    ss << "overlapGenerator: " << ini.overlapGenerator << n;
    ss << "suppress: " << ini.suppress << n;
    ss << "targetZ_at_order1: " << ini.targetZ_at_order1 << n;
    ss << "suppressionOrder: " << ini.suppressionOrder << n;
    ss << "updateCovariance: " << ini.updateCovariance << n;
    ss << "updateVariance: " << ini.updateVariance << n;
    ss << "varianceOrder: " << ini.varianceOrder << n;

    return ss.str();
  }

  std::string write_params(const ReadZProps_t &ini)
  {
    std::stringstream ss;

    ss << "sourceZ:";
    for(int i = 0; i < ini.zsource_r.size(); i++)
      ss << " " << std::complex<double>(ini.zsource_r[i],ini.zsource_i[i]);
    ss << n;

    ss << "sinkZ:";
    for(int i = 0; i < ini.zsink_r.size(); i++)
      ss << " " << std::complex<double>(ini.zsink_r[i],ini.zsink_i[i]);
    ss << n;

    return ss.str();
  }

  std::string write_params(const StateProps_t &ini)
  {
    std::stringstream ss;

    ss << "sameOp: " << ini.sameOp << n;
    ss << "readZ: " << ini.readZ << n;
    ss << "mProps: " << write_params(ini.mProps) << n;
    ss << "zProps: " << write_params(ini.zProps) << n;
    if(ini.readZ)
      ss << "readZProps: " << write_params(ini.readZProps); 

    return ss.str();
  }


  std::string write_params(const MatElemProps_t &ini)
  {
    std::stringstream ss;
    ss << "off: " << ini.off << n;
    ss << "diag: " << ini.diag << n;
    ss << "left_target: " << ini.left_target << n;
    ss << "right_target: " << ini.right_target << n;

    return ss.str();
  }

  std::string write_params(const TimeProps_t &ini)
  {
    std::stringstream ss;
    ss << "tsource: " << ini.tsource << n;
    ss << "tsink: " << ini.tsink << n;

    return ss.str();
  }

  std::string write_params(const pProp_t &ini)
  {
    std::stringstream ss;
    ss << "momSource: " << ini.momSource[0] << " " << ini.momSource[1] << " " << ini.momSource[2] << n;
    ss << "momSink:   " << ini.momSink[0] << " " << ini.momSink[1] << " " << ini.momSink[2] << n;

    return ss.str();
  }

  std::string write_params(const DataProps_t &ini)
  {
    std::stringstream ss;

    ss << "ncfg: " << ini.ncfg << n;
    ss << "hel_source: ";
    for(int i = 0; i < ini.hel_source.size(); i++)
      ss << ini.hel_source[i] << " ";
    ss << n;

    ss << "hel_sink: ";
    for(int i =0; i < ini.hel_sink.size(); i++)
      ss << ini.hel_sink[i] << " ";
    ss << n;

    for(int i = 0; i < ini.momenta.size(); i++)
      ss << "elem: \n" << write_params(ini.momenta[i]) << n;

    return ss.str();
  }  

  std::string write_params(const DispersionProps_t &ini)
  {
    std::stringstream ss;

    ss << "a_t_inverse: " << ini.a_t_inverse << n;
    ss << "xi: " << ini.xi << n;
    ss << "L_s: " << ini.L_s << n;

    return ss.str();
  }

  std::ostream& operator<<(std::ostream &o , const FakeDataIni_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const StateProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const MProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const ZProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const ReadZProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const MatElemProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const TimeProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const DataProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const pProp_t &ini)
  {
    o << write_params(ini);
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const DispersionProps_t &ini)
  {
    o << write_params(ini);
    return o;
  }

}
