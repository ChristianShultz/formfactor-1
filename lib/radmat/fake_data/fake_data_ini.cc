// fake_data_ini.cc -
//
// Wednesday, June 20 2012
//

#include"fake_data_ini.h"
#include "io/adat_xmlio.h"
//#include "radmat/ff/ll_solver.h"
#include "radmat/utils/splash.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace ADATXML;
//using namespace ADATIO;
using namespace std;

namespace radmat
{

  void read(XMLReader &xml, const std::string &path, FakeDataIni_t &prop)
  {
    XMLReader ptop(xml,path);
    const int version = 0;
    
    //check version
    if(ptop.count("version") > 0)
      read(ptop,"version",prop.version);
    else
      {
	SPLASH("No version present, trying to proceed");
	prop.version = version;
      }

    if(version != prop.version)
      std::cout << "warning: version " << prop.version << " is not supported. The supported version is " << version << std::endl;
	
    if(ptop.count("stateProps") > 0)
      read(ptop,"stateProps",prop.stateProps);
    else
      {
	SPLASH("No stateProps present, exiting");
	exit(1);
      }

    if(ptop.count("suppressionProps") > 0)
      read(ptop,"suppressionProps",prop.suppressionProps);
    else
      {
	SPLASH("No suppressionProps present, exiting");
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
  }

  void read(XMLReader &xml, const std::string &path, StateProps_t &prop)
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

    if(ptop.count("overlapGenerator") > 0)
      read(ptop,"overlapGenerator",prop.overlapGenerator);
    else
      prop.overlapGenerator = std::string("unitary");
  } 
  

  void read(XMLReader &xml, const std::string & path, SuppressionProps_t &prop)
  {
    XMLReader ptop(xml,path);

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
  }


  void read(XMLReader &xml, const std::string &path, MatElemProps_t &prop)
  {
    XMLReader ptop(xml,path);

    if(ptop.count("upper") > 0 )
      read(ptop,"upper",prop.upper);
    else
      {
	SPLASH("need to provide decompositions");
	exit(1);
      }

    if(ptop.count("lower") > 0)
      read(ptop,"lower",prop.lower);
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

    if(ptop.count("sameOp") > 0)
      read(ptop,"sameOp",prop.sameOp);
    else
      prop.sameOp = true;

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

    if(ptop.count("momenta") > 0)
      read(ptop,"momenta",prop.momenta);
    else
      {
	SPLASH("no momenta provided, barfing");
	exit(1);
      }
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
    ss << "suppressionProps: " << n << write_params(ini.suppressionProps) << n;
    ss << "matElemProps: " << n << write_params(ini.matElemProps) << n;
    ss << "timeProps: " << n << write_params(ini.timeProps) << n;
    ss << "dataProps: " << n << write_params(ini.dataProps) << n;
    return ss.str();
  }


  std::string write_params(const StateProps_t &ini)
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
    ss << "overlapGenerator: " << ini.overlapGenerator << n;

    return ss.str();
  }

  std::string write_params(const SuppressionProps_t &ini)
  {
    std::stringstream ss;
    ss << "suppress: " << ini.suppress << n;
    ss << "targetZ_at_order1: " << ini.targetZ_at_order1 << n;
    ss << "suppressionOrder: " << ini.suppressionOrder << n;
    
    return ss.str();
  }

  std::string write_params(const MatElemProps_t &ini)
  {
    std::stringstream ss;
    ss << "upper: " << ini.upper << n;
    ss << "lower: " << ini.lower << n;
    ss << "diag: " << ini.diag << n;
    ss << "left_target: " << ini.left_target << n;
    ss << "right_target: " << ini.right_target << n;
    ss << "sameOp: " << ini.sameOp << n;
    ss << "updateCovariance: " << ini.updateCovariance << n;
    ss << "updateVariance: " << ini.updateVariance << n;
    ss << "varianceOrder: " << ini.varianceOrder << n;

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

    for(int i = 0; i < ini.momenta.size(); i++)
      ss << "elem: \n" << write_params(ini.momenta[i]) << n;

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

  std::ostream& operator<<(std::ostream &o, const SuppressionProps_t &ini)
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


}
