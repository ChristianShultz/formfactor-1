#ifndef FAKE_DATA_INI_H_H_GUARD
#define FAKE_DATA_INI_H_H_GUARD

#include "io/adat_xmlio.h"
//#include "radmat/ff/ll_solver.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace ADATXML;
//using namespace ADATIO;

namespace radmat
{
  // fwd
  struct FakeDataIni_t;
  struct StateProps_t;
  struct ZProps_t;
  struct MProps_t;
  struct MatElemProps_t;
  struct TimeProps_t;
  struct pProp_t;
  struct DataProps_t;
  struct DispersionProps_t;
  struct ReadZProps_t;


  void read(XMLReader &, const std::string &, FakeDataIni_t &);
  void read(XMLReader &, const std::string &, StateProps_t &);
  void read(XMLReader &, const std::string &, ZProps_t &);
  void read(XMLReader &, const std::string &, MProps_t &);
  void read(XMLReader &, const std::string &, MatElemProps_t &);
  void read(XMLReader &, const std::string &, TimeProps_t &);
  void read(XMLReader &, const std::string &, pProp_t &);
  void read(XMLReader &, const std::string &, DataProps_t &);
  void read(XMLReader &, const std::string &, DispersionProps_t &);
  void read(XMLReader &, const std::string &, ReadZProps_t &);


  std::string write_params(const FakeDataIni_t &);
  std::string write_params(const StateProps_t &); 
  std::string write_params(const ZProps_t &);
  std::string write_params(const MProps_t &);
  std::string write_params(const MatElemProps_t &);
  std::string write_params(const TimeProps_t &);
  std::string write_params(const pProp_t &);
  std::string write_params(const DataProps_t &);
  std::string write_params(const DispersionProps_t &);
  std::string write_params(const ReadZProps_t &);

  std::ostream& operator<<(std::ostream &, const FakeDataIni_t &);
  std::ostream& operator<<(std::ostream &, const StateProps_t &);
  std::ostream& operator<<(std::ostream &, const MProps_t &);
  std::ostream& operator<<(std::ostream &, const MatElemProps_t &);
  std::ostream& operator<<(std::ostream &, const TimeProps_t &);
  std::ostream& operator<<(std::ostream &, const pProp_t &);
  std::ostream& operator<<(std::ostream &, const DataProps_t &);
  std::ostream& operator<<(std::ostream &, const DispersionProps_t &);
  std::ostream& operator<<(std::ostream &, const ReadZProps_t &);


  //impl

  struct MProps_t
  {
    Array<double> sourceMasses;
    double sourceVarO;
    bool sourceUpdateCovariance;
    bool sourceUpdateVariance;
    Array<double> sinkMasses;
    double sinkVarO;
    bool sinkUpdateCovariance;
    bool sinkUpdateVariance;
  };

  struct ZProps_t
  {
    std::string overlapGenerator;
    bool suppress;
    bool targetZ_at_order1;
    double suppressionOrder;
    bool updateCovariance; // do we want all the Z to have the same timeslice covariance (probably not)
    bool updateVariance;   // do we want all the Z to have the same variance (propbably not)
    double varianceOrder;  // gets multiplied by a vector of normal rands [0,1]
  };

  struct ReadZProps_t
  { 
    Array<double> zsource_r;
    Array<double> zsource_i;
    Array<double> zsink_r;
    Array<double> zsink_i;
  };

  struct StateProps_t 
  {
    bool sameOp;         // is it a diag corr 
    bool readZ;
    MProps_t mProps;
    ZProps_t zProps;
    ReadZProps_t readZProps;
  };

  struct MatElemProps_t
  {
    std::string upper;   // think of the things that go into fake data as a sum over a matrix
    std::string lower;   // Sum_n,m   factor(time) * <n| j | m> -- the decomps are time dependent
    std::string diag;
    int left_target;
    int right_target;
  };

  struct TimeProps_t 
  {
    int tsource;
    int tsink;
  };

  struct pProp_t
  {
    Array<int> momSource;
    Array<int> momSink;
  };

  struct DataProps_t
  {
    int ncfg;
    Array<int> hel_source;
    Array<int> hel_sink;
    Array<pProp_t> momenta;
  };

  struct DispersionProps_t
  {
    double a_t_inverse;
    double xi;
    double L_s;
  };

  struct FakeDataIni_t
  {
    int version;
    StateProps_t stateProps;
    MatElemProps_t matElemProps;
    TimeProps_t timeProps;
    DataProps_t dataProps;
    DispersionProps_t dispersionProps;
  };

}


#endif
