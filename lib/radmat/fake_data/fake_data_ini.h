#ifndef FAKE_DATA_INI_H_H_GUARD
#define FAKE_DATA_INI_H_H_GUARD

#include "io/adat_xmlio.h"
#include "radmat/ff/LLSQ_solver.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace ADATXML;
using namespace ADATIO;

namespace radmat
{
  // fwd
  struct FakeDataIni_t;
  struct StateProps_t;
  struct SuppressionProps_t;
  struct MatElemProps_t;
  struct TimeProps_t;
  struct pProp_t;
  struct DataProps_t;


  void read(XMLReader &, const std::string &, FakeDataIni_t &);
  void read(XMLReader &, const std::string &, StateProps_t &);
  void read(XMLReader &, const std::string &, SuppressionProps_t &);
  void read(XMLReader &, const std::string &, MatElemProps_t &);
  void read(XMLReader &, const std::string &, TimeProps_t &);
  void read(XMLReader &, const std::string &, pProp_t &);
  void read(XMLReader &, const std::string &, DataProps_t &);


  std::string write_params(const FakeDataIni_t &);
  std::string write_params(const StateProps_t &); 
  std::string write_params(const SuppressionProps_t &);
  std::string write_params(const MatElemProps_t &);
  std::string write_params(const TimeProps_t &);
  std::string write_params(const pProp_t &);
  std::string write_params(const DataProps_t &);

  std::ostream& operator<<(std::ostream &, const FakeDataIni_t &);
  std::ostream& operator<<(std::ostream &, const StateProps_t &);
  std::ostream& operator<<(std::ostream &, const SuppressionProps_t &);
  std::ostream& operator<<(std::ostream &, const MatElemProps_t &);
  std::ostream& operator<<(std::ostream &, const TimeProps_t &);
  std::ostream& operator<<(std::ostream &, const pProp_t &);
  std::ostream& operator<<(std::ostream &, const DataProps_t &);

  //impl
  struct StateProps_t 
  {
    int nStates;
    std::string spectrumGenerator;
    std::string overlapGenerator;
    double ladderSpacing;
  };
  
  struct SuppressionProps_t
  {
    bool suppress;
    bool targetZ_at_order1;
    double suppressionOrder;
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
    Array<pProp_t> momenta;
  };

  struct FakeDataIni_t
  {
    int version;
    StateProps_t stateProps;
    SuppressionProps_t suppressionProps;
    MatElemProps_t matElemProps;
    TimeProps_t timeProps;
    DataProps_t dataProps;
  };


}


#endif
