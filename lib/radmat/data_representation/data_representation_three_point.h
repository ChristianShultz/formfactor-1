#ifndef DATA_REPRESENTATION_THREE_POINT_H
#define DATA_REPRESENTATION_THREE_POINT_H 


#include "data_representation_primitive_rep.h"
#include "data_representation_factory.h"
#include "data_representation_cubic_groups.h"
#include "data_representation_lorentz_groups.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/stringify.h"
#include <string>


namespace radmat
{

  // where does a matrix element live 
  struct DataRep3pt
  {
    DataRep3pt( const std::string &ll,
        const std::string &gg, 
        const std::string &rr)
      : l(ll) , r(rr) , g(gg) 
    { }

    virtual ~DataRep3pt() {}

    virtual rHandle<Rep_p> lefty() const { return call_rep_prim(l); }
    virtual rHandle<Rep_p> righty() const { return call_rep_prim(r); }
    virtual rHandle<Rep_p> gamma() const { return call_rep_prim(g); }

    virtual rHandle<Rep_p> call_rep_prim(const std::string &s) const 
    {
      return DataRepresentationFactoryEnv::callFactory(s); 
    }

    virtual rHandle<LorentzRep> call_rep_lor(const std::string &s) const
    {
      return LorentzRepresentationFactoryEnv::callFactory(s); 
    }

    virtual rHandle<CubicRep> call_rep_cub(const std::string &s) const 
    {
      return CubicRepresentationFactoryEnv::callFactory(s); 
    }

    bool is_cubic(const std::string &s)
    {
      rHandle<Rep_p> r = call_rep_prim(s); 
      return r->rep_type() == Stringify<CubicRep_t>(); 
    }

    bool is_lorentz(const std::string &s)
    {
      rHandle<Rep_p> r = call_rep_prim(s); 
      return r->rep_type() == Stringify<LorentzRep_t>(); 
    }

    std::string l,r,g; 
  };

}


#endif /* DATA_REPRESENTATION_THREE_POINT_H */
