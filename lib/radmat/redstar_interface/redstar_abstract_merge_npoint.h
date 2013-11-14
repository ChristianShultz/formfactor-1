#ifndef REDSTAR_ABSTRACT_MERGE_NPOINT_H
#define REDSTAR_ABSTRACT_MERGE_NPOINT_H 

#include "redstar_npoint_function_xml.h"
#include "redstar_abstract_block_base.h"
#include "radmat/utils/obj_expr_t.h"
#include "radmat/utils/stringify.h"
#include "hadron/hadron_npart_npt_corr.h"
#include <algorithm>


namespace radmat
{

  // a coefficient times an n-point function
  //    this basically corresponds to some euclidean 
  //    continuum n-point now 
  typedef ListObjExpr_t<ENSEM::Complex,
          Hadron::KeyHadronNPartNPtCorr_t> EnsemRedstarNPtBlock; 


  // the merged npoint data store
  struct AbsRedstarMergeNPtData_t
  {
    AbsRedstarMergeNPtData_t(void) {}
    AbsRedstarMergeNPtData_t( const AbsRedstarMergeNPtData_t &o) 
    {
      npoint = o.npoint; 
      input = o.input; 
    }

    AbsRedstarMergeNPtData_t& operator=(const AbsRedstarMergeNPtData_t &o)
    {
      if (this != &o) 
      { 
        AbsRedstarMergeNPtData_t foo(o); 
        std::swap(npoint,foo.npoint); 
        std::swap(input,foo.input); 
      }

      return *this; 
    }


    std::vector<EnsemRedstarNPtBlock> npoint; 
    std::vector<std::vector< ADAT::Handle< AbsRedstarInput_t > > > input; 

  };



  // an abstract polymorphic base to tell us what to do with the NPointXML
  struct AbsRedstarMergeNPt
  {
    virtual std::string 
      type(void) const = 0; 

    // do something with my_npt and populate my_data
    virtual void 
      do_work(void) = 0; 


    virtual std::vector<EnsemRedstarNPtBlock> 
      xml(void) const {return my_data.npoint;} 
    virtual AbsRedstarMergeNPtData_t 
      data(void) const {return my_data;}
    virtual void 
      read(ADATXML::XMLReader &xml, const std::string &path) 
      {
        ::radmat::read(xml,path,my_npt); 
      }

    // what timeslices do each guy sit on
    virtual ADATXML::Array<int> timeslice_info(void) const
    {
      ADATXML::Array<int> t(my_npt.N);
      for(int i = 0; i < my_npt.N; ++i)
        t[i] = my_npt.npoint[i].param->timeslice(); 
      return t;
    } 

    NPointXML my_npt;
    AbsRedstarMergeNPtData_t my_data; 
  }; 


} // radmat 


#endif /* REDSTAR_ABSTRACT_MERGE_NPOINT_H */
