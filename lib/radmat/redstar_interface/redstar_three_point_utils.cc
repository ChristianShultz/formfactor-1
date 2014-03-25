/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_three_point_utils.cc

 * Purpose :

 * Creation Date : 21-03-2014

 * Last Modified : Fri 21 Mar 2014 01:58:23 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_three_point_utils.h"


namespace radmat
{
  namespace
  {
    // returns true if momentum is conserved
    bool check_mom(const BlockData &l, const BlockData &g, const BlockData &r)
    {
      ADATXML::Array<int> ll,gg,rr;
      ll = l.data.begin()->m_obj.irrep.mom;
      gg = g.data.begin()->m_obj.irrep.mom;
      rr = r.data.begin()->m_obj.irrep.mom;

      int ls(1),gs(1),rs(1); 

      if( !!! l.data.begin()->m_obj.irrep.creation_op )
        ls = - 1 ; 

      if( !!! g.data.begin()->m_obj.irrep.creation_op )
        gs = - 1 ; 

      if( !!! r.data.begin()->m_obj.irrep.creation_op )
        rs = - 1 ; 
  
      for(int i =0; i < 3; ++i)
        if( ls*ll[i] + gs*gg[i] + rs*rr[i] != 0 )
         return false;  
    
      return true; 
    }

    EnsemRedstarNPtBlock 
      merge_ensem_blocks( const EnsemRedstarBlock &lefty,
          const EnsemRedstarBlock &gamma,
          const EnsemRedstarBlock &righty, 
          const std::string &ensemble)
      { 
        EnsemRedstarNPtBlock ret; 
        EnsemRedstarBlock::const_iterator l,g,r;

        for( l = lefty.begin(); l != lefty.end(); ++l)
         for( g = gamma.begin(); g != gamma.end(); ++g)
          for( r = righty.begin(); r != righty.end(); ++r)  
          {
            ENSEM::Complex coeff;
            Hadron::KeyHadronNPartNPtCorr_t npt; 
            coeff = l->m_coeff * g->m_coeff * r->m_coeff; 
            npt.npoint.resize(3); 

            // arrays are FORTRAN style 
            npt.npoint[1] = l->m_obj;
            npt.npoint[2] = g->m_obj;
            npt.npoint[3] = r->m_obj;
            npt.ensemble = ensemble; 

            ret = ret + EnsemRedstarNPtBlock::ListObj_t(coeff,npt); 
          }

        return ret; 
      } 


    ThreePointData
      make_data(const BlockData &l, 
          const BlockData &g, 
          const BlockData &r, 
          const std::string &ensemble)
      {
        ThreePointData ret; 
        ret.origin_rep = DataRep3pt(l.origin_rep,
            g.origin_rep,
            r.origin_rep);
        ret.data_rep = DataRep3pt(l.data_rep,
            g.data_rep,
            r.data_rep);
        ret.left_row = l.row;
        ret.gamma_row = g.row;
        ret.right_row = r.row; 

        ret.data = merge_ensem_blocks( l.data,g.data,r.data,ensemble);

        return ret; 
      }


  } // anonomyous 

  std::vector<ThreePointData>
    merge_blocks(const std::vector<BlockData> &lefty, 
        const std::vector<BlockData> &gamma,
        const std::vector<BlockData> &righty,
        const std::string &ensemble)
    {
      std::vector<ThreePointData> ret;
      std::vector<BlockData>::const_iterator l,r,g;

      // some poor estimate of the size, lots of them are removed via 
      // momentum conservation
      ret.reserve( lefty.size() * righty.size() / 2 ); 

      for( l = lefty.begin(); l != lefty.end(); ++l)
        for( g = gamma.begin(); g != gamma.end(); ++g)
          for(r = righty.begin(); r != righty.end(); ++r)
          {
            if( check_mom(*l,*g,*r) )
              all.push_back( make_data(*l,*g,*r,ensemble) ); 
          }

      return ret; 
    }


}

