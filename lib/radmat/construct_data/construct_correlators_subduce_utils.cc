/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : construct_correlators_subduce_utils.cc

 * Purpose :

 * Creation Date : 14-03-2014

 * Last Modified : Mon 17 Mar 2014 03:20:45 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "construct_correlators_subduce_utils.h"
#include "radmat/utils/printer.h"
#include "radmat/ff_interface/formfactor_subduced_formfactors.h"
#include "hadron/ensem_filenames.h"
#include <algorithm>

namespace radmat
{

  namespace 
  {
    //
    //  a bunch of printers, turn them on to peek at the algorithm 
    //

    struct avail_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << msg << std::endl; }
    };


    struct subduce_print_1
    {
      static void print(const std::string &msg)
      { std::cout << "print1 " << msg << std::endl;}
    };

    struct subduce_print_2
    {
      static void print(const std::string &msg)
      { std::cout << "print2 " << msg << std::endl;}
    };

    //
    // real work starts herrr
    //

    // what are all the continuum matrix elements we were passed? 
    std::map<std::string,TaggedEnsemRedstarNPtBlock>
      make_continuum_map(const std::vector<TaggedEnsemRedstarNPtBlock> & cont)
      {
        std::map<std::string,TaggedEnsemRedstarNPtBlock> ret; 
        std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator it; 
        EnsemRedstarNPtBlock::const_iterator block; 

        for(it = cont.begin(); it != cont.end(); ++it)
        {
          printer_function<avail_printer>(it->continuum_tag.unique_tag());
          ret.insert(std::make_pair(it->continuum_tag.unique_tag(),*it)); 
        }

        return ret; 
      }

    // string magic  -- split the thing on delim like in perl 
    std::vector<std::string> tokenize(const std::string &s,
        const std::string &delim, 
        const bool keep_empty_tokens=false)
    {
      std::vector<std::string> tokens;

      if(delim.empty())
        return std::vector<std::string>(1,s);

      std::string::const_iterator front,back;
      front = s.begin();
      const std::string::const_iterator end = s.end();

      while(true)
      {

        back = std::search(front,end,delim.begin(),delim.end());
        std::string token(front,back);
        if(keep_empty_tokens || !!!token.empty())
          tokens.push_back(token);

        if(back == end)
          break;

        front = back + delim.size();
      }

      return tokens; 
    }


    // a quickie for a continuum representation 
    struct cont_rep
    {
      cont_rep() {}
      cont_rep(const std::string &s, const ADATXML::Array<int> &p, const int r)
        : rep(s) , mom(p) , row(r)
      { }

      std::string rep;
      ADATXML::Array<int> mom;
      int row; 
    };


    // hold three of these guys corresponding to each 
    // operator in the three point
    struct rep_holder
    { 
      std::map<std::string,cont_rep> left,right,gamma; 
    };

    // the keys of a map are unique, this automatically 
    // pulls only the unique combinations 
    //
    // we are looping over all of the cont strings and pulling out something 
    // that looks like J1p,p000,r1 
    rep_holder
      get_reps(const std::map<std::string,TaggedEnsemRedstarNPtBlock> &b)
      {
        rep_holder ret; 

        std::map<std::string,TaggedEnsemRedstarNPtBlock>::const_iterator block_it; 
        for(block_it = b.begin(); block_it != b.end(); ++block_it)
        {
          std::vector<std::string> tokens = tokenize(block_it->first,".");
          if( tokens.size() != 3 )
          {
            std::cout << __PRETTY_FUNCTION__ << ": token error" << std::endl;
            exit(1);
          }

          ret.left[tokens[0]] = cont_rep( tokens[0] , 
              block_it->second.continuum_tag.p_f,
              block_it->second.continuum_tag.hf);
          ret.gamma[tokens[1]] = cont_rep( tokens[1] ,
              block_it->second.continuum_tag.q,
              block_it->second.continuum_tag.jmu);
          ret.right[tokens[2]] = cont_rep( tokens[2] ,
              block_it->second.continuum_tag.p_i,
              block_it->second.continuum_tag.hi);
        }

        return ret; 
      }



    // look into our subduction table map and see if we know about this fool  
    std::vector<const SubduceTableMap::irrep_sub_table*> 
      possible_subductions(const std::string &J, const ADATXML::Array<int> &mom)
      {
        std::vector<const SubduceTableMap::irrep_sub_table*> ret; 
        typedef TheSmarterSubduceTableMap mmap;
        SubduceTableMap::map_t::const_iterator it; 
        for(it = mmap::Instance().mappy.begin(); it != mmap::Instance().mappy.end(); ++it)
        {
          RepHandler foo; 
          rHandle<Rep_p> mom_rep = foo.gen_rep(mom); 

          // this is a key from the subduce table map , it 
          // looks something like J1mT1 or J0pH0D2A1
          // for the case of rest or flight respectively 
          std::string key = it->first;  

          // we are matching on the J#p part of the key 

          // this is like =~ in perl 
          if(key.find(J) != std::string::npos) 
          {
            // now another check, 
            //
            // do they live in the same symmetry group (Oh, D2 etc)
            //   subduce table also has the irrep embedded in it 
            //   so we can deal with it one level higher and loop over the 
            //   possible stuff (rows and irreps)
            const SubduceTableMap::irrep_sub_table* tab = it->second; 
            if( tab->lat->rep_g() == mom_rep->id() )
              ret.push_back(tab); 
          }
        }

        return ret;
      }

    // a quick struct for the subduced rep 
    struct subduce_rep 
    {
      subduce_rep() {}
      subduce_rep(const std::string &t, 
          const std::string &c, 
          const ADATXML::Array<int> &p, 
          const int r)
        : cont_token(t) , cubic_rep(c) , mom(p) , row(r)
      { }

      std::string cont_token; 
      std::string cubic_rep; 
      ADATXML::Array<int> mom; 
      int row; 
    }; 


    // hold three of them 
    struct subduce_rep_holder
    {
      std::map<std::string,subduce_rep> left,right,gamma; 
    };



    // loop over the full set of a cont_rep and get all possible subductions
    std::map<std::string,subduce_rep> 
      do_subduction( const std::map<std::string, cont_rep> &cont )
      {
        std::map<std::string,subduce_rep> ret; 
        std::map<std::string,cont_rep>::const_iterator it; 
        for(it = cont.begin(); it != cont.end(); ++it)
        {
          // split something like J1p,p000,r1 into 3 tokens on ","
          std::vector<std::string> tokens = tokenize(it->first,",");
          if( tokens.size() != 3 )
          {
            std::cout << __PRETTY_FUNCTION__ << ": token error" << std::endl;
            exit(1);
          }

          std::string cont_token = tokens[0]; 
          std::vector<const SubduceTableMap::irrep_sub_table*> sub_table; 
          std::vector<const SubduceTableMap::irrep_sub_table*>::const_iterator subduce_it; 
          ADATXML::Array<int> mom = it->second.mom; 
          sub_table = possible_subductions(cont_token,mom); 
        
          for(subduce_it = sub_table.begin(); subduce_it != sub_table.end(); ++subduce_it)
          {
            rHandle<CubicRep_p> lat_rep = (*subduce_it)->lat; 
            for(int row = 1; row <= lat_rep->dim(); ++row)
            {
              std::stringstream id; 
              id << lat_rep->rep_id() << ",p" << mom[0] << mom[1] << mom[2] << ",r" << row;
              std::string tag = id.str(); 
              ret[tag] = subduce_rep(cont_token,lat_rep->rep_id(),mom,row); 
              printer_function<subduce_print_1>(tag); 
            }

          } // close subduce loop 

        } // close cont map loop 

        return ret; 
      }


    subduce_rep_holder 
      get_subduce_reps(const rep_holder &r)
      {
        subduce_rep_holder ret; 
        ret.left = do_subduction(r.left);
        ret.right = do_subduction(r.right);
        ret.gamma = do_subduction(r.gamma); 
        return ret; 
      }


  } // anonomyous 


  // re subduce the continuum variants, this is hackey and unsatisfactory but 
  // the design parameters changed when we realized we also wanted to do 
  // unstable particles, can no longer use a helicity representation, must 
  // do the calculation in terms of cubic irreps. 
  std::vector<TaggedEnsemRedstarNPtBlock> 
    retag_subduced_lattice_xml( 
        const std::vector<TaggedEnsemRedstarNPtBlock> &cont_variant)
    {
      std::vector<TaggedEnsemRedstarNPtBlock> ret; 

      std::map<std::string,TaggedEnsemRedstarNPtBlock> cont_npts; 
      cont_npts = make_continuum_map(cont_variant); 

      rep_holder cont_reps = get_reps(cont_npts); 
      subduce_rep_holder sub_reps = get_subduce_reps(cont_reps); 


      // dont forget -- we need to do O(3) + 1 on the lorentz index
      // to put the photon into a helicity basis 
      //    -- can piggy back this off the one inside of ff_interface 

      return ret; 
    }


} // radmat 
