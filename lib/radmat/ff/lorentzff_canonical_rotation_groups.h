#ifndef LORENTZFF_CANONICAL_ROTATION_GROUPS_H
#define LORENTZFF_CANONICAL_ROTATION_GROUPS_H 






/*
 *
 *  Roughly this guy picks out canonical frames
 *  for evaluating matrix elements in
 *
 *  when we go through and evaluate helicity 
 *  matrix elements the manner in which we 
 *  build them attaches an extra phase
 *  corresponding to a rotation about the
 *  momentum direction of one of the particles
 *
 *  this class picks out which frames we want 
 *  and the canonical rotation classes then
 *  work together to apply a convention 
 *  ( an extra rotation ) to kill the phase
 *
 *  the choice of convention is arbitrary 
 *  the phase only appears when both paricles 
 *  are spin-1 or higher 
 *
 */



#include "lorentzff_canonical_rotations_utils.h"
#include "radmat/utils/stringify.h"
#include "io/adat_xmlio.h"
#include "hadron/irrep_utils.h"
#include <exception>
#include <sstream>
#include <string>
#include <map>
#include <vector>

namespace radmat
{

  struct Rep_p
  {
    std::string id(void) const = 0; 
  };

  template<typename T> 
    struct Rep 
    : public Rep_p 
    {
      std::string id(void) {return Stringify<T>();} 
    };

  struct Oh : public Rep<Oh> {};
  struct D2 : public Rep<D2> {};
  struct D3 : public Rep<D3> {};
  struct D4 : public Rep<D4> {};

  REGISTER_STRINGIFY_TYPE(Oh);
  REGISTER_STRINGIFY_TYPE(D2);
  REGISTER_STRINGIFY_TYPE(D3);
  REGISTER_STRINGIFY_TYPE(D4);

  template<int MomMax, int InsMax>
    struct RotationGroupGenerator
    {
      struct RepPair
      {
        RepPair(const mom_t &left, const mom_t &right)
          : l(left) , r(right), lefty(gen_rep(l)) , righty(gen_rep(r))
        { }

        mom_t l;
        mom_t r; 
        rHandle<Rep> lefty;
        rHandle<Rep> righty; 
      }; 

      // utility 
      rHandle<Rep_p> gen_rep(const mom_t &p)
      {
        std::string LG = Hadron::generateLittleGroup(p); 
        if(LG == "Oh")
          return rHandle<Rep_p>(new Oh() ); 
        else if (LG == "D2")
          return rHandle<Rep_p>(new D2() ); 
        else if (LG == "D3")
          return rHandle<Rep_p>(new D3() ); 
        else if (LG == "D4")
          return rHandle<Rep_p>(new D4() ); 
        else
        {
          std::cout << "LG " << LG << " not supported" << std::endl;
          exit(1); 
        }
        exit(1); 
      }

      mom_t rest_specialization(void)
      {
        return gen_mom<0,0,1>(); 
      }

      std::string label(const mom_t &p, const rHandle<Rep> &r)
      {
        std::stringstream ss; 
        mom_t lab = p;
        if(r->id() == Stringify<Oh>())
          lab = rest_specialization(); 

        ss << r->id() << "_" << lab[0] << lab[1] << lab[2];
        return ss.str(); 
      }


      // specialize rest to have orientation z-axis
      std::pair<mom_t,mom_t> frame(const RepPair &rep)
      {
        std::pair<mom_t,mom_t> r;
        r.first = rep.l;
        r.second = rep.r; 
        if( rep.lefty->id() == Stringify<Oh>())
          r.first = rest_specialization();
        if( rep.righty->id() == Stringify<Oh>())
          r.second = rest_specialization();
        return r; 
      }

      std::pair<mom_t,mom_t> momentum(const RepPair &rep)
      {
        std::pair<mom_t,mom_t> r;
        r.first = rep.l;
        r.second = rep.r; 
        return r; 
      }

      std::string frame_label(const RepPair &rep)
      {
        std::stringstream ss; 
        ss << "lefty_" << label(rep.l,rep.lefty);
        ss << ".righty_" << label(rep.r,rep.righty);
        return ss.str(); 
      }

      bool registerAll(void)
      {
        initialize_frames(); 
        return true; 
      }

      mom_t mmom(const int a, const int b, const int c) const
      {
        mom_t foo;
        foo.resize(3);
        foo[0] = a;
        foo[1] = b;
        foo[2] = c; 
        return foo; 
      }

      int find_frame(const RepPair &r)
      {
        int pos(-1);
        std::pair<mom_t,mom_t> f = frame(r); 

        for(int i = 0; i < frames.size(); ++i)
        {
          std::pair<mom_t,mom_t> c = frame(frames[i]); 
          if(related_by_rotation(c.first,c.second,f.first,f.second,false)i)
            pos = i; 
        }

        return i; 
      }

      void insert_can_frame(const RepPair &r)
      {
        frames.push_back(r); 
        can_frame_map.insert(std::pair<std::string,std::string>(frame_label(r),frame_label(r))); 
      }

      void insert_related_frame(const int pos, const RepPair &r)
      {
        std::string id = frame_label(r);
        std::string can_id = frame_label(frames[pos]); 
        if(can_frame_map.find(id) != can_frame_map.end())
          throw std::string("frames double reg error"); 
        can_frame_map.insert(std::pair<std::string,std::string>(id,can_id)); 
      }

      void register_frame(const RepPair &r)
      {
        std::string id = frame_label(r);
        if(frame_label_map.find(id) != frame_label_map.end())
          throw std::string("frames double reg error"); 
        frame_label_map.insert(std::pair<std::string,RepPair>(id,r)); 
      }


      void insert(const RepPair &r)
      {
        int pos = find_frame(r); 

        if( pos == -1 ) 
          insert_can_frame(r); 
        else
          insert_related_frame(pos,r); 

        register_frame(r); 
      }


      void initialize_frames(void)
      {
        // the loop ordering here sets the canonical names, 
        // use positive then lower so it is human readable
        for(int i = MOM_MAX; i >= -MOM_MAX; --i)
          for(int j = MOM_MAX; j >= -MOM_MAX; --j)
            for(int k = MOM_MAX; k >= -MOM_MAX; --k)
              for(int l = MOM_MAX; l >= -MOM_MAX; --l)
                for(int m = MOM_MAX; m >= -MOM_MAX; --m)
                  for(int n = MOM_MAX; n >= -MOM_MAX; --n)
                  {
                    if( (i-l)*(i-l) + (j-m)*(j-m) + (k-n)*(k-n) > INSERTION_MAX_SQ )
                      continue;
                    if( i*i + j*j + k*k > MOM_MAX * MOM_MAX )
                      continue; 
                    if( l*l + m*m + n*n > MOM_MAX * MOM_MAX )
                      continue; 

                    insert( RepPair( mmom(i,j,k),mmom(l,m,n) ) ); 
                  }
      }
  
      template<typename T> 
        void do_exit(const std::string &s, const T &t)
        {
          std::cout << __PRETTY_FUNCTION__ << s << t << std::endl;
          exit(1); 
        }

      // Query
      /////////// 

      std::string get_can_frame_string(const mom_t &l, const mom_t &r)
      {
        std::string lab = frame_label(RepPair(l,r)); 
        std::map<std::string,std::string>::const_iterator it = can_frame_map.find(lab); 
        if(it == can_frame_map.end())
          do_exit("missing label" , lab);
        return it->second; 
      }

      std::pair<mom_t,mom_t> get_can_frame(const mom_t &l, const mom_t &r) const
      {
        std::string can_frame = get_can_frame_string(l,r); 
        return get_frame_momentum(can_frame); 
      }

      std::pair<mom_t,mom_t> get_frame_momentum(const std::string &id) const
      {
        std::map<std::string,RepPair>::const_iterator it; 
        it = frame_label_map.find(id); 
        if(it == frame_label_map.end())
          do_exit("missing label " , id);
        return momentum(it->second); 
      }

      std::vector<std::string> RotationGroupGenerator_untemp::get_related_frames(const mom_t &l, const mom_t &r) const
      {
        std::string can = get_can_frame_string(l,r); 
        std::vector<std::string> ret; 
        std::map<std::string,std::string>::const_iterator it; 
        for(it = can_frame_map.begin(); it != can_frame_map.end(); ++it)
          if( it->second == can )
            ret.push_back(it->first); 

        return ret; 
      }

      std::vector<std::string> RotationGroupGenerator_untemp::unique_frames(void) const
      {
        std::vector<RepPair>::const_iterator it; 
        std::vector<std::string> ret; 
        for( it =  frames.begin(); it != frames.end(); ++it)
          ret.push_back(frame_label(*it)); 

        return ret; 
      }

      // key is frame label value is the momentum pair
      std::map<std::string,RepPair> frame_label_map;
      // key is frame label value is canonical frame label
      std::map<std::string,std::string> can_frame_map; 
      // a record of canonical frames 
      std::vector<RepPair> frames; 
    };


} // radmat 



#endif /* LORENTZFF_CANONICAL_ROTATION_GROUPS_H */
