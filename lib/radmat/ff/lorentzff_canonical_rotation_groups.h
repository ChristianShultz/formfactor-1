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
 *  NB: this picks out an orientation, in the case
 *  of a particle at rest we use the convention 
 *  in adat that helicity and J_z states overlap
 *  (polarized along the z-axis) - this is done 
 *  by pretending a particle at rest is in the 
 *  LG D4 w/ one unit of momentum along the z-axis 
 *  
 */



#include "lorentzff_canonical_rotations_utils.h"
#include "lorentzff_cubic_reps.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/stringify.h"
#include "io/adat_xmlio.h"
#include "hadron/irrep_util.h"
#include <exception>
#include <sstream>
#include <string>
#include <map>
#include <vector>

namespace radmat
{
  template<int MOM_MAX, int INSERTION_MAX_SQ>
    struct RotationGroupGenerator
    {

      private:

      // how do we treat rest internally
      mom_t rest_specialization(void) const
      {
        return gen_mom<0,0,1>(); 
      }

      // what do we check for before we talk to the world
      mom_t invert_rest_specialization(const mom_t &p) const
      {
        mom_t f = rest_specialization(); 
        if ((p[0] != f[0]) || (p[1] != f[1]) || (p[2] != f[2]) )
        {
          std::cout << __func__ << ": rot error - "
           << p[0] << p[1] << p[2] << std::endl;
          throw std::string("rest rotation error"); 
          exit(1); 
        }
        return gen_mom<0,0,0>(); 
      }

      // half of a map label 
      std::string label(const mom_t &p, const rHandle<Rep_p> &r) const
      {
        std::stringstream ss; 
        mom_t lab = p;
        if(r->rep_id() == Stringify<Oh>())
          lab = rest_specialization(); 

        ss << r->rep_id() << "_" << lab[0] << lab[1] << lab[2];
        return ss.str(); 
      }


      // specialize rest to have orientation z-axis
      std::pair<mom_t,mom_t> frame(const RepPair &rep) const
      {
        std::pair<mom_t,mom_t> r;
        r.first = rep.l;
        r.second = rep.r; 
        if( rep.lefty->rep_id() == Stringify<Oh>())
          r.first = rest_specialization();
        if( rep.righty->rep_id() == Stringify<Oh>())
          r.second = rest_specialization();
        return r; 
      }

      // undo specialize rest to have orientation z-axis
      std::pair<mom_t,mom_t> momentum(const RepPair &rep) const
      {
        std::pair<mom_t,mom_t> r;
        r.first = rep.l;
        r.second = rep.r; 
        if( rep.lefty->rep_id() == Stringify<Oh>())
          r.first = invert_rest_specialization(rep.l);
        if( rep.righty->rep_id() == Stringify<Oh>())
          r.second = invert_rest_specialization(rep.r); 
        return r; 
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

      bool same_reps(const RepPair &r, const RepPair &rr) const
      {
        return ( (r.lefty->rep_id() == rr.lefty->rep_id()) 
            && (r.righty->rep_id() == rr.righty->rep_id())); 
      }

      // frames in which one particle is at rest have the 
      // z axis chosen as a canonical direction 
      int find_frame(const RepPair &r) const
      {
        int pos(-1);

        std::pair<mom_t,mom_t> f = frame(r); 

        for(int i = 0; i < frames.size(); ++i)
        {
          if(!!! same_reps(frames[i],r) )
            continue; 

          std::pair<mom_t,mom_t> c = frame(frames[i]); 
          if(related_by_rotation(c.first,c.second,f.first,f.second,false))
            pos = i; 
        }
        return pos; 
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


      // inserted rest frames are polarized along the z-axis
      void insert(const RepPair &actual)
      {
        RepPair sort(actual.l,actual.r); 
        std::pair<mom_t,mom_t> p = frame(actual); 

        sort.l = p.first;
        sort.r = p.second; 

        int pos = find_frame(sort); 

        if( pos == -1 ) 
          insert_can_frame(sort); 
        else
          insert_related_frame(pos,sort); 

        register_frame(sort); 
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
        void do_exit(const std::string &s, const T &t) const
        {
          std::cout << __PRETTY_FUNCTION__ << s << t << std::endl;
          exit(1); 
        }

      // public access
      public:

      bool registerAll(void)
      {
        initialize_frames(); 
        return true; 
      }

      // Query
      /////////// 

      std::string frame_label(const RepPair &rep) const
      {
        std::stringstream ss; 
        ss << "lefty_" << label(rep.l,rep.lefty);
        ss << ".righty_" << label(rep.r,rep.righty);
        return ss.str(); 
      }

      std::string frame_label(const mom_t &l, const mom_t &r) const
      {
        return frame_label(RepPair(l,r)); 
      }

      std::string get_can_frame_string(const mom_t &l, const mom_t &r) const
      {
        std::string lab = frame_label(RepPair(l,r)); 
        std::map<std::string,std::string>::const_iterator it = can_frame_map.find(lab); 
        if(it == can_frame_map.end())
          do_exit("missing label" , lab);

        return it->second; 
      }

      std::pair<mom_t,mom_t> get_can_frame_orientation(const mom_t &l, const mom_t &r) const
      {
        std::string can_frame = get_can_frame_string(l,r); 
        return get_frame_orientation(can_frame); 
      }

      std::pair<mom_t,mom_t> get_frame_orientation(const std::string &id) const
      {
        typename std::map<std::string,RepPair>::const_iterator it; 
        it = frame_label_map.find(id); 
        if(it == frame_label_map.end())
          do_exit("missing label " , id);

        return frame(it->second); 
      }

      std::pair<mom_t,mom_t> get_frame_orientation(const mom_t&l, const mom_t &r) const
      {
        return frame( RepPair(l,r) ); 
      }

      // in the case of rest this is not the same as the orientation!!!
      std::pair<mom_t,mom_t> get_frame_momentum(const std::string &id) const
      {
        typename std::map<std::string,RepPair>::const_iterator it; 
        it = frame_label_map.find(id); 
        if(it == frame_label_map.end())
          do_exit("missing label " , id);

        return momentum(it->second); 
      }

      std::vector<std::string> get_related_frames(const mom_t &l, const mom_t &r) const
      {
        std::string can = get_can_frame_string(l,r); 
        std::vector<std::string> ret; 
        std::map<std::string,std::string>::const_iterator it; 
        for(it = can_frame_map.begin(); it != can_frame_map.end(); ++it)
          if( it->second == can )
            ret.push_back(it->first); 

        return ret; 
      }

      std::vector<std::string> unique_frames(void) const
      {
        typename std::vector<RepPair>::const_iterator it; 
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
