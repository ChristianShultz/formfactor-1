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
#include "io/adat_xmlio.h"
#include <exception>
#include <sstream>
#include <string>
#include <map>
#include <vector>

namespace radmat
{



  template<int MOM_MAX, int INSERTION_MAX_SQ> 
    struct
    RotationGroupGenerator
    {
      typedef ADATXML::Array<int> mom_t; 

      RotationGroupGenerator()
      {}

      bool BlastOff(void)
      {
        std::cout << __PRETTY_FUNCTION__ << ": initializing" << std::endl;
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

      std::string get_can_frame_string(const mom_t &l, const mom_t &r) const
      {
        std::string id = reg_id(l,r); 
        std::map<std::string,std::string>::const_iterator it; 
        it = can_frame_map.find(id); 
        if(it == can_frame_map.end())
        {
          std::cout << "error: id " << id << "not present" << std::endl;
          exit(1); 
        }
        return it->second; 
      } 

      std::pair<mom_t,mom_t> get_can_frame(const mom_t &l, const mom_t &r) const
      {
        std::string frame = get_can_frame_string(l,r); 
        std::map<std::string,std::pair<mom_t,mom_t> >::const_iterator it; 
        it = frame_label_map.find(frame); 
        if(it == frame_label_map.end())
        {
          std::cout << "error: frame " << frame << "not present" << std::endl;
          exit(1); 
        }
        return it->second; 
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
        std::vector<std::pair<mom_t,mom_t> >::const_iterator it; 
        std::vector<std::string> ret; 
        for( it =  frames.begin(); it != frames.end(); ++it)
          ret.push_back(reg_id(it->first,it->second)); 

        return ret; 
      }

      private:

      void insert(const int a, const int b, const int c, 
          const int aa, const int bb, const int cc)
      {
        insert( mmom(a,b,c) , mmom(aa,bb,cc) ); 
      }

      std::string reg_id(const mom_t &l, const mom_t &r) const
      {
        std::stringstream ss; 
        ss << "lefty_" << Hadron::generateLittleGroup(l) << "_p" << l[0] << l[1] << l[2];
        ss << ".righty_" << Hadron::generateLittleGroup(r) << "_p" << r[0] << r[1] << r[2]; 

        return ss.str(); 
      }

      void insert(const mom_t &l, const mom_t &r)
      {
        int pos(-1); 
        for(int i = 0; i < frames.size(); ++i)
          if( related_by_rotation(frames[i].first,frames[i].second,l,r,false) )
            pos = i; 

        if( pos == -1 ) 
        {
          frames.push_back(std::pair<mom_t,mom_t>(l,r)); 
          can_frame_map.insert(std::pair<std::string,std::string>(
                reg_id(l,r),reg_id(l,r) )); 
        }
        else
        {
          std::string lab = reg_id(l,r); 
          std::string can_lab = reg_id(frames[pos].first,frames[pos].second); 
          if( can_frame_map.find(lab) != can_frame_map.end() )
          {
            std::cout << __func__ << ": double insert error" 
              << std::endl;
            throw std::string("canonical frame double insert error"); 
          }

          can_frame_map.insert(std::pair<std::string,std::string>(lab,can_lab)); 
        }

        frame_label_map.insert(std::pair<std::string,std::pair<mom_t,mom_t> >(
              reg_id(l,r) , std::pair<mom_t,mom_t>(l,r))); 
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

                    insert(i,j,k,l,m,n);  
                  }
      }

      // unique frames
      std::vector<std::pair<mom_t,mom_t> > frames; 
      // key is rotated frame id, value is canonical id
      std::map<std::string,std::string> can_frame_map; 
      // key is frame id, value is momentum pair
      std::map<std::string,std::pair<mom_t,mom_t> > frame_label_map;
    }; 



} // radmat 


#endif /* LORENTZFF_CANONICAL_ROTATION_GROUPS_H */
