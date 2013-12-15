/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_rotated_levi_civita_tensor.cc

 * Purpose :

 * Creation Date : 13-12-2013

 * Last Modified : Fri 13 Dec 2013 02:27:07 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "lorentzff_rotated_levi_civita_tensor.h"
#include "lorentzff_canonical_rotations.h"
#include "radmat/utils/levi_civita.h"


namespace radmat
{

  namespace RotatedLeviCivitaTensorEnv
  {

    namespace
    {

      // this is far from optimal -- do it once then clone the result
      Tensor<double,4> rotate(const rHandle<RotationMatrix_t> R)
      {
        Tensor<double,4> levi = levi_civita<double,4>(); 

        return 
          contract(
              contract(
                contract(
                  contract(levi , *R, 3, 1 ),
                  *R , 2 , 1 ),
                *R , 1 , 1 ),
              *R , 0 , 1); 
      }

      // determine the lattice rotation  
      Tensor<double,4> gen_levi(const mom_t &l , const mom_t &r)
      {
        rHandle<RotationMatrix_t> R = radmat::LatticeRotationEnv::get_right_rotation(l,r);
        return rotate(R); 
      }

      // how we tag them -- use rotation group generator conventions 
      //  
      //  if we don't then Oh can mix up with D4 since we choose 
      //  to use D4 along the z-axis to orient Oh frames
      //
      std::string gen_id(const mom_t &l, const mom_t &r)
      {
        typedef radmat::LatticeRotationEnv::TheRotationGroupGenerator G; 
        return G::Instance().frame_label(l,r); 
      }

      void dump_keys(void)
      {
        std::map<std::string,Tensor<double,4> >::const_iterator it; 
        for (it = TheRotatedLeviCivitaTensorGenerator::Instance().begin();
            it != TheRotatedLeviCivitaTensorGenerator::Instance().end();
            ++it)
          std::cout << it->first << std::endl;
      }


      // reg function 
      bool do_reg(const mom_t &l, const mom_t &r)
      {
        typedef radmat::LatticeRotationEnv::TheRotationGroupGenerator G; 
        std::string id = gen_id(l,r);

        if( TheRotatedLeviCivitaTensorGenerator::Instance().find(id)
            != TheRotatedLeviCivitaTensorGenerator::Instance().end())
        {
          std::cout << __func__ << " double reg error for levi civita" 
           << "\nid: " << id <<  std::endl;
          dump_keys(); 
          throw std::string("double reg error TheRotatedLeviCivitaTensorEnv"); 
        }

        TheRotatedLeviCivitaTensorGenerator::Instance().insert(
            std::pair<std::string,Tensor<double,4> >(id,gen_levi(l,r))); 

        return true; 
      }

      // reg from a pair
      bool do_reg(const std::pair<mom_t,mom_t> &p)
      {
        return do_reg(p.first,p.second); 
      }

      // use the clone method to avoid continuously contracting 
      // over the indicies of a rank 4 tensor
      Tensor<double,4>* 
        query_factory(const mom_t &l, const mom_t &r)
        {
          std::string id = gen_id(l,r); 
          std::map<std::string,Tensor<double,4> >::const_iterator it; 
          it = TheRotatedLeviCivitaTensorGenerator::Instance().find(id);

          if( it == TheRotatedLeviCivitaTensorGenerator::Instance().end())
          {
            std::cout << __func__ << ": missing tag " << id << std::endl;
            exit(1); 
          }

          return it->second.clone(); 
        }

    } // anonomyous 

    // reg once 
    namespace
    {
      bool local_registration = false; 
    }

    // reg all frames in the frame map
    bool registerAll(void)
    {
      bool success = true; 
      if(!!! local_registration)
      {
        typedef radmat::LatticeRotationEnv::TheRotationGroupGenerator G; 
        std::vector<std::string> u_frames; 
        std::vector<std::string>::const_iterator u_it; 
        u_frames = G::Instance().unique_frames(); 

        for(u_it = u_frames.begin(); u_it != u_frames.end(); ++u_it)
        {
          std::vector<std::string> r_frames; 
          std::vector<std::string>::const_iterator r_it; 
          std::pair<mom_t,mom_t> can = G::Instance().get_frame_momentum(*u_it); 
          r_frames = G::Instance().get_related_frames(can.first,can.second); 

          for(r_it = r_frames.begin(); r_it != r_frames.end(); ++r_it)
            success &= do_reg( G::Instance().get_frame_momentum(*r_it) ); 
        }

        local_registration = true;
      }
      return success; 
    }


    // only interface
    rHandle<Tensor<double,4> >
      call_factory(const mom_t &l, const mom_t &r)
      {
        return rHandle<Tensor<double,4> >(query_factory(l,r)); 
      }

    // piggy back 
    rHandle<Tensor<double,4> >
      call_factory(const Tensor<double,1> &l , const Tensor<double,1> &r, const double kick)
      {
        mom_t ll = get_space_mom(l,kick);
        mom_t rr = get_space_mom(r,kick); 

        return call_factory(ll,rr); 
      }

  } // RotatedLeviCivitaTensorEnv

} // radmat



