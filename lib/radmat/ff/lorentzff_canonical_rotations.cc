/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_canonical_rotations.cc

 * Purpose :


 * Last Modified : Wed 18 Dec 2013 11:37:25 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "lorentzff_canonical_rotations.h"
#include "lorentzff_canonical_rotations_utils.h"
#include "formfac/formfac_qsq.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/pow2assert.h"
#include "itpp/itbase.h"
#include <sstream>
#include <exception>
#include <omp.h>
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"



namespace radmat
{


  namespace LatticeRotationEnv
  {

    namespace
    {
      void check_exit_LG(const mom_t &p, const std::string &LG)
      {
        if( Hadron::generateLittleGroup(p) != LG )
        {
          std::cout << __func__ << ": error, expected LG " 
            << LG << ", got momentum " << p[0] 
            << p[1] << p[2] << "  -->  " << Hadron::generateLittleGroup(p)
            << ": exiting" << std::endl;
          exit(1); 
        }
      }


      std::string mom_string(const mom_t &p)
      {
        std::stringstream ss; 
        ss << p[0] << p[1] << p[2]; 
        return ss.str(); 
      }

      bool
        check_rotation_invariant(const RotationMatrix_t *R,
            const mom_t &l)
        {
          mom_t chk = gen_mom<0,0,0>();

          for(int i = 0; i < 3; ++i)
          {
            double res(0.); 
            for(int j = 0; j < 3; ++j)
              if ( fabs( (*R)[i+1][j+1] ) > 1e-6 )
                res += (*R)[i+1][j+1] * double(l[j]); 

            chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
          }

          if( (chk[0] != l[0]) || (chk[1] != l[1]) || (chk[2] != l[2]) )
          {
            std::cout << __func__ << ": Rotation error," 
              << " in " << l[0] << l[1] << l[2] 
              << " out " << chk[0] << chk[1] << chk[2]
              << "\n R:\n" << *R << std::endl;

            return false; 
          }

          return true; 
        } 



      //  Apply the transformation to 
      //  particle at rest to maintain a
      //  consistent frame
      RotationMatrix_t *
        rest_specialization(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &il, 
            const mom_t &ir) 
        {
          // this kind of sucks..
          mom_t l = TheRotationGroupGenerator::Instance().get_frame_orientation(il,ir).first;
          mom_t r = TheRotationGroupGenerator::Instance().get_frame_orientation(il,ir).second;

          std::string LGl = Hadron::generateLittleGroup(cl);
          std::string LGr = Hadron::generateLittleGroup(cr);

          check_exit_LG(l,LGl);
          check_exit_LG(r,LGr); 

          // if both are at rest return identity
          if( is_rest(cr) && is_rest(cl) )
          {
            RotationMatrix_t *I4= new Tensor<double,2>((TensorShape<2>())[4][4],0.);
            I4->lower_index(1); 
            (*I4)[0][0] = 1.;
            (*I4)[1][1] = 1.;
            (*I4)[2][2] = 1.;
            (*I4)[3][3] = 1.;
          }

          RotationMatrix_t *R; 
          if(is_rest(ir))
            R = generate_frame_transformation(l,cl); 
          else if(is_rest(il))
            R = generate_frame_transformation(r,cr); 
          else
          {
            throw std::string("rotation error"); 
            exit(1); 
          }

          R->lower_index(1); 
          clean_up_rot_mat(R); 

          return R; 
        }

      // generates the frame transformation between cr and r
      // 
      //  cr = R r
      //
      RotationMatrix_t *
        generate_right_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          std::string LGl = Hadron::generateLittleGroup(cl);
          std::string LGr = Hadron::generateLittleGroup(cr);

          check_exit_LG(l,LGl);
          check_exit_LG(r,LGr); 

          if( is_rest(cl) || is_rest(cr) )
            return rest_specialization(cl,cr,l,r); 

          RotationMatrix_t *R = generate_frame_transformation(r,cr); 
          R->lower_index(1); 
          clean_up_rot_mat(R); 

          if( !!! check_total_frame_transformation(R,l,r,cl,cr,true) )
          {
            std::cout << __func__ << ": error rotating frame " 
              << std::endl;
            throw std::string("right rotation error"); 
          }

          return R; 
        }


      // attaches the phase to the transformation between 
      // cl and l that respects cr to r
      RotationMatrix_t *
        generate_left_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          std::string LGl = Hadron::generateLittleGroup(cl);
          std::string LGr = Hadron::generateLittleGroup(cr);

          check_exit_LG(l,LGl);
          check_exit_LG(r,LGr); 

          if( is_rest(cl) || is_rest(cr) )
            return rest_specialization(cl,cr,l,r); 


          // the rotation from cl to l
          RotationMatrix_t *R = generate_frame_transformation(l,cl); 
          R->lower_index(1); 
          clean_up_rot_mat(R); 

          if( !!! check_total_frame_transformation(R,l,r,cl,cr,true) )
          {
            std::cout << __func__ << ": error rotating frame " 
              << std::endl;
            throw std::string("right rotation error"); 
          }

          return R; 
        }


      RotationMatrix_t * 
        primitive_rotation(const mom_t &cl,
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r,
            const std::string &id)
        {
          if( is_rest(l) || is_rest(r) )
            return rest_specialization(cl,cr,l,r); 
          else if( id == "lefty")
            return generate_left_rotation(cl,cr,l,r); 
          else if( id == "righty")
            return generate_right_rotation(cl,cr,l,r); 
          else
          {
            std::cout << __func__ << ": error " << std::endl;
            exit(1); 
          }
          exit(1); 
        }


    } // anonomyous




    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    namespace
    {
      bool local_registration = false; 

      RotationMatrix_t *
        left_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          return generate_left_rotation(cl,cr,l,r); 
        }

      RotationMatrix_t *
        right_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          return generate_right_rotation(cl,cr,l,r); 
        }

    } // anonomyous
    



    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////


    bool registerAll(void)
    {
      bool success = true; 

      if ( !!! local_registration)
      {
        success &= TheRotationGroupGenerator::Instance().registerAll(); 
        local_registration = true; 
      }

      if( !!! success )
      {
        throw std::string("reg error in LatticeRotationEnv");
      }

      return success; 
    }

    //////////////////////////////////////////////////////////////////////

    std::string
      rotation_group_label(const mom_t &l, const mom_t &r)
      {
        return TheRotationGroupGenerator::Instance().get_can_frame_string(l,r);  
      }
  
    //////////////////////////////////////////////////////////////////////

    FrameOrientation_t 
      get_frame_orientation(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> can = TheRotationGroupGenerator::Instance().get_can_frame_orientation(l,r); 
        std::pair<mom_t,mom_t> f = TheRotationGroupGenerator::Instance().get_frame_orientation(l,r); 

        return FrameOrientation_t(can.first,can.second,f.first,f.second); 
      }

    //////////////////////////////////////////////////////////////////////

    rHandle<RotationMatrix_t> 
      get_left_rotation(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> can; 
        can = TheRotationGroupGenerator::Instance().get_can_frame_orientation(l,r); 
        return rHandle<RotationMatrix_t>(primitive_rotation(can.first,can.second,l,r,"lefty")); 
      }

    //////////////////////////////////////////////////////////////////////

    rHandle<RotationMatrix_t> 
      get_right_rotation(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> can; 
        can = TheRotationGroupGenerator::Instance().get_can_frame_orientation(l,r); 
        return rHandle<RotationMatrix_t>(primitive_rotation(can.first,can.second,l,r,"righty")); 
      }

    //////////////////////////////////////////////////////////////////////
    
    rHandle<RotationMatrix_t> 
      get_left_can_frame_rotation(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> can; 
        can = TheRotationGroupGenerator::Instance().get_can_frame_orientation(l,r); 
        return rHandle<RotationMatrix_t>(radmat::CanonicalRotationEnv::call_factory(can.first)); 
      }

    //////////////////////////////////////////////////////////////////////
    
    rHandle<RotationMatrix_t> 
      get_right_can_frame_rotation(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> can; 
        can = TheRotationGroupGenerator::Instance().get_can_frame_orientation(l,r); 
        return rHandle<RotationMatrix_t>(radmat::CanonicalRotationEnv::call_factory(can.second)); 
      }

    //////////////////////////////////////////////////////////////////////

    RotationMatrix_t*
      get_left_rotation(const mom_t &cl, 
          const mom_t &cr, 
          const mom_t &l, 
          const mom_t &r) 
      {
        return left_rotation(cl,cr,l,r); 
      }

    //////////////////////////////////////////////////////////////////////

    RotationMatrix_t*
      get_right_rotation(const mom_t &cl, 
          const mom_t &cr, 
          const mom_t &l, 
          const mom_t &r) 
      {
        return right_rotation(cl,cr,l,r); 
      }


  } // LatticeRotationEnv

} // radmat




