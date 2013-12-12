/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_canonical_rotations_utils.cc

 * Purpose :

 * Creation Date : 10-12-2013

 * Last Modified : Thu 12 Dec 2013 09:10:41 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "lorentzff_canonical_rotations_utils.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include <sstream>

namespace radmat
{

  // Ax = b
  bool
    check_frame_transformation(const RotationMatrix_t *A,
        const mom_t &x, 
        const mom_t &b)
    {
      mom_t chk = gen_mom<0,0,0>();

      for(int i = 0; i < 3; ++i)
      {
        double res(0.); 
        for(int j = 0; j < 3; ++j)
          if ( fabs( (*A)[i+1][j+1] ) > 1e-6 )
            res += (*A)[i+1][j+1] * double(x[j]); 

        chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
      }

      if( (chk[0] != b[0]) || (chk[1] != b[1]) || (chk[2] != b[2]) )
        return false; 

      return true; 
    } 

  namespace
  {
    std::string string_mom(const mom_t &p)
    {
      std::stringstream ss;
      ss << p[0] << p[1] << p[2];
      return ss.str(); 
    }
  }


  RotationMatrix_t*
    generate_frame_transformation(const mom_t &p, const mom_t &pp)
    {
      if( dot_mom_t(p,p) != dot_mom_t(pp,pp) )
      {
        std::cout << __func__ << ": error, frames are unrelated " 
          << "\np" << string_mom(p) << " pp " << string_mom(pp) 
          << std::endl;
        exit(1); 
      }
     
      RotationMatrix_t *R = new Tensor<double,2>( (TensorShape<2>())[4][4], 0.); 
      RotationMatrix_t *Rp = radmat::CanonicalLatticeRotationEnv::call_factory(p);
      RotationMatrix_t *Rpp = radmat::CanonicalLatticeRotationEnv::call_factory(pp);
      R->lower_index(1); 

      for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
        {
          if( fabs( (*Rp)[i][j] ) < 1e-6 )
            continue;

          for(int k = 0; k < 4; ++k)
          (*R)[i][k] += (*Rp)[i][j] * (*Rpp)[k][j]; 
        }

      delete Rp;
      delete Rpp; 

      if( !!! check_frame_transformation(R,pp,p) )
      {
        std::cout << __func__ << ": error transforming frames "
          << "\np" << string_mom(p) << " pp " << string_mom(pp) 
          << " p = Rpp\nR:" << *R << std::endl;
        delete R; 
        exit(1); 
      }

      return R; 
    }


  int dot_mom_t(const mom_t &left, const mom_t &right)
  {
    return left[0]*right[0] + left[1]*right[1] + left[2]*right[2]; 
  }

  mom_t cross_product(const mom_t &l, const mom_t &r)
  {
    mom_t c = gen_mom<0,0,0>(); 
    c[0] = l[1]*r[2] - l[2]*r[1];
    c[1] = -(l[0]*r[2] - l[2]*r[0]);
    c[2] = l[0]*r[1] - l[1]*r[0]; 
    return c; 
  } 

  bool is_rest(const mom_t &p)
  {
    return (dot_mom_t(p,p) == 0);
  }

  namespace
  {
    bool same_length(const mom_t &l, const mom_t &ll)
    {
      return ( dot_mom_t(l,l) == dot_mom_t(ll,ll) ); 
    }

    double cos_theta(const mom_t &l, const mom_t &ll)
    {
      if(is_rest(l) || is_rest(ll) )
        return 0.;

      return double(dot_mom_t(l,ll))/sqrt(double(dot_mom_t(l,l) * dot_mom_t(ll,ll))); 
    }

    int volume_element(const mom_t &l, const mom_t &r)
    {
      return ( dot_mom_t(l,cross_product(l,r)) ); 
    }

  }

  bool related_by_Oh_rotation(const mom_t &left, const mom_t &lleft)
  {
    return (
        (Hadron::generateLittleGroup(left) == Hadron::generateLittleGroup(lleft)) 
        && (dot_mom_t(left,left) == dot_mom_t(lleft,lleft))
        ); 
  }


  // test the dot product and volume element are the same
  bool related_by_rotation(const mom_t &left, const mom_t &right, 
      const mom_t &lleft, const mom_t &rright,
      const bool allow_flip)
  {
    // knock out left or right is at rest first 
    if( is_rest(left) && is_rest(lleft) )
      return related_by_Oh_rotation(right,rright); 

    if( is_rest(right) && is_rest(rright) )
      return related_by_Oh_rotation(left,lleft); 


    // easy checks first
    bool success = true; 
    success &= (cos_theta(left,right) == cos_theta(lleft,rright));
    success &= (volume_element(left,right) == volume_element(lleft,rright));

    if( allow_flip ) 
    {
      if(dot_mom_t(left,left) == dot_mom_t(lleft,lleft))
        success &= (dot_mom_t(right,right) == dot_mom_t(rright,rright));  
      else if(dot_mom_t(left,left) == dot_mom_t(rright,rright)) 
        success &= (dot_mom_t(right,right) == dot_mom_t(lleft,lleft));  
      else
        return false; 
    }
    else
    {
      success &= (dot_mom_t(left,left) == dot_mom_t(lleft,lleft)); 
      success &= (dot_mom_t(right,right) == dot_mom_t(rright,rright));  
    }

    if( !!! success )
      return false; 


    // this is a bit clunky 
    RotationMatrix_t *Rr = radmat::CanonicalLatticeRotationEnv::call_factory(right);
    RotationMatrix_t *Rrr = radmat::CanonicalLatticeRotationEnv::call_factory(rright);
    RotationMatrix_t *R = new Tensor<double,2>((TensorShape<2>())[4][4],0.); 

    //  Christopher was nice enough to program all of the euler 
    //  angles, lets make use of them and test rather than search
    //  for a possible rotation relating the two frames 
    //
    //  this rotation takes me from the frame defined by rright
    //  to pref then from pref to right
    // 
    //  right = Rr * Rrr^T rright
    // 
    //  if they are related by a rotaion then 
    //  we should get left from the same op on lleft
    //
    for(int i = 1; i < 4; ++i)
      for(int j = 1; j < 4; ++j)
        for(int k = 1; k < 4; ++k)
          (*R)[i][k] += (*Rr)[i][j] * (*Rrr)[k][j] ; 

    success &= check_frame_transformation(R,lleft,left); 

    delete R; 
    delete Rr; 
    delete Rrr; 

    return success;
  }

  itpp::Vec<double> normalize(const mom_t &m)
  {
    itpp::Vec<double> r(3); 
    double norm = double(dot_mom_t(m,m)); 
    r[0] = double(m[0])/norm;
    r[1] = double(m[1])/norm;
    r[2] = double(m[2])/norm;
    return r; 
  }


}
