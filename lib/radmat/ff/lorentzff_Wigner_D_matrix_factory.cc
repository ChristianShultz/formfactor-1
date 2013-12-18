/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_Wigner_D_matrix_factory.cc

 * Purpose :

 * Creation Date : 14-12-2013

 * Last Modified : Wed 18 Dec 2013 03:06:37 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "lorentzff_Wigner_D_matrix_factory.h"
#include "lorentzff_canonical_rotations_utils.h"
#include "lorentzff_formfac_utils.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "semble/semble_meta.h"
#include <sstream>
#include <exception> 


namespace radmat
{

  namespace
  {
    
    WignerMatrix_t gen_wigner_matrix(const mom_t &p, const int J)
    {
      int bound = 2*J + 1; 
      WignerMatrix_t W( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.) ); 
      Hadron::CubicCanonicalRotation_t eul = Hadron::cubicCanonicalRotation(p); 

      for(int m1 = -J; m1 <= J; ++m1)
        for(int m2 = -J; m2 <= J; ++m2)
        {
          std::complex<double> cd =  SEMBLE::toScalar( 
              Hadron::Wigner_D( 2*J , 2*m1 , 2*m2 , eul.alpha, eul.beta, eul.gamma) );  
            W[J-m1][J-m2] = round_to_zero(cd,1e-6); 
        }

      return W; 
    }

    std::string gen_id(const mom_t &p, const int J)
    {
      std::stringstream ss; 
      ss << "J_" << J << "__" << Hadron::generateLittleGroup(p) 
        << "_p"<< p[0] << p[1] << p[2] << std::endl;  
      return ss.str(); 
    }

    typedef radmat::WignerDMatrixEnv::TheWignerDMatrixFactory WDFac; 

    bool reg_wigner_matrix(const WignerMatrix_t &W, const std::string &id)
    {
      if( WDFac::Instance().find(id) != WDFac::Instance().end() )
        throw std::string("WignerDMatrixEnv double reg error"); 

      WDFac::Instance().insert(std::pair<std::string,WignerMatrix_t>(id,W)); 
      return true; 
    }

    WignerMatrix_t* query_factory(const mom_t &p, const int J)
    {
      std::string id = gen_id(p,J); 
      std::map<std::string,WignerMatrix_t>::const_iterator it; 
      it = WDFac::Instance().find(id); 
      if( it == WDFac::Instance().end())
      {
        std::cout << __func__ << ": WignerDMatrixEnv, missing id " << id; 
        throw std::string("WignerDMatrixEnv missed key"); 
      }

      return it->second.clone(); 
    }

    template<int X, int Y, int Z> 
      bool do_reg(const int J)
      {
        mom_t mom = gen_mom<X,Y,Z>(); 
        std::string id = gen_id(mom,J); 
        WignerMatrix_t W = gen_wigner_matrix(mom,J);  
        return reg_wigner_matrix(W,id); 
      }


    bool do_mom_reg(const int J)
    {
      bool success = true; 

      // D4 
      success &= do_reg<1,0,0>(J); 
      success &= do_reg<0,1,0>(J); 
      success &= do_reg<0,0,1>(J); 
      success &= do_reg<-1,0,0>(J); 
      success &= do_reg<0,-1,0>(J); 
      success &= do_reg<0,0,-1>(J); 

      success &= do_reg<2,0,0>(J); 
      success &= do_reg<0,2,0>(J); 
      success &= do_reg<0,0,2>(J); 
      success &= do_reg<-2,0,0>(J); 
      success &= do_reg<0,-2,0>(J); 
      success &= do_reg<0,0,-2>(J); 

      // D2      
      success &= do_reg<1,1,0>(J); 
      success &= do_reg<0,1,1>(J); 
      success &= do_reg<1,0,1>(J); 
      success &= do_reg<1,-1,0>(J); 
      success &= do_reg<0,1,-1>(J); 
      success &= do_reg<-1,0,1>(J); 
      success &= do_reg<-1,1,0>(J); 
      success &= do_reg<0,-1,1>(J); 
      success &= do_reg<1,0,-1>(J); 
      success &= do_reg<-1,-1,0>(J); 
      success &= do_reg<0,-1,-1>(J); 
      success &= do_reg<-1,0,-1>(J); 


      // D3
      success &= do_reg<1,1,1>(J); 
      success &= do_reg<-1,1,1>(J); 
      success &= do_reg<1,-1,1>(J); 
      success &= do_reg<1,1,-1>(J); 
      success &= do_reg<-1,-1,1>(J); 
      success &= do_reg<1,-1,-1>(J); 
      success &= do_reg<-1,1,-1>(J); 
      success &= do_reg<-1,-1,-1>(J); 

      return success; 
    }

  }  // anonomyous 



  namespace WignerDMatrixEnv
  {

    namespace
    {
      bool local_registration = false; 
    }

    bool registerAll(const int Jmax)
    {
      bool success = true; 

      if(!!! local_registration)
      {
        for(int J =0 ; J <= Jmax; ++J)
         success &= do_mom_reg(J);  

        local_registration = true; 
      }
      return success; 
    }

  WignerMatrix_t* 
    call_factory(const mom_t &p, const int J)
    {
      return query_factory(p,J); 
    }

  } // WignerDMatrixEnv


} // radmat

