/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_canonical_rotations_checker.cc

 * Purpose :

 * Creation Date : 24-12-2013

 * Last Modified : Mon 24 Mar 2014 03:29:24 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "lorentzff_canonical_rotations_checker.h"
#include "lorentzff_canonical_rotations_utils.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_JJlist_ff.h"
#include "semble/semble_semble.h"
#include <sstream>
#include <string> 
#include <complex>
#include <map>
#include <vector>
#include <exception>
#include "ensem/ensem.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/plot.h"
#include <math.h>
#include <iostream>
#include <complex>
#include "adat/handle.h"


namespace radmat
{


  namespace 
  {
    ////////////////////////////////////////////////
    std::string 
      string_mom(const mom_t &p) 
      {
        std::stringstream ss; 
        ss << p[0] << p[1] << p[2] ;
        return ss.str(); 
      }

    ////////////////////////////////////////////////
    mom_t 
      mul_tran(const RotationMatrix_t* R, const mom_t &x)
      {
        mom_t chk = gen_mom<0,0,0>();

        for(int i = 0; i < 3; ++i)
        {
          double res(0.); 
          for(int j = 0; j < 3; ++j)
            if ( fabs( (*R)[j+1][i+1] ) > 1e-6 )
              res += (*R)[j+1][i+1] * double(x[j]); 

          chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
        }
        return chk; 
      }


    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    struct KEY
    {
      KEY(const mom_t &l, const mom_t &r, const int hl, const int hr, const int mu)
        : left(l) , right(r) , h_left(hl) , h_right(hr), jmu(mu)
      {}

      std::string key(void) const
      {
        std::stringstream ss; 
        ss << "lefty_p" << string_mom(left) << ",H" <<h_left;
        ss << ".V" << jmu; 
        ss << ".righty_p" << string_mom(right) << ",H" <<h_right;
        return ss.str();
      }

      mom_t left;
      mom_t right; 
      int h_left; 
      int h_right; 
      int jmu; 
    };


    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    struct map_holder
    {
      typedef ENSEM::EnsemVectorComplex DATA; 

      void insert(const KEY &k , const DATA &d)
      {
        std::string key = k.key(); 
        if(keys_map.find(key) != keys_map.end())
        {
          std::cout << __func__ << ": k " << key << std::endl; 
          throw std::string("dub key err"); 
        }

        if(data_map.find(key) != data_map.end())
        {
          std::cout << __func__ << ": k " << key << std::endl; 
          throw std::string("dub dat err"); 
        }
        keys_map.insert(std::pair<std::string,KEY>(key,k)); 
        data_map.insert(std::pair<std::string,DATA>(key,d)); 
      }

      void init_map( const rHandle<LLSQLatticeMultiData> &d)
      {
        std::vector<ThreePointDataTag> tags = d->tags(); 

        for(int i =0; i < tags.size(); ++i)
        {
          KEY foo(tags[i].left_mom, tags[i].right_mom, tags[i].left_row, tags[i].right_row,tags[i].gamma_row); 
          DATA e = d->get_row_ensem(i); 
          insert(foo,e); 
        }
      } 

      DATA query(const KEY &k) const
      {
        std::string key = k.key(); 
        std::map<std::string,DATA>::const_iterator it; 
        it = data_map.find(key); 

        if( it == data_map.end())
          throw std::string("missing key"); 

        return it->second; 
      }

      std::map<std::string,KEY> keys_map; 
      std::map<std::string,DATA> data_map; 
    }; 

    ////////////////////////////////////////////////
    KEY can_key(const KEY &k, const int hl, const int hr, const int jmu)
    {
      DMatrixManager D;
      RotationMatrix_t *Rtriad; 
      Rtriad = D.triad_rotation_matrix(k.left,k.right);
      mom_t ll = mul_tran(Rtriad,k.left); 
      mom_t rr = mul_tran(Rtriad,k.right); 
      delete  Rtriad; 

      return KEY(ll,rr,hl,hr,jmu); 
    } 

    ////////////////////////////////////////////////
    Tensor<double,1>
      cook_p_tens(const double mass,
          const mom_t &p,
          const double kick)
      {
        Tensor<double,1> ret( (TensorShape<1>())[4], 0.); 
        ret[0] = sqrt( mass*mass + kick*kick*( p[0]*p[0] + p[1]*p[1] + p[2]*p[2])); 
        ret[1] = kick * p[0];
        ret[2] = kick * p[1];
        ret[3] = kick * p[2];
        return ret; 
      }

    ////////////////////////////////////////////////
    bool 
      test_equivalence( const ENSEM::EnsemVectorComplex &A, 
          const ENSEM::EnsemVectorComplex &B)
      {
        ENSEM::EnsemVectorComplex C = A - B; 

        ENSEM::EnsemVectorReal real, imag;

        double const_real, const_imag; 

        real = ENSEM::real(C);
        imag = ENSEM::imag(C);

        std::vector<double> t;
        for(int i = 0; i < C.numElem(); ++i)
          t.push_back(double(i)); 

        int thigh = t.size() * 0.85; 
        int tlow = t.size() * 0.15; 

        EnsemData ereal(t,real),eimag(t,imag); 
        ereal.hideDataAboveX(thigh - 0.1); 
        eimag.hideDataBelowX(tlow -0.1); 

        ADAT::Handle<FitFunction> freal(new ThreePointConstant), fimag(new ThreePointConstant);  
        JackFit fit_real(ereal,freal), fit_imag(eimag,fimag); 

        fit_real.runAvgFit(); 
        fit_imag.runAvgFit(); 
        const_real = fit_real.getAvgFitParValue(0);
        const_imag = fit_imag.getAvgFitParValue(0);

        double const_real_var = fit_real.getAvgFitParError(0);
        double const_imag_var = fit_imag.getAvgFitParError(0);

        // zero vy value
        if( fabs(const_real) < 1e-3 ) 
          if( fabs(const_imag) < 1e-3)
            return true; 

        // inverse covariance was singular
        if( isnan(const_real) )
          if( isnan(const_imag) )
            return true; 

        // zero w/ in err
        if( fabs(const_real) - 3.*sqrt(fabs(const_real_var)) < 0.)
          if( fabs(const_imag) - 3.*sqrt(fabs(const_imag_var)) < 0.)
            return true; 

        std::cout << __func__ << ": reporting failure, C = A-B (should be zero) " << std::endl;
        std::cout << __func__ << ": real = " << const_real << " +/- " << const_real_var << std::endl;
        std::cout << __func__ << ": imag = " << const_imag << " +/- " << const_imag_var << std::endl;
        ENSEM::write( std::string("C.jack") , C ); 


        return false; 
      }

  

    ////////////////////////////////////////////////
    
    template<int Jl, int Jr>
    bool 
      do_work(const map_holder *mappy, 
          const KEY &k)
      {
        JJFFimpl<Jl,Jr> foo; 
        Tensor<canIdx_t,1> tens = foo.canonicalize(
            std::make_pair( cook_p_tens(0.216,k.left,0.11), k.h_left), 
            std::make_pair( cook_p_tens(0.216,k.right,0.11), k.h_right),
            0.11);

        canIdx_t chk = tens[ k.jmu ]; 

        map_holder::DATA k_frame =  mappy->query(k); 
        map_holder::DATA can_sum;
        can_sum = SEMBLE::toScalar( std::complex<double>(0.,0.) )* k_frame; 
        canIdx_t::const_iterator it; 
        std::stringstream ss; 

        for( it = chk.begin(); it != chk.end(); ++it)
        {
          int hl = it->m_obj.h_left;
          int hr = it->m_obj.h_right;
          int jmu = it->m_obj.idx;
          KEY q_key = can_key(k,hl,hr,jmu);
          ss << q_key.key() << std::endl;
          map_holder::DATA tmp = mappy->query(q_key);
          can_sum = can_sum + SEMBLE::toScalar( it->m_coeff ) * tmp;  
        }

        bool success =  test_equivalence( k_frame, can_sum ); 

        if( !!! success )
        {
          ENSEM::write( k.key() , k_frame ) ;
          ENSEM::write( k.key() + std::string(".can_sum") , can_sum ); 

          std::cout << k.key() << " - ingredients:\n" <<  ss.str() << std::endl; 
        }

        return success; 
      }

    bool  do_work(const map_holder *mappy, const KEY &k, const int Jl, const int Jr)
    {
      if( (Jl == 0) && (Jr == 0) )
        return do_work<0,0>(mappy,k); 
      else if( (Jl == 1) && (Jr == 0) )
        return do_work<1,0>(mappy,k); 
      else if( (Jl == 2) && (Jr == 0) )
        return do_work<2,0>(mappy,k); 
      else if( (Jl == 3) && (Jr == 0) )
        return do_work<3,0>(mappy,k); 
      else if( (Jl == 0) && (Jr == 1) )
        return do_work<0,1>(mappy,k); 
      else if( (Jl == 1) && (Jr == 1) )
        return do_work<1,1>(mappy,k); 
      else if( (Jl == 2) && (Jr == 1) )
        return do_work<2,1>(mappy,k); 
      else if( (Jl == 3) && (Jr == 1) )
        return do_work<3,1>(mappy,k); 
      else if( (Jl == 0) && (Jr == 2) )
        return do_work<0,2>(mappy,k); 
      else if( (Jl == 1) && (Jr == 2) )
        return do_work<1,2>(mappy,k); 
      else if( (Jl == 2) && (Jr == 2) )
        return do_work<2,2>(mappy,k); 
      else if( (Jl == 3) && (Jr == 2) )
        return do_work<3,2>(mappy,k); 
      else if( (Jl == 0) && (Jr == 3) )
        return do_work<0,3>(mappy,k); 
      else if( (Jl == 1) && (Jr == 3) )
        return do_work<1,3>(mappy,k); 
      else if( (Jl == 2) && (Jr == 3) )
        return do_work<2,3>(mappy,k); 
      else if( (Jl == 3) && (Jr == 3) )
        return do_work<3,3>(mappy,k); 
      else
      {
        std::cerr << "unsuported spin" << std::endl;
        exit(1); 
      }
    }



    ////////////////////////////////////////////////
    void
      do_work(const rHandle<LLSQLatticeMultiData> &d, 
          const int Jl,
          const int Jr)
      {
        map_holder mappy;
        try
        {
          mappy.init_map(d); 
        }
        catch(std::string &s)
        {
          std::cout << __func__ << ": caught " << s << std::endl;
          exit(1); 
        }

        std::map<std::string,KEY>::const_iterator it; 
        for (it = mappy.keys_map.begin(); it != mappy.keys_map.end(); ++it)
        {
          if( !!! do_work(&mappy,it->second,Jl,Jr) )
          {
            std::cout << it->first << ": failed!!!" << std::endl;
          }
        }

      }

  } // anonomyous 


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  void
    LatticeRotationRelationChecker::check(
        const rHandle<LLSQLatticeMultiData> &d, 
        const int Jl, 
        const int Jr) const
    {
      do_work(d,Jl,Jr); 
    }



}

