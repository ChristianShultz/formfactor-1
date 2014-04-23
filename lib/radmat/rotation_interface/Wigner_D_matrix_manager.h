#ifndef WIGNER_D_MATRIX_MANAGER_H
#define WIGNER_D_MATRIX_MANAGER_H 



#include "Wigner_D_matrix_factory.h"
#include "rotation_utils.h"
#include "hadron/irrep_util.h"
#include "radmat/utils/tensor.h"
#include <complex>
#include <sstream>

namespace radmat
{

  // most positive helicity is mapped to zero
  struct DMatrixManager
  {
    virtual ~DMatrixManager(void) {}

    virtual std::string 
      string_mom(const mom_t &p) const
    {
      std::stringstream ss; 
      ss << p[0] << p[1] << p[2] ;
      return ss.str(); 
    }

    virtual std::pair<mom_t,mom_t> 
      get_frame(const mom_t &l, 
          const mom_t &r) const ;

    virtual void 
      check_throw_frame_err(const RotationMatrix_t* Rtriad, 
          const std::pair<mom_t,mom_t> &f, 
          const std::pair<mom_t,mom_t> &can) const; 

    virtual WignerMatrix_t* 
      get_can_mat(const mom_t &p, 
          const int J) const;

    virtual void 
      conjugate(WignerMatrix_t * D) const;

    virtual void 
      transpose(WignerMatrix_t *D) const; 

    virtual void
      dagger(WignerMatrix_t *D) const; 

    virtual void 
      clean(WignerMatrix_t *D, 
          const double thresh=1e-6) const;  

    virtual RotationMatrix_t*
      rotation_matrix(const mom_t &l, 
          const mom_t &r) const; 

    virtual WignerMatrix_t*
      wigner_matrix(const RotationMatrix_t *R, 
          const int J) const;

    virtual WignerMatrix_t*
      left_wigner_matrix(const RotationMatrix_t *R, 
          const mom_t &l, 
          const mom_t &r, 
          const int J,
          bool use_map=true) const; 

    virtual WignerMatrix_t*
      right_wigner_matrix(const RotationMatrix_t *R,
          const mom_t &l, 
          const mom_t &r, 
          const int J,
          bool use_map=true) const; 
  }; 

  namespace WignerThreadMapEnv
  {
    bool registerAll(); 
  }

} // radmat



#endif /* WIGNER_D_MATRIX_MANAGER_H */