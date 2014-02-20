#ifndef FORMFACTOR_KINEMATIC_FACTORS_H
#define FORMFACTOR_KINEMATIC_FACTORS_H 

#include "formfactor_abs_base_cfg.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include "semble/semble_algebra.h"
#include "semble/semble_meta.h"
#include <utility>
#include <sstream>
#include <string>
#include <list>

namespace radmat
{

  // the kinematic invariants that go into the calculation
  ////////////////////////////////////////////////////////
  struct FFSingleKinematicInvariants
  {
    // optional extra string field -- not currently in use
    FFSingleKinematicInvariants( const ENSEM::EnsemReal &energy, 
        const Array<int> &mom,
        const int rrow,
        const double mom_factor,
        const std::string &id)
    {
      init(energy,mom,rrow,mom_factor,id);
    }

    // no extra field
    FFSingleKinematicInvariants( const ENSEM::EnsemReal &energy, 
        const Array<int> &mom,
        const int rrow,
        const double mom_factor)
    {
      init(energy,mom,rrow,mom_factor);
    }

    // default id field
    void init( const EnsemReal &energy, 
        const Array<int> &mom,
        const int rrow,
        const double mom_factor,
        const std::string id="default")
    {
      E = energy;  
      p.resize(3);
      p[0] = mom_factor*mom[0]; 
      p[1] = mom_factor*mom[1];
      p[2] = mom_factor*mom[2]; 
      row = rrow; 
      LG_id = id;
    }

    Array<double> p; 
    ENSEM::EnsemReal E; 
    std::string LG_id;  
    int row;          
  }; 


  // the kinematic invariants that go into the calculation
  ////////////////////////////////////////////////////////
  struct FFKinematicInvariants
  {
    // squish two singles together 
    FFKinematicInvariants( const FFSingleKinematicInvariants &llefty,
        const FFSingleKinematicInvariants &rrighty,
        const double mom_factorr)
      : lefty(llefty) , righty(rrighty) , mom_factor(mom_factorr)
    { }

    // what is the true Q2 here?
    ENSEM::EnsemReal Q2(void)
    {
      ENSEM::EnsemReal q0,a,b,c; 
      q0 = lefty.E - righty.E; 
      a = q0;
      b = q0; 
      c = q0; 
      a = SEMBLE::toScalar( lefty.p[0] - righty.p[0] );
      b = SEMBLE::toScalar( lefty.p[1] - righty.p[1] );
      c = SEMBLE::toScalar( lefty.p[2] - righty.p[2] );
      
      return (-q0*q0 + a*a + b*b + c*c);
    }

    // for strings
    double Q2_doub(void)
    {
      return SEMBLE::toScalar(ENSEM::mean(Q2())); 
    }

    void set_mom_factor(const double p_fac) {mom_factor = p_fac;}

    FFSingleKinematicInvariants lefty,righty;
    double mom_factor; 
  };



  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  struct FFKinematicFactors_t
  {
    // save some typing
    typedef rHandle<FFAbsBase_t > FFBase_h;
    typedef SEMBLE::SembleMatrix<std::complex<double> > KinematicFactorMatrix;

    // the only available constructor
    FFKinematicFactors_t(const FFBase_h &KFacGen)
      : m_KFacGen(KFacGen) 
    {  }

    ~FFKinematicFactors_t(void) {} // handle cleans itself up

    // basically generate the 4 X (n multipole) matrix of kinematic factors, 
    // the row index is the lorentz index of the lattice matrix element
    KinematicFactorMatrix genFactors(const FFKinematicInvariants &inv)
    {
      FFSingleKinematicInvariants lefty(inv.lefty),righty(inv.righty); 

      int nfacs = m_KFacGen->nFacs();
      int nbins = righty.E.size();

      // this also checks size of righty vs size of lefty
      POW2_ASSERT_DEBUG( (nbins == lefty.E.size()) && (nfacs > 0) );

      // the matrix to be returned
      KinematicFactorMatrix KF(nbins,4,nfacs);
      KF.zeros();

      // scale down the energies, momentum have zero variance 
      ENSEM::rescaleEnsemDown( lefty.E );
      ENSEM::rescaleEnsemDown( righty.E ); 

      // loop over cfgs and use the wrapper to operator()
      for(int bin = 0; bin < nbins; bin++)
        KF[bin] = m_KFacGen->wrapper( lefty.E.elem(bin) , lefty.p , lefty.row,
            righty.E.elem(bin) , righty.p , righty.row , inv.mom_factor);

      // scale up return matrix
      KF.rescaleSembleUp();

      return KF;
    }

    int nFacs(void) {return m_KFacGen->nFacs();}

    // hide ctor    
    private:
    FFKinematicFactors_t(void);
    FFKinematicFactors_t(const FFKinematicFactors_t &o);
    FFKinematicFactors_t& operator=(const FFKinematicFactors_t &o);

    private:
    FFBase_h m_KFacGen;
  };

} // radmat

#endif /* FORMFACTOR_KINEMATIC_FACTORS_H */
