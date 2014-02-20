#ifndef FORMFACTOR_ABS_BASE_CFG_H
#define FORMFACTOR_ABS_BASE_CFG_H 

#include "formfactor_utils.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include <complex>
#include <utility>
#include <sstream>
#include <string>
#include <list>


namespace radmat
{

  enum FFMODE
  {
    HELICITY, 
    CUBIC
  };


  typedef std::pair< Tensor<double,1> , int > MomHelPair_t; 

  // this is a base class to set up list of handles 
  // to polymorphic functors to construct a generalization
  // of a matrix element
  //
  // actual classes will implement the operator() 
  //
  // a decomposition is a list of functors

  struct FFAbsBlockBase_t;
  REGISTER_STRINGIFY_TYPE(FFAbsBlockBase_t); 

  struct FFAbsBlockBase_t
  {
    virtual ~FFAbsBlockBase_t() {}
    virtual std::string ff() const {return Stringify<FFAbsBlockBase_t>();}


    virtual Tensor<std::complex<double> , 1> 
      operator()( const MomHelPair_t &lefty, 
          const MomHelPair_t &righty,
          const double mom_fac) const = 0; 
  };

  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////


  // generate some linear least squares system w/o having to know very much about the 
  // matrix element invariants. These are still at the configuration level. Just 
  // shove this into a SembleMatrix bin by bin to deal with the ensemble stats.

  // this will be polymorphic with different named classes corresponding to 
  // different quantum numbers, ie: PiPi  <-->   < 0+ | j_mu | 0+ >

  struct FFAbsBase_t; 
  REGISTER_STRINGIFY_TYPE( FFAbsBase_t );

  struct FFAbsBase_t
  {
    // save some typing
    typedef FFAbsBlockBase_t BBType;
    typedef rHandle< FFAbsBlockBase_t > BBHandle_t;
    typedef std::list< BBHandle_t > FFAbs_list;

    // this will be useful when we derive
    FFAbsBase_t(const FFAbs_list& list)
      : m_list(list) 
    {  }

    FFAbsBase_t& operator=(const FFAbsBase_t &o)
    {
      if(this != &o)
      {
        m_list = o.m_list;
      }
      return *this;
    }

    FFAbsBase_t(const FFAbsBase_t &o)
      : m_list(o.m_list)
    {  }

    // needs to be present and virtual b/c we are using pointers to derived
    virtual ~FFAbsBase_t(void) {}

    // useful higher up
    virtual int nFacs(void) {return m_list.size();}

    // generate some tex code corresponding to thestring of stuff we think this is making
    std::string ff(void) const 
    {
      std::stringstream ss;
      FFAbs_list::const_iterator it;
      for (it = m_list.begin(); it != m_list.end(); it++)
        ss << (*it)->ff() << "  ";
      return ss.str();
    }

    // matrix elem id
    virtual std::string id() const { return Stringify<FFAbsBase_t>(); }
    virtual int left_spin() const = 0; 
    virtual int right_spin() const = 0; 

    // generate the linear system based on the available set of kinematic factors
    virtual itpp::Mat<std::complex<double> > 
      operator()( const MomHelPair_t &lefty, 
          const MomHelPair_t &righty,
          const double mom_fac) const
      {
        itpp::Mat<std::complex<double> > ret;

        FFAbs_list::const_iterator it;
        for (it = m_list.begin(); it != m_list.end(); it++)
          ret.append_col(toItpp<std::complex<double> >((**it)(lefty,righty,mom_fac)));

        return ret;
      }

    MomHelPair_t to_mom_hel_pair( const double E, 
        const Array<double> &p, 
        const int h)
    {
      Tensor<double,1> mom( (TensorShape<1>())[4] , 0. );
      mom[0] = E;
      mom[1] = p[0];
      mom[2] = p[1];
      mom[3] = p[2];

      return MomHelPair_t(mom,h); 
    }

    // wrap the call to operator() so we have a condensed
    // version of the information later
    virtual itpp::Mat<std::complex<double> >
      wrapper(const ENSEM::Real El, 
        const Array<double> &pl,
        const int hl, 
        const ENSEM::Real Er, 
        const Array<double> &pr,
        const int hr,
        const double mom_fac)
    {
      return this->operator()(to_mom_hel_pair(ENSEM::toDouble(El),pl,hl),
          to_mom_hel_pair(ENSEM::toDouble(Er),pr,hr),
          mom_fac); 
    }

    protected:  // hide ctor
    FFAbsBase_t(void);

    // data store
    FFAbs_list m_list;
  };


} // radmat 

#endif /* FORMFACTOR_ABS_BASE_CFG_H */
