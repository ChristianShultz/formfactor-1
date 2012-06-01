#ifndef POLARISATION_TENSORS_H_H_GUARD
#define POLARISATION_TENSORS_H_H_GUARD

#include "tensor.h"
#include "aux.h"
#include "adat/handle.h"
#include "ensem/ensem.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "xml_array.h"
#include <utility>
#include <map>

namespace radmat
{

  typedef std::pair<idx_t, short> pKey_t; // J, helicity
  
  // fwd of a recursive template pattern for generating polarisation tensors
  template<idx_t J>
  struct genPolTens;

  /*
    A genPolTens is a recursive way of generating polarisation tensors in a
    mildly efficient maner.  Use of handles allows sharing of data between the 
    parent instance (rank N) and its children (rank N-1 .. 1) so that we only have
    to generate the intermediate tensors one time at most.
   */

  // J = 0 specialization
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<>
  struct genPolTens<0>
  {  }; // empty

  // J = 1 specialization
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  // gcc is frustrating .. it won't let a many to one friendship pattern without
  // totally rewriting the templates so I'm leaving genPolTens<1> completly public
  template<>
  struct genPolTens<1>
  {
    // save some typing
    typedef ADAT::Handle<Tensor<double,2> > R_ij;
    typedef XMLArray::Array<int> mom_t;
    typedef std::map<pKey_t,TensorBase*> map_t;
    typedef ADAT::Handle<map_t> map_handle;
    typedef ADAT::Handle<mom_t> mom_handle;

    genPolTens(void); // hidden

    //! constructor taking momentum -- tells us how to do the rotations via adat
    genPolTens(const mom_t &_p)
    : h_map(new map_t()) , R(genRotationMatrix(_p).clone()) , h_mom(new mom_t(_p))
    { }

    //! get a polarisation tensor for the inputs in the direction of p
    Tensor<std::complex<double>,1> operator()(const double E, const short hel)
    {
      return make(E,hel);
    }

    //!  do not use!
    genPolTens(map_handle &_h_map, R_ij _R, mom_handle & _h_mom )
      : h_map(_h_map) , R(_R) , h_mom(_h_mom)
    { }

    //! do not use!
    Tensor<std::complex<double> , 1> get(const double E, const short hel)
    {
      map_t::const_iterator it;
      it = h_map->find(pKey_t(1,hel));
      if(it != h_map->end())
	return m_downcast(it->second);
   
      h_map->insert(map_t::value_type(pKey_t(1,hel),make(E,hel).clone()));
      return get(E,hel);
    }

  protected:
    //! do not use!
    Tensor<std::complex<double> , 1> m_downcast(TensorBase* bar)
    {
      Tensor<std::complex<double> , 1> * foo = dynamic_cast<Tensor<std::complex<double>, 1 >* >(bar);
      POW2_ASSERT(foo);
      return *foo;
    }

    //! make a J = 1 polarisation tensor
    Tensor<std::complex<double>,1> make(const double E, const short hel)
    {
      Tensor<std::complex<double> ,1 > eps(TensorShape<1>()[4],0.);
      bool in_bounds(false);
      if(hel == 1)
	{
	  in_bounds = true;
	  eps[1] = -1./sqrt(2.);
	  eps[2] = std::complex<double>(0.,-1./sqrt(2.));
	}
      if(hel == -1)
	{
	  in_bounds = true;
	  eps[1] = 1./sqrt(2.);
	  eps[2] = std::complex<double>(0.,-1./sqrt(2.));
	}
      if(hel == 0)
	{
	  in_bounds = true;
	  double m_p = sqrt( (*h_mom)[0] * (*h_mom)[0] 
			     +  (*h_mom)[1] * (*h_mom)[1] 
			     +  (*h_mom)[2] * (*h_mom)[2] );
	  double m_m = sqrt(E * E - m_p * m_p);

	  // sanity and nan check , in lattice units m_pi should obviously be larger than this hardwire
	  POW2_ASSERT(m_m > 1e-14);
	  eps[0] = m_p / m_m;
	  eps[3] = E / m_m;
	}
      
      // sanity
      POW2_ASSERT(in_bounds);

      // rotate then return
      return *R*eps; 
    }

    // shared data store    
    map_handle h_map;
    R_ij R;
    mom_handle h_mom;
  };

  // recursive defs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<idx_t J>
  struct genPolTens
  {
    typedef ADAT::Handle<Tensor<double,2> > R_ij;
    typedef XMLArray::Array<int> mom_t;
    typedef std::map<pKey_t,TensorBase*> map_t;
    typedef ADAT::Handle<map_t> map_handle;
    typedef ADAT::Handle<mom_t> mom_handle;

    genPolTens(void); // hidden

    //! constructor taking momentum -- tells us how to do the rotations via adat
    genPolTens(const mom_t &_p)
    : h_map(new map_t()) , R(genRotationMatrix(_p).clone()) , h_mom(new mom_t(_p))
    { }

    //! return a polarisation tensor of rank J for the inputs
    Tensor<std::complex<double>, J> operator()(const double E, const short hel)
    {
      refresh_map();
      return get(E,hel);
    }

  private:
    // clean up the map each time we get a new input
    void refresh_map(void)
    {
      map_t::iterator it;
      for(it = h_map->begin(); it != h_map->end(); it++)
	{
	  delete it->second;
	  it->second = NULL;
	}
      h_map->clear();
    }

    // wrap make in a sensible form
    Tensor<std::complex<double> , J> get(const double E, const short hel)
    {
      map_t::const_iterator it;
      it = h_map->find(pKey_t(J,hel));
      if(it != h_map->end())
	return m_downcast(it->second);

      make(E,hel);
      return get(E,hel);
    }

    // downcast the base pointers to the derived type and return an actual obj
    Tensor<std::complex<double>, J> m_downcast(TensorBase * bar)
    {
      Tensor<std::complex<double> , J> * foo = dynamic_cast<Tensor<std::complex<double>, J >* >(bar);
      POW2_ASSERT(foo);
      return *foo;
    }

    // couple them together using the SU(2) CG found in adat
    void make(const double E, const short hel)
    {
      genPolTens < J - 1 > Jm(h_map, R, h_mom);
      genPolTens<1> J1(h_map, R, h_mom);
      double factor;

      Tensor<std::complex<double> , J> ptensor(std::vector<idx_t>(J, 4), 0.);

      for(short big_h = -(J-1); big_h < J ; ++big_h)
	for(short small_h = -1; small_h < 2; ++small_h)
	  {
	    factor = Hadron::clebsch(2 * (J - 1), 2 * big_h, 2, 2 * small_h, 2 * J, 2 * hel);

	    if(factor != 0.)
	      ptensor += (factor * (Jm.get(E,big_h) ^ J1.get(E,small_h)));
	  }

      h_map->insert(map_t::value_type(pKey_t(J, hel), ptensor.clone()));
    }

    //! private constructor to share parent's information with children 
    genPolTens(map_handle &_h_map, R_ij _R, mom_handle &_h_mom)
      : h_map(_h_map) , R(_R) , h_mom(_h_mom)
    { }

    friend struct genPolTens<J+1>;

    map_handle h_map;
    R_ij R;
    mom_handle h_mom;
  };

  
} // radmat

#endif
