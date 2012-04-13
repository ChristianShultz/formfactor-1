#ifndef POLARISATION_FACTORY_H_H_GUARD
#define POLARISATION_FACTORY_H_H_GUARD

#include "hadron/clebsch.h"
#include "tensor/tensorbase.h"
#include "polarisation_factory_keys.h"
#include "polarisation_factory_inventory.h"
#include "utils/rodrigues_rotation_matrix.h"
#include "itpp/itbase.h"
#include <complex>

/**
   @file polarisation_factory.h
   @brief contains the factory and helper classes
 */


namespace polarisation
{
  /**
     @brief the factory, this class boils down to a container for member functions

     @details the factory, this class boils down to a container for member functions to 
     interface easily with the factory inventory and make syntax more readable
    
     The factory is only responsible for providing an interface to get the data. The actual creation 
     of polarisation tensors is managed by the Coupler<J> class.  

     An example to illustrate (helicity sums are implicit):

     1) we call factory ( key<J=3> ) on an empty factory.
     factory( key<J=3> ) calls coupler<J=3> which in turn calls
     factory(key<J=1>) and factory(key<J=2>) 

     2.a) factory(key<j=1>) calls the base template coupler<J=1> which makes the 
     polarisation vector for whatever helicity was in the key

     2.b) the call to factory(key<J=2>) calls two coupler<J=1> and combines them in the appropriate manner (summing
     on the helicity) to  make a J=2 tensor (of a given helicity).  This is returned to the original call in 1) and combined
     with the result from 2.a) to make the J=3 tensor

     as all of the recursive calls were made each intermediate tensor was registered to a global static repository so that 
     if we later call factory( key<J=4> ) we already have saved versions of J=3,2,1  and don't have to 
     repeate the recursive steps through to make all of the sub objects.
   */
  struct pFac
  {  
    pFac(void);                                                              //! do nothing constructor
    tensor::TensorImplBase* get(const pFacKey &k);                           //! get the data for a key
    tensor::TensorImplBase* operator()(const pFacKey &k) {return get(k);}    //! get the data for a key
    void dumpFactory(void) const;                                            //! print contents of factory (testing)
  };




  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /**
     @brief template recursive coupler to make spin J = N by coupling J = 1 to J = N-1
     @details this is the general case, there is a specialization below for J = 1
  */
  template<idx_t J>
  struct Coupler
  {
    //! allow only pFac to use this class since the constructor is private
    friend struct pFac;
  private:
    pFacKey k;      //! the key we need to make
    pFac * factory; //! ptr to a factory for intermediary registration

    //! private constructor
    Coupler(const pFacKey &_k, pFac *_factory)
      : k(_k) , factory(_factory)
    {   }

    //! make spin J = N by coupling J = 1 to J = N-1
    tensor::Tensor<std::complex<double> , J >* operator()(void)
    {
      // set up some keys for recursion 
      pFacKey k_Jm1_h(k), k_1_h1(k);
      k_Jm1_h.J = J - 1;
      k_1_h1.J = 1;
      
      // set up the tensor that we need to make
      tensor::Tensor<std::complex<double> , J> ret;
      std::vector<idx_t> dim(J,4);
      ret.create(&dim[0]);

      // sum on the helicity components of J-1 and 1 to make total spin J, helicity k.helicity
      for(short h = -(J-1); h < J; ++h)
	{
	  k_Jm1_h.helicity = h;
	  for(short h1 = -1; h1 < 2; ++h1)
	    {
	      k_1_h1.helicity = h1;
	      double cg = Hadron::clebsch(2,2*h1,2*(J-1),2*h,2*J,2*k.helicity);
	      
	      // obviously don't bother if the coefficient is zero 
	      if(cg != 0.)
		{
		  // call the factory to get the tensors we need, if it doesnt exist then the 
		  // factory will recall this function and will recurse to make the bits we want
		  // until it gets to J = 1 which is specialized. operator^(Tensor,Tensor) is 
		  // overloaded to be a tensor product operator

		  ret+= tensor::castAndDeref<std::complex<double> , J - 1>(factory->get(k_Jm1_h))
		    ^ tensor::castAndDeref<std::complex<double> , 1>(factory->get(k_1_h1))
		    * cg;
		}
	    }
	}
    
      // return a pointer to a new tensor for registration by the factory   
      return ret.clone();
    }
  };  // close Coupler<J>



  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /**
     @brief J = 1 base case, make a J = 1 polarisation 4-vector along a given direction
     @details this guy is public for testing purposes so we can explicitly check that this 
     stupidly constructed design actually gets the correct answer
  */
  template<>
  struct Coupler<1>
  {
    //! doesn't actually need to be here
    friend struct pFac;

    pFacKey k;        //! the key for J = 1 that we want to make
    pFac * factory;   //! pointer to the factory for intermediary registration

    //! constructor
    Coupler(const pFacKey &_k, pFac *_factory)
      : k(_k) , factory(_factory)
    {   }
    
    /**
       @brief make the polarisation 4-vector along a given momentum direction
       @details will break for massless particles and helicity = 0, it starts by 
       making a polarisation vector along the positive z-axis then rotating the spatial compontents 
       to 'point' along the direction of the specified momentum. Its possible that we may wish to 
       make the vector along some direction other than the z-axis first, in this case we will register
       the intermediary vector along the z-axis for future factory use.       
     */
    tensor::Tensor<std::complex<double> , 1> * operator()(void)
    {
      pFacKey k_z(k);
      
      double modp = sqrt(k[1]*k[1] +
			 k[2]*k[2] +
			 k[3]*k[3]
			 );
      k_z[1] = 0.;
      k_z[2] = 0.;
      k_z[3] = modp;
      
      tensor::Tensor<std::complex<double> , 1> eps_z;

      if(polarisation::pFacInv::pFacInvIsRegistered(k_z))  // check for eps_z -- do we have the k_z data, skip some work if we do 
	{
	  eps_z = tensor::castAndDeref<std::complex<double> , 1>(polarisation::pFacInv::pFacInvGetData(k_z));
	}
      else   	  // if not we will create and possibly register it
	{
	  POW2_ASSERT(k_z.J == 1);                         // sanity check that we didn't get here through stupidity
	  short helicity = k_z.helicity;
	  eps_z.create(std::vector<idx_t>(1,4));

	  // the definition of the polarisation vector for each of the three possible helicities
	  if(helicity == -1)
	    {
	      eps_z[0] = 0.;
	      eps_z[1] = 1./sqrt(2.);
	      eps_z[2] = std::complex<double>(0.,-1.);
	      eps_z[3] = 0.;
	    }
	  else if(helicity == 0)
	    {
	      double p = k_z[3];
	      double E = k_z[0];
	      double m = sqrt(E*E - p*p);

	      POW2_ASSERT(m > 1e-15); // massless particles only have two polarization states h = +/-
	                              // the masses in the calculationsare expected to be order 1 
	      eps_z[0] = p/m;
	      eps_z[1] = 0.;
	      eps_z[2] = 0.;
	      eps_z[3] = E/m;
	    }
	  else if(helicity == 1)
	    {
	      eps_z[0] = 0.;
	      eps_z[1] = - 1./sqrt(2.);
	      eps_z[2] = std::complex<double>(0.,-1.);
	      eps_z[3] = 0.;
	    }
	  else
	    POW2_ASSERT(false); // the helicity was not possible for J = 1 -- something went horribly wrong
	}
      
      // return early if we wanted k_z 
      if(k == k_z)
	return eps_z.clone();
      
      // register the k_z guy to save work later
      polarisation::pFacInv::pFacInvRegister(k_z,eps_z.clone());

      // determine the rotation matrix between the two momentum vectors
      itpp::Vec<double> pz(3),p(3);
      itpp::Mat<double> R;
      pz.zeros();
      p.zeros();
      pz[2] = modp;
      p[0] = k[1];
      p[1] = k[2];
      p[2] = k[3];

      // solve for the rotation matrix between p and pz      
      R = rodRotMat(pz,p);

      // set up the polarization vector along p 
      tensor::Tensor<std::complex<double>,1> eps;
      eps = eps_z;

      // rotate the spatial components of eps_z to get eps
      for(idx_t i = 0; i < 3; ++i)
	{
	  eps[i+1] = 0.;
	  eps[i+1] += R(i,0) * eps_z[1];
	  eps[i+1] += R(i,1) * eps_z[2];
	  eps[i+1] += R(i,2) * eps_z[3];
	}
      
      // return a pointer to new obj for registration by the factory
      return eps.clone();
    } // close operator()

  }; // close Coupler<J=1>

}


#endif
