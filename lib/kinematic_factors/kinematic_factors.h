#ifndef KINFAC_ABC_H_H_GUARD
#define KINFAC_ABC_H_H_GUARD

#include "kinfac_fwd.h"
#include "kinfac_ini.h"
#include "kinfac_key.h"
#include "ensem/ensem.h"
#include "utils/pow2assert.h"
#include "semble/semble_vector.h"

/**
   @file kinematic_factors.h
   @brief contains the classes for creating the linear systems
 */

namespace kinfac
{
    /**
       @brief holds what will be the 'rows' of the linear system
       @details actually holds all of the lorentz indicies for a
       given row of the linear system
     */

    template<typename T>
    struct EnsembleRowHolder
    {
        //! construct and allocate storage for a 4 elem rank 1 tensor
        EnsembleRowHolder(void)
        {
            data.create(std::vector<tensor::idx_t>(1, 4));
        }

      //! do nothing destructor
      ~EnsembleRowHolder(void) { }

        //! put data into the holder
        void setLorentzIndex(const tensor::idx_t idx, 
			     SEMBLE::SembleVector<T> &ensembleOfKinematicFactors)
        {
            data[idx] = ensembleOfKinematicFactors;
        }

        //! get a reference to a lorentz index
        typename SEMBLE::SembleVector<T> &operator[](const tensor::idx_t idx)
        {
            return getLorentzIndex(idx);
        }

        //! get a const reference to a lorentz index
        const SEMBLE::SembleVector<T> &operator[](const tensor::idx_t idx) const
        {
            return getLorentzIndex(idx);
        }

        //! get a reference to a lorentz index
        typename SEMBLE::SembleVector<T> &getLorentzIndex(const tensor::idx_t idx)
        {
            return data[idx];
        }

        //! get a const reference to a lorentz index
        const typename SEMBLE::SembleVector<T> &getLorentzIndex(const tensor::idx_t idx) const
        {
            return data[idx];
        }

        //! resize the std::vectors to the length of the ensemble
      void reDim(const tensor::idx_t ncfgs, const tensor::idx_t nelem)
        {
	  data[0].reDim(ncfgs,nelem);
	  data[1].reDim(ncfgs,nelem);
	  data[2].reDim(ncfgs,nelem);
	  data[3].reDim(ncfgs,nelem);
        }

      void jackRescaleUp(void)
      {
	data[0].rescaleSembleUp();
	data[1].rescaleSembleUp();
	data[2].rescaleSembleUp();
	data[3].rescaleSembleUp();
      }

      void jackRescaleDown(void)
      {
	data[0].rescaleSembleDown();
	data[1].rescaleSembleDown();
	data[2].rescaleSembleDown();
	data[3].rescaleSembleDown();
      }

    private:
        /** @brief tensor index is lorentz index,
        @detials vector index is cfg (ensemble) index,
        itpp::Vec<T> is the 'row' of kinematic factors
        */
      tensor::Tensor<SEMBLE::SembleVector<T>, 1 > data;
    };


    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////



    //! we might want to change this later for some reason so stick in a typedef
    template<typename T>
    struct LinSysRetType
    {
        typedef EnsembleRowHolder<T> type;
    };


    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////


    /**
       @brief the abstract base class for creating linear systems
     */

    template <typename T>
    struct KinFacABC
    {
        virtual typename LinSysRetType<T>::type operator()(const KinKey &key) = 0;
    };


    // 0 -> 0
    ////////////////////////////////////////////////////////////////////////////

    /**
       @brief \f$ \langle 0^+ | j^{\mu} | 0^+ \rangle \f$
    */

    template<typename T>
    struct kf0p0p : public KinFacABC<T>
    {
        //!  \f$ \langle 0^+ | j^{\mu} | 0^+ \rangle \f$ kinematic factors
        /*!  \f$ M^{\mu} \equiv  \langle 0^+(p_f) | j^{\mu} | 0^+(p_i) \rangle \f$,
          the general decomposition is

          \f{equation*} M^{\mu} = F_0 p_+^{\mu} + F_1 p_-^{\mu} \f}
          case A: \f$ m_i = m_f \f$
          \f{eqnarray*}{
          q_{\mu} M^{\mu} &=& 0 \\
          &\rightarrow & F_0(m_i^2 - m_f^2) + F_1q^2 = 0 \\
          &\therefore & F_1 = 0 \\
          M^{\mu} &=& F_0 p_+^{\mu}
          \f}
          case B:  \f$ m_i \neq m_f \f$
          \f{align*}{
          q_{\mu}M^{\mu} &= 0 \\
          &\rightarrow F_0(m_i^2 - m_f^2) + F_1q^2 = 0 \\
          &\therefore F_0 = F_1 \frac{Q^2}{m_i^2 -m_f^2} \\
          M^{\mu} &= F_1 \left(p_+^{\mu} \frac{Q^2}{m_i^2 -m_f^2} + p_-^{\mu} \right)
          \f}
        */
      
      kf0p0p(void); //! not implemented
      kf0p0p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 0^+ | j^{\mu} | 0^- \rangle \f$
    */

    template<typename T>
    struct kf0p0m : public KinFacABC<T>
    {
      kf0p0m(void); //! not implemented
      kf0p0m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    // 1 -> 1
    ////////////////////////////////////////////////////////////////////////////

    /**
       @brief \f$ \langle 1^+ | j^{\mu} | 1^+ \rangle \f$
    */

    template<typename T>
    struct kf1p1p : public KinFacABC<T>
    {
      kf1p1p(void); //! not implemented
      kf1p1p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 1^+ | j^{\mu} | 1^- \rangle \f$
    */

    template<typename T>
    struct kf1p1m : public KinFacABC<T>
    {
      kf1p1m(void); //! not implemented
      kf1p1m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    // 2 -> 2
    ////////////////////////////////////////////////////////////////////////////

    /**
       @brief \f$ \langle 2^+ | j^{\mu} | 2^+ \rangle \f$
    */

    template<typename T>
    struct kf2p2p : public KinFacABC<T>
    {     
      kf2p2p(void); //! not implemented
      kf2p2p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 2^+ | j^{\mu} | 2^- \rangle \f$
    */

    template<typename T>
    struct kf2p2m : public KinFacABC<T>
    {
      kf2p2m(void); //! not implemented
      kf2p2m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };





    // 0 -> 1
    ////////////////////////////////////////////////////////////////////////////

    /**
       @brief \f$ \langle 0^+ | j^{\mu} | 1^+ \rangle \f$
    */

    template<typename T>
    struct kf0p1p : public KinFacABC<T>
    {
      kf0p1p(void); //! not implemented
      kf0p1p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 0^+ | j^{\mu} | 1^- \rangle \f$
    */

    template<typename T>
    struct kf0p1m : public KinFacABC<T>
    {
      kf0p1m(void); //! not implemented
      kf0p1m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;

        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 0^- | j^{\mu} | 1^- \rangle \f$
    */

    template<typename T>
    struct kf0m1m : public KinFacABC<T>
    {
      kf0m1m(void); //! not implemented
      kf0m1m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };

    /**
       @brief \f$ \langle 0^- | j^{\mu} | 1^+ \rangle \f$
    */

    template<typename T>
    struct kf0m1p : public KinFacABC<T>
    {
      kf0m1p(void); //! not implemented
      kf0m1p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    // 0 -> 2
    ////////////////////////////////////////////////////////////////////////////

    /**
       @brief \f$ \langle 0^+ | j^{\mu} | 2^+ \rangle \f$
    */

    template<typename T>
    struct kf0p2p : public KinFacABC<T>
    {
      kf0p2p(void); //! not implemented
      kf0p2p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 0^+ | j^{\mu} | 2^- \rangle \f$
    */

    template<typename T>
    struct kf0p2m : public KinFacABC<T>
    {
      kf0p2m(void); //! not implemented
      kf0p2m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 0^- | j^{\mu} | 2^- \rangle \f$
    */

    template<typename T>
    struct kf0m2m : public KinFacABC<T>
    {
      kf0m2m(void); //! not implemented
      kf0m2m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };

    /**
       @brief \f$ \langle 0^- | j^{\mu} | 2^+ \rangle \f$
    */

    template<typename T>
    struct kf0m2p : public KinFacABC<T>
    {
      kf0m2p(void); //! not implemented
      kf0m2p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    // 1 -> 2
    ////////////////////////////////////////////////////////////////////////////

    /**
       @brief \f$ \langle 1^+ | j^{\mu} | 2^+ \rangle \f$
    */

    template<typename T>
    struct kf1p2p : public KinFacABC<T>
    {
      kf1p2p(void); //! not implemented
      kf1p2p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 1^+ | j^{\mu} | 2^- \rangle \f$
    */

    template<typename T>
    struct kf1p2m : public KinFacABC<T>
    {
      kf1p2m(void); //! not implemented
      kf1p2m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    /**
       @brief \f$ \langle 1^- | j^{\mu} | 2^- \rangle \f$
    */

    template<typename T>
    struct kf1m2m : public KinFacABC<T>
    {
      kf1m2m(void); //! not implemented
      kf1m2m(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };

    /**
       @brief \f$ \langle 1^- | j^{\mu} | 2^+ \rangle \f$
    */

    template<typename T>
    struct kf1m2p : public KinFacABC<T>
    {
      kf1m2p(void); //! not implemented
      kf1m2p(kinIni &_ini) //! give ini information
      : ini(_ini)
      {    }

      kinIni ini;


        typename LinSysRetType<T>::type operator()(const KinKey &key);
    };


    // impl
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    typename LinSysRetType<T>::type kf0p0p<T>::operator()(const KinKey &key)
    {
      typename LinSysRetType<T>::type ret;
        int ncfg = key.src.E.size();
	
	POW2_ASSERT_DEBUG(key.sink.E.size() == ncfg);
        ret.reDim(static_cast<tensor::idx_t>(ncfg), 1);

	SourceKey source = key.src;
	SinkKey sink = key.sink; 
	
	double Esource = ENSEM::toDouble(ENSEM::mean(source.E));
	double Esink = ENSEM::toDouble(ENSEM::mean(source.E));
	double msource2 = (Esource*Esource
			   -source.mom[0]*source.mom[0] 
			   -source.mom[1]*source.mom[1] 
			   -source.mom[2]*source.mom[2]);

	double msink2 = (Esink*Esink
			 -sink.mom[0]*sink.mom[0]
			 -sink.mom[1]*sink.mom[1]
			 -sink.mom[2]*sink.mom[2]);

	itpp::Vec<T> dum(1);

	// sanity
	POW2_ASSERT(msource2 > 0.);
	POW2_ASSERT(msink2 > 0.);
	
	  if( fabs(msource2 - msink2) < ini.massCut) // same mass
	  {
	    tensor::Tensor<double,1> pmu_plus;
	    pmu_plus.create(std::vector<tensor::idx_t>(1,4));
	    EnsemReal E_plus = ENSEM::rescaleEnsemDown(source.E) 
	      + ENSEM::rescaleEnsemDown(sink.E);
	    pmu_plus[1] = source.mom[0] + sink.mom[0];
	    pmu_plus[2] = source.mom[1] + sink.mom[1];
	    pmu_plus[3] = source.mom[2] + sink.mom[2];

	    for(tensor::idx_t cfg = 0; cfg < ncfg; ++cfg)
	      {
		pmu_plus[0] = ENSEM::toDouble(E_plus.elem(cfg));
		
		for(tensor::idx_t lorentz = 0; lorentz < 4; ++lorentz)
		  {
		    dum[0] = pmu_plus[lorentz];
		    ret[lorentz][cfg] = dum;
		  }
	      }

	  }
	else                               // different mass
	  {
	    tensor::Tensor<double,1> pmu_plus, pmu_minus;
	    pmu_plus.create(std::vector<tensor::idx_t>(1,4));
	    pmu_minus = pmu_plus;
	    EnsemReal E_plus = ENSEM::rescaleEnsemDown(source.E) 
	      + ENSEM::rescaleEnsemDown(sink.E);
	    EnsemReal E_minus = ENSEM::rescaleEnsemDown(source.E) 
	      - ENSEM::rescaleEnsemDown(sink.E);

	    pmu_plus[1] = source.mom[0] + sink.mom[0];
	    pmu_plus[2] = source.mom[1] + sink.mom[1];
	    pmu_plus[3] = source.mom[2] + sink.mom[2];

	    pmu_minus[1] = source.mom[0] - sink.mom[0];
	    pmu_minus[2] = source.mom[1] - sink.mom[1];
	    pmu_minus[3] = source.mom[2] - sink.mom[2];
	    

	    for(tensor::idx_t cfg = 0; cfg < ncfg; ++cfg)						
	      {
		pmu_plus[0] = ENSEM::toDouble(E_plus.elem(cfg));
		pmu_minus[0] = ENSEM::toDouble(E_minus.elem(cfg));

		double Q2 = -(pmu_minus[0]*pmu_minus[0]
			      -pmu_minus[1]*pmu_minus[1]
			      -pmu_minus[2]*pmu_minus[2] 
			      -pmu_minus[3]*pmu_minus[3]);
		
		for(tensor::idx_t lorentz = 0; lorentz < 4; ++lorentz)
		  {
		    dum[0] =  ( pmu_plus[lorentz]*Q2/(msource2 - msink2)
				+ pmu_minus[lorentz]
				);

		    ret[lorentz][cfg] = dum;
		  }
	      }
	    
	  }

	// rescale the ensemble of kinematic factors back up
	ret.jackRescaleUp();

        return ret;
    }



}

#endif
