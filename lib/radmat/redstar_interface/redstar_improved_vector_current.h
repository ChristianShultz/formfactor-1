#ifndef REDSTAR_IMPROVED_VECTOR_CURRENT_H
#define REDSTAR_IMPROVED_VECTOR_CURRENT_H 

#include "redstar_abstract_block_base.h"
#include "radmat/utils/stringify.h"
#include "io/adat_xmlio.h"
#include <utility>

// follows arXiv:hep-lat/0103026v1 (2001) Harada, Kronfeld
//        Matsufuru, Nakajima, Onogi

namespace radmat
{
  //
  //
  // this is a representation of an improved photon fragment
  struct RedstarImprovedVectorCurrentPFrag
  {
    std::string op_name;          // this is the photon stem, something like rho_rhoxD0_J0__J1 
    double op_coeff;              // this is the coefficient it gets
    std::string i_name;           // this is the stem of the improvement -- "no_improvement" is a breakout flag 
    double i_coeff_r;             // this is a real coefficient that it gets
    double i_coeff_i;             // this is a imag coefficient that it gets
    double m_qxa_t;               // temporal spacing times m_q, (-0.0840)
  };

  //
  //
  // improvement parameters
  //    -- everyone gets a copy since it is tiny
  //       relative to the size of a single correlator
  struct RedstarImprovedVectorCurrentParPack
  {
    double        xi;           // anisotropy
    int           L_s;          // spatial lattice spacing
    double        d_1;          // the parameter from paper
    double        gamma_f;      // bare fermion anisotropy
    double        mu_tilde_t;   // temporal tadpole factor
    double        msnk;         // rest mass of sink
    double        msrc;         // rest mass of source
  };

  //
  //
  // this is the input structure for improved vector currents
  // 
  //
  struct RedstarImprovedVectorCurrentInput;
  REGISTER_STRINGIFY_TYPE(RedstarImprovedVectorCurrentInput);

  struct RedstarImprovedVectorCurrentInput
    : public AbsRedstarInput_t
  {
    typedef RedstarImprovedVectorCurrentPFrag pfrag; 
    typedef std::vector<pfrag>::const_iterator const_iterator; 
    typedef RedstarImprovedVectorCurrentParPack ImprovementPack; 

    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarImprovedVectorCurrentInput>()); 
    }

    virtual std::string write(void) const; 

    virtual rHandle<AbsRedstarInput_t> clone(void) const 
    {
      return rHandle<AbsRedstarInput_t>(
          new RedstarImprovedVectorCurrentInput(*this) ); 
    }

    const_iterator begin(void) const {return photons.begin();}
    const_iterator end(void) const {return photons.end();}

    int lorentz;                            // 1,2,3,4 -- euclidean time is 4
    ADATXML::Array<int> mom;                // insertion momentum 
    ADATXML::Array<int> psrc;               // sink momentum 
    ADATXML::Array<int> psnk;               // source momentum 
    std::vector<pfrag> photons;             // photons and improvement term
    bool creation_op;                       // momentum conventions
    bool smearedP;                          // is it smeared
    int t_slice;                            // where does it live 
    ImprovementPack ipack;                  // a bunch of pure numbers
  }; 

  //
  //
  //  the improved vector current functor
  //
  //
  struct RedstarImprovedVectorCurrentBlock;
  REGISTER_STRINGIFY_TYPE(RedstarImprovedVectorCurrentBlock);

  struct RedstarImprovedVectorCurrentBlock
    : public AbsRedstarBlock_t
  {

    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarImprovedVectorCurrentBlock>()); 
    }

    virtual EnsemRedstarBlock
      operator()(const AbsRedstarInput_t *ptr2derived) const ; 
  }; 


  //
  // an actual insertion (temporal or spatial)
  //
  struct RedstarImprovedVectorCurrentInsertion
  {
    bool active; 
    bool creation_op; 
    bool smearedP; 
    ADATXML::Array<RedstarImprovedVectorCurrentPFrag> photons; 
  };

  //
  //
  //  the outside world xml interface
  //
  //
  struct RedstarImprovedVectorCurrentXML;
  REGISTER_STRINGIFY_TYPE(RedstarImprovedVectorCurrentXML); 

  struct RedstarImprovedVectorCurrentXML
    : public AbsRedstarXMLInterface_t
  {
    typedef RedstarImprovedVectorCurrentParPack ImprovementPack; 
    typedef RedstarImprovedVectorCurrentInsertion Insertion; 

    virtual void read(ADATXML::XMLReader &xml, const std::string &path); 
    virtual std::string type(void) const 
    {
      return std::string(Stringify<RedstarImprovedVectorCurrentXML>()); 
    }
    virtual std::string write(void) const; 
    virtual void write(ADATXML::XMLWriter &xml, const std::string &pth) const; 
    virtual int timeslice(void) const {return t_slice;}

    int pmin;
    int pmax;
    int t_slice;
    ImprovementPack ipack;
    Insertion time;
    Insertion space;
  }; 

  // xml readers
  void read(ADATXML::XMLReader &xml,
      const std::string &path,
      RedstarImprovedVectorCurrentPFrag &);

  void read(ADATXML::XMLReader &xml,
      const std::string &path,
      RedstarImprovedVectorCurrentInsertion &);

  void read(ADATXML::XMLReader &xml,
      const std::string &path,
      RedstarImprovedVectorCurrentParPack &);

  // xml writers
  void write(ADATXML::XMLWriter &xml,
      const std::string &path,
      const RedstarImprovedVectorCurrentPFrag &);

  void write(ADATXML::XMLWriter &xml,
      const std::string &path,
      const RedstarImprovedVectorCurrentInsertion &);

  void write(ADATXML::XMLWriter &xml,
      const std::string &path,
      const RedstarImprovedVectorCurrentParPack &);
} // radmat 


#endif /* REDSTAR_IMPROVED_VECTOR_CURRENT_H */
