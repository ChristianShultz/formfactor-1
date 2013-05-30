#ifndef G_PARITY_WORLD_GENERATE_REDSTAR_XML_H
#define G_PARITY_WORLD_GENERATE_REDSTAR_XML_H


#include "g_parity_world.h"
#include "radmat/utils/obj_expr_t.h"

#include "hadron/hadron_npart_npt_corr.h"   
#include "ensem/ensem.h" 

namespace radmat
{

  // use one of these GParityState Thingys to generate a list of xml with subduction 
  // coefficients to get the operator 
  struct redstarSubducedOperator
  {
    typedef Hadron::KeyHadronNPartNPtCorr_t::NPoint_t redKey;  
    typedef ObjExpr_t<ENSEM::Complex,redKey> irrepOperator; 
    typedef ListObjExpr_t<ENSEM::Complex,redKey> listIrrepOperator;    

    void generate_subduction_list(const gParityWorld::GParityState &, const int t_slice); 
    listIrrepOperator subduced; 
    gParityWorld::GParityState orig;    
  };


  // now specialize to the GParityHelicityMatrixElement::State
  struct redstarSubducedState
  {
    typedef redstarSubducedOperator::redKey redKey;
    typedef redstarSubducedOperator::irrepOperator irrepOperator;
    typedef redstarSubducedOperator::listIrrepOperator listIrrepOperator; 

    void generate_subduction_list(const gParityWorld::GParityHelicityMatrixElement::State &); 

    listIrrepOperator subduced;
    gParityWorld::GParityHelicityMatrixElement::State orig; 
  };

  // specialize to the GParityInsertion -- collapse the old list over the flavor 
  // types down to a single list
  struct redstarSubducedPhoton
  {
    typedef redstarSubducedOperator::redKey redKey;
    typedef redstarSubducedOperator::irrepOperator irrepOperator;
    typedef redstarSubducedOperator::listIrrepOperator listIrrepOperator; 

    void generate_subduction_list(const gParityWorld::GParityInsertion::photon &, 
        const int t_slice); 

    listIrrepOperator subduced; 
    gParityWorld::GParityInsertion::photon orig; 
  };

  // generate the subduced insertion map
  struct redstarSubducedHelicityInsertion
  {
    typedef redstarSubducedOperator::redKey redKey;
    typedef redstarSubducedOperator::irrepOperator irrepOperator;
    typedef redstarSubducedOperator::listIrrepOperator listIrrepOperator; 
    typedef std::map<std::string,listIrrepOperator> map_t; 

    void generate_subduction_list(const gParityWorld::GParityInsertion &e); 

    gParityWorld::GParityInsertion orig; 
    map_t insertion_map; 
  };


  // join three listIrrepOperators into a list of 3pt xml keys with 
  // their coefficients .. ie multiply and sum 
  struct mergeSubducedOperators3pt
  {
    typedef redstarSubducedOperator::redKey redKey;
    typedef redstarSubducedOperator::irrepOperator irrepOperator;
    typedef redstarSubducedOperator::listIrrepOperator listIrrepOperator; 

    typedef Hadron::KeyHadronNPartNPtCorr_t NPointKey;  
    typedef ObjExpr_t<ENSEM::Complex,NPointKey> objNPointKey; 
    typedef ListObjExpr_t<ENSEM::Complex,NPointKey> listNPointKey;  

    void merge_subduction_lists(const listIrrepOperator &sink,
        const listIrrepOperator &insertion,
        const listIrrepOperator &source,
        const std::string &ensemble);

    listNPointKey subduced; 
    listIrrepOperator m_sink,m_insertion,m_source;
  };


  // take the internal representation of a matrix element and 
  // use it to generate redstar xml
  struct generateCircularRedstarXML
  {
    typedef redstarSubducedOperator::redKey redKey;
    typedef redstarSubducedOperator::irrepOperator irrepOperator;
    typedef redstarSubducedOperator::listIrrepOperator listIrrepOperator; 

    typedef mergeSubducedOperators3pt::NPointKey NPointKey;
    typedef mergeSubducedOperators3pt::objNPointKey objNPointKey; 
    typedef mergeSubducedOperators3pt::listNPointKey listNPointKey; 

    typedef redstarSubducedHelicityInsertion::map_t map_t; 

    typedef std::pair<bool,listNPointKey> data_t;

    generateCircularRedstarXML() 
    {
      time = data_t(false,listNPointKey());
      plus = data_t(false,listNPointKey());
      zero = data_t(false,listNPointKey());
      minus = data_t(false,listNPointKey());
    } 

    void fill(const std::string &id, const map_t& , data_t &, 
        const listIrrepOperator &sink, const listIrrepOperator &source, 
        const std::string &ensemble);

    void generate_subduction_list(const gParityWorld::GParityHelicityMatrixElement &); 

    // lorentz components in a helicity basis
    data_t time;
    data_t plus;
    data_t zero;
    data_t minus; 
    gParityWorld::GParityHelicityMatrixElement orig; 
  };

  // take the internal representation and turn it into the cartesian representation
  // via the intermediary Circular class above
  struct generateCartesianRedstarXML
  {
    typedef generateCircularRedstarXML::data_t data_t; 
    typedef mergeSubducedOperators3pt::NPointKey NPointKey;
    typedef mergeSubducedOperators3pt::objNPointKey objNPointKey; 
    typedef mergeSubducedOperators3pt::listNPointKey listNPointKey;

    generateCartesianRedstarXML(const gParityWorld::GParityHelicityMatrixElement &elem)
    { generate_subduction_list(elem); }

    void generate_subduction_list(const gParityWorld::GParityHelicityMatrixElement &); 

    data_t combine_duplicates(const data_t &);

    void run_g_parity_symmetry(void) {std::cout << __func__ << "implement me!!" << std::endl;}
    void run_cubic_symmetry(void) {std::cout << __func__ << "implement me !!" << std::endl;}

    // lorentz components in a cartesian basis
    data_t t,x,y,z;
    gParityWorld::GParityHelicityMatrixElement orig; 
  };





} // radmat




#endif /* G_PARITY_WORLD_GENERATE_REDSTAR_XML_H */
