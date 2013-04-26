#ifndef G_PARITY_WORLD_H
#define G_PARITY_WORLD_H

#include <string>
#include <iostream>
#include <vector>
#include <map>


#include "g_parity_world_xml.h"
#include "radmat/utils/obj_expr_t.h"

// this is the internal representation of a matrix element 


namespace radmat
{

  namespace gParityWorld
  {

    //! this is for internal storage -- represents a state (meson) 
    struct GParityState
    {
      GParityState() : name("un-init") {}

      int J;
      int H;
      bool parity;
      ADATXML::Array<int> mom;
      int twoI_z;
      std::string name;
      bool creation_op;
      bool smearedP;
      bool isProjected;
    };


    //! write it to a string for error
    std::string toString(const GParityState &);

    //! stream a ContinuumExprPrim
    std::ostream& operator<<(std::ostream& , const GParityState &);

    std::string stateFileName(const GParityState &s); 

    // internal representation of a set of insertions in the helicity basis (another set of mesons)
    struct GParityInsertion
    {
      typedef ListObjExpr_t<double,GParityState> photon;
      typedef ObjExpr_t<double,GParityState> photon_frag; 
      typedef std::map<std::string,photon> map_t; 
    
      GParityInsertion() {}

      void insert(const std::string &s, const photon &p)
      {
        if( insertion_map.find(s) == insertion_map.end())
          insertion_map.insert(map_t::value_type(s,p));
        else
        {
          insertion_map[s] = insertion_map[s] + p;
        }
      }

      void insert(const std::string &s, const photon_frag &p)
      {
        insert(s,photon(p)); 
      }

      map_t::const_iterator begin(void) const {return insertion_map.begin();}
      map_t::const_iterator end(void) const {return insertion_map.end();}

      void set_create(const bool &b) {creation_op = b;}
      bool create(void) const {return creation_op;}

      int t_slice; 
      bool creation_op;
      ADATXML::Array<int> mom;  // usefull to have sitting around
      map_t insertion_map; // keys t , m , 0 , p 
    };

    //! stream
    std::string toString(const GParityInsertion &);
    std::ostream& operator<<(std::ostream& , const GParityInsertion &);


    // internal representation of a helicity matrix element <meson (meson_list) meson>
    struct GParityHelicityMatrixElement
    {
      struct State
      {
        State(void) {}
        State(const GParityState &s , const int t) : state(s) , t_slice(t) {}
        GParityState state; 
        int t_slice; 
      };

      GParityHelicityMatrixElement() : ensemble("fred") {}

      State sink;
      GParityInsertion insertion;
      State source;

      std::string ensemble;  
    };

    //! State
  
    //! stream 
    std::string toString(const GParityHelicityMatrixElement::State &);
    std::ostream& operator<<(std::ostream&, const GParityHelicityMatrixElement::State &);
    std::string stateFileName(const GParityHelicityMatrixElement::State &);

    // GParityHelicityMatrixElement
    std::string toString(const GParityHelicityMatrixElement &); 
    std::ostream& operator<<(std::ostream & , const GParityHelicityMatrixElement &); 


    std::vector<GParityHelicityMatrixElement> 
      getGParityHelicityMatrixElementFromXML(const GParityContinuumMatElemXML &);

  } // gParityWorld


} // radmat



#endif /* G_PARITY_WORLD_H */
