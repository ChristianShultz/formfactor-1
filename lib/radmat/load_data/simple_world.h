#ifndef SIMPLE_WORLD_H
#define SIMPLE_WORLD_H


#include "hadron/hadron_npart_irrep.h"
#include "io/adat_xmlio.h"

#include "radmat/load_data/invert_subduction.h"

#include <string>
#include <iostream>
#include <vector>
#include <map>


namespace radmat
{


  namespace simpleWorld
  {

    //! this is for internal storage
    struct ContinuumStatePrimitive
    {
      int J;
      int H;
      bool parity;
      ADATXML::Array<int> mom;
      int twoI_z;
      std::string name;
      bool creation_op;
      bool smearedP;
    };


    //! write it to a string for error
    std::string toString(const ContinuumStatePrimitive &);

    //! stream a ContinuumExprPrim
    std::ostream& operator<<(std::ostream& , const ContinuumStatePrimitive &);

    std::string stateFileName(const ContinuumStatePrimitive &s); 


    //------------------------------------------------------------------------------------------


    //! this is for an XML interface 
    struct ContinuumStateXML
    {
      int J;
      ADATXML::Array<int> H; 
      bool parity;
      bool fill_star;
      ADATXML::Array<ADATXML::Array<int> > mom;
      int twoI_z;
      std::string op_stem; 
      bool creation_op;
      bool smearedP; 
    };

    //! write to a string for error
    std::string toString(const ContinuumStateXML &);

    //! stream this thing
    std::ostream& operator<<(std::ostream& , const ContinuumStateXML &);

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &, ContinuumStateXML &);

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumStateXML &);

    //-------------------------------------------------------------------------------------------   



    /// here

    // internal representation of the guys 
    struct ContinuumInsertion
    {
      typedef ContinuumStatePrimitive op_insertion;   
   
      int t_slice; 

      void insert(const std::string &s, const op_insertion &o)
      {
        insertion_map.insert(std::map<std::string,op_insertion>::value_type(s,o));
      }

      std::map<std::string,op_insertion> insertion_map;  // keys are "t" , "p", "m", "0" .. not found == don't use
    };


    std::string toString(const ContinuumInsertion &);
    std::ostream& operator<<(std::ostream &, const ContinuumInsertion &); 

    //-------------------------------------------------------------------------------------------   


    //// here


    struct ContinuumInsertionXML
    {
      struct Insertion
      {
        int J;
        ADATXML::Array<int> H; // take note of circular basis!!!
        bool parity;
        int twoI_z;
        std::string op_stem; 
        bool creation_op;
        bool smearedP; 
      };

      int pmax;    // we currently dont have the meson elementals for 
      int t_slice; 
      Insertion time; 
      Insertion space; 
    };

    //! write state to a string
    std::string toString(const ContinuumInsertionXML::Insertion &);

    //! stream a state
    std::ostream& operator<<(std::ostream& , const ContinuumInsertionXML::Insertion &);

    //! state xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumInsertionXML::Insertion &);

    //! state xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, ContinuumInsertionXML::Insertion &);



    //! write to a string for error
    std::string toString(const ContinuumInsertionXML &);

    //! stream this thing
    std::ostream& operator<<(std::ostream& , const ContinuumInsertionXML &);

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &, ContinuumInsertionXML &);

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumInsertionXML &);

    //-------------------------------------------------------------------------------------------   


    /// here

    //! internal representation of the matrix element we're sandwiching an op between
    struct ContinuumMatElem 
    {
      struct State
      {
        State(void) {}
        State(const ContinuumStatePrimitive &s, const int t) : state(s), t_slice(t) {}

        ContinuumStatePrimitive state;
        int t_slice;
      };

      State source;
      ContinuumInsertion insertion; 
      State sink;

      std::string ensemble;
    };

    //! write state to a string
    std::string toString(const ContinuumMatElem::State &);

    //! stream a state
    std::ostream& operator<<(std::ostream &, const ContinuumMatElem::State &);

    std::string stateFileName(const ContinuumMatElem::State &); 

    //! write it to a string for error
    std::string toString(const ContinuumMatElem &);

    //! stream a ContinuumMatElem
    std::ostream& operator<<(std::ostream& , const ContinuumMatElem &);


    //----------------------------------------------------------------------------------------

    /// here

    struct ContinuumMatElemXML
    {

      struct State
      {
        ContinuumStateXML state;
        int t_slice;
      };

      State source;
      ContinuumInsertionXML insertion; 
      State sink;

      std::string ensemble; 
    }; 

    //! write state to a string
    std::string toString(const ContinuumMatElemXML::State &);

    //! stream a state
    std::ostream& operator<<(std::ostream& , const ContinuumMatElemXML::State &);

    //! state xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumMatElemXML::State &);

    //! state xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, ContinuumMatElemXML::State &);

    //! write it to a string
    std::string toString(const ContinuumMatElemXML &);

    //! stream it
    std::ostream& operator<<(std::ostream&, const ContinuumMatElemXML &);

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumMatElemXML &);

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumMatElemXML &);



    //! transform from xml to something we can loop on 
    std::vector<ContinuumMatElem> getContinuumMatElemFromXML(const ContinuumMatElemXML &);


  } // namespace simpleWorld



} // namespace radmat






#endif /* SIMPLE_WORLD_H */
