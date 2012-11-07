#ifndef SIMPLE_WORLD_H
#define SIMPLE_WORLD_H


#include "hadron/hadron_npart_irrep.h"
#include "io/adat_xmlio.h"

#include <string>
#include <iostream>

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


    //------------------------------------------------------------------------------------------


    //! this is for an XML interface 
    struct ContinuumStateXML
    {
      int J;
      ADATXML::Array<int> H; 
      bool parity;
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


    //! internal representation of the states we're sandwiching an op between
    struct ContinuumMatElem 
    {

      struct State
      {
        ContinuumStatePrimitive state;
        int t_slice;
      };

      State source;
      ADATXML::Array<State> insertion;
      State sink;

      std::string ensemble;
    };

    //! write state to a string
    std::string toString(const ContinuumMatElem::State &);

    //! stream a state
    std::ostream& operator<<(std::ostream &, const ContinuumMatElem::State &);

    //! write it to a string for error
    std::string toString(const ContinuumMatElem &);

    //! stream a ContinuumMatElem
    std::ostream& operator<<(std::ostream& , const ContinuumMatElem &);


    //----------------------------------------------------------------------------------------


    struct ContinuumLorentzMatElem
    {

      struct State
      {
        ContinuumStateXML state;
        int t_slice;
      };

      State source;
      ADATXML::Array<State> lorentz;  // indexed by lorentz value
      State sink;

    }; 

    //! write state to a string
    std::string toString(const ContinuumLorentzMatElem::State &);

    //! stream a state
    std::ostream& operator<<(std::ostream& , const ContinuumLorentzMatElem::State &);

    //! state xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumLorentzMatElem::State &);

    //! state xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, ContinuumLorentzMatElem::State &);

    //! write it to a string
    std::string toString(const ContinuumLorentzMatElem &);

    //! stream it
    std::ostream& operator<<(std::ostream&, const ContinuumLorentzMatElem &);

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumLorentzMatElem &);

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumLorentzMatElem &);


  } // namespace simpleWorld



} // namespace radmat






#endif /* SIMPLE_WORLD_H */
