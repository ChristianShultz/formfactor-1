#ifndef G_PARITY_WORLD_XML_H
#define G_PARITY_WORLD_XML_H


#include "io/adat_xmlio.h"

#include <string>
#include <iostream>


// these are the external xml interface that allows us to specify a matrix element


namespace radmat
{

  namespace gParityWorld
  {


    //! this is for an XML interface 
    struct GParityContinuumStateXML
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
      bool isProjected;  
    };

    //! write to a string for error
    std::string toString(const GParityContinuumStateXML &);

    //! stream this thing
    std::ostream& operator<<(std::ostream& , const GParityContinuumStateXML &);

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &, GParityContinuumStateXML &);

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityContinuumStateXML &);


    // external representation of an insertion is a list of "photon bits"
    // and a coefficient
    struct GParityInsertionXML
    { 
      struct GParityListElement
      {
        std::string op_stem;
        double charge_coefficient;
        int twoI_z;
        std::string op_manip; // useful to split off improvement terms
      };

      struct Insertion
      {
        bool active;           // do we want to do time/space/ mix both of them together? 
        int J;
        ADATXML::Array<int> H; // take note of circular basis!!!        
        bool parity;
        bool creation_op;
        bool smearedP; 

        ADATXML::Array<GParityListElement> photons; 
      };      

      // we may be missing elementals 
      int pmax;
  
      // may only want to look at one thing
      int pmin;

      int t_slice; 

      //  one for spatial currents, one for temporal
      Insertion time;
      Insertion space; 
    };

    // GParityListElement

    // stream
    std::string toString(const GParityInsertionXML::GParityListElement &);   
    std::ostream& operator<<(std::ostream &o, const  GParityInsertionXML::GParityListElement &);

    //!  xml reader/writer
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityInsertionXML::GParityListElement &);
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityInsertionXML::GParityListElement &);

    // Insertion

    // stream
    std::string toString(const GParityInsertionXML::Insertion &);
    std::ostream& operator<<(std::ostream &o, const  GParityInsertionXML::Insertion &);

    //!  xml reader/writer
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityInsertionXML::Insertion &);
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityInsertionXML::Insertion &);

    // GParityInsertionXML

    // stream
    std::string toString(const GParityInsertionXML &);
    std::ostream& operator<<(std::ostream &o, const  GParityInsertionXML &);

    //!  xml reader/writer
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityInsertionXML &);
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityInsertionXML &);


    struct GParityContinuumMatElemXML
    {
      struct State
      {
        GParityContinuumStateXML state;
        int t_slice;
      };

      State sink;
      GParityInsertionXML insertion; 
      State source;


      std::string ensemble; 
    };


    //! write state to a string
    std::string toString(const GParityContinuumMatElemXML::State &);

    //! stream a state
    std::ostream& operator<<(std::ostream& , const GParityContinuumMatElemXML::State &);

    //! state xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityContinuumMatElemXML::State &);

    //! state xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityContinuumMatElemXML::State &);

    //! write it to a string
    std::string toString(const GParityContinuumMatElemXML &);

    //! stream it
    std::ostream& operator<<(std::ostream&, const GParityContinuumMatElemXML &);

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityContinuumMatElemXML &);

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityContinuumMatElemXML &);


  } // gParityWorld


} // radmat






#endif /* G_PARITY_WORLD_XML_H */
