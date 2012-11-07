/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : simple_world.cc

 * Purpose :

 * Creation Date : 19-10-2012

 * Last Modified : Wed Oct 24 15:05:15 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "simple_world.h"
#include <sstream>


using namespace ENSEM;
using namespace ADATXML;
using namespace Hadron;
using namespace radmat::simpleWorld;

namespace radmat
{

  namespace simpleWorld
  {


    namespace
    {


      template<typename T>
        void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
        {
          if(ptop.count(path) > 0)
            read(ptop,path,place);
          else
          {
            std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
              << ", path was empty, exiting" << std::endl;
            exit(1);
          }
        }

    } // namespace anonomyous 




    //! write it to a string for error
    std::string toString(const ContinuumStatePrimitive &op)
    {
      std::stringstream ss;
      ss << "J = " << op.J << " H = " << op.H << " parity " << op.parity 
        << " mom = " << op.mom[0] << " " << op.mom[1] << " " << op.mom[2] 
        << " twoI_z " << op.twoI_z << " name " << op.name << " create " 
        << op.creation_op << " smear " << op.smearedP;
      return ss.str();
    }

    //! stream a ContinuumExprPrim
    std::ostream& operator<<(std::ostream& o , const ContinuumStatePrimitive &op)
    {
      o << toString(op);
      return o;
    }

    //! write to a string for error
    std::string toString(const ContinuumStateXML &op)
    {
      std::stringstream ss;
      ss << "J = " << op.J << " H = {";
      for(int i = 0; i < op.H.size(); ++i)
        ss << op.H[i] << " ";
      ss << "} parity = " << op.parity << " mom = { ";
      for(int i = 0; i < op.mom.size(); ++i)
        ss <<  "(" <<op.mom[i][0] << " " << op.mom[i][1] << " " << op.mom[i][2] << ") ";
      ss << "} twoI_z = " << op.twoI_z << " op_stem " << op.op_stem << " create = "
        << op.creation_op << " smeared = " << op.smearedP;
      return ss.str();
    }

    //! stream this thing
    std::ostream& operator<<(std::ostream& o , const ContinuumStateXML &op)
    {
      o << toString(op);
      return o;
    }

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumStateXML &op)
    {   
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"J",op.J,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"H",op.H,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"parity",op.parity,__PRETTY_FUNCTION__);      
      doXMLRead(ptop,"mom",op.mom,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"twoI_z",op.twoI_z,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"op_stem",op.op_stem,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"creation_op",op.creation_op,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"smearedP",op.smearedP,__PRETTY_FUNCTION__);
    } 

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumStateXML &op)
    {
      ADATXML::push(xml,path);
      write(xml,"J",op.J);
      write(xml,"H",op.H);
      write(xml,"parity",op.parity);
      write(xml,"mom",op.mom);
      write(xml,"twoI_z",op.twoI_z);
      write(xml,"op_stem",op.op_stem);
      write(xml,"creation_op",op.creation_op);
      write(xml,"smearedP",op.smearedP);
      ADATXML::pop(xml);
    }



    //! write state to a string
    std::string toString(const ContinuumMatElem::State &s)
    {
      std::stringstream ss;
      ss << "state = " << s.state << " t_slice = " << s.t_slice;
      return ss.str();
    }

    std::ostream& operator<<(std::ostream& o, const ContinuumMatElem::State &s)
    {
      o << toString(s);
      return o;
    }



    //! write it to a string for error
    std::string toString(const ContinuumMatElem & s)
    {
      std::stringstream ss;
      ss << " source = " << s.source << "\n";
      for(int ins = 0; ins < s.insertion.size(); ++ins)
        ss << "ins_" << ins << " = " << s.insertion[ins] << "\n";
      ss << "sink = " << s.sink;
      return ss.str();
    }

    //! stream a ContinuumMatElem
    std::ostream& operator<<(std::ostream& o, const ContinuumMatElem & s)
    {
      o << toString(s);
      return o;
    }



    //! write it to a string
    std::string toString(const ContinuumLorentzMatElem &s)
    {
      std::stringstream ss;
      ss << "source = " << s.source << "\n";
      for(int i = 0; i < 4; ++i)
        ss << " mu_" << i << " = " << s.lorentz[i] << "\n";
      ss << "sink = " << s.sink;
      return ss.str();
    }

    //! stream it
    std::ostream& operator<<(std::ostream& o, const ContinuumLorentzMatElem &s)
    {
      o << toString(s);
      return o;
    }


    //! write state to a string
    std::string toString(const ContinuumLorentzMatElem::State &s)
    {
      std::stringstream ss;
      ss << "state = " << s.state << " t_slice = " << s.t_slice;
      return ss.str(); 
    }

    std::ostream& operator<<(std::ostream& o , const ContinuumLorentzMatElem::State &s)
    {
      o << toString(s);
      return o;
    }


    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumLorentzMatElem::State &op)
    {   
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"t_slice",op.t_slice,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"state",op.state,__PRETTY_FUNCTION__);

    }

    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumLorentzMatElem::State &op)
    { 
      ADATXML::push(xml,path);
      write(xml,"state",op.state);
      write(xml,"t_slice",op.t_slice);
      ADATXML::pop(xml);
    }



    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumLorentzMatElem &op)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"source",op.source,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"lorentz",op.lorentz,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"sink",op.sink,__PRETTY_FUNCTION__);
    }

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumLorentzMatElem &op)
    {
      ADATXML::push(xml,path);
      write(xml,"source",op.source);
      write(xml,"lorentz",op.lorentz);
      write(xml,"sink",op.sink);
      ADATXML::pop(xml);
    }


  } // namespace simpleWorld

} // namespace radmat


