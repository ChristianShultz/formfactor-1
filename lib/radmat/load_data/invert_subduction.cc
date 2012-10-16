/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 12-10-2012

 * Last Modified : Tue Oct 16 09:22:41 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/




#include "radmat/load_data/invert_subduction.h"
#include "radmat/utils/pow2assert.h"

#include "hadron/subduce_tables_factory.h"
#include "semble/semble_meta.h"

#include <sstream>
#include <vector>
#include <algorithm>
#include <list>
#include <complex>


namespace
{
  bool registered_subduction_tables = false;

  void doRegistration(void)
  {
    if(!!!registered_subduction_tables)
      POW2_ASSERT(Hadron::SubduceTableEnv::registerAll());

    registered_subduction_tables = true;
  }

  std::vector<std::string> getSubTableKeys(void)
  {
    doRegistration();
    return Hadron::TheSubduceTableFactory::Instance().keys();
  }

  bool match_string(const std::string &pattern, const std::string &text)
  {
    return std::search(text.begin(), text.end(), pattern.begin(), pattern.end()) != text.end();
  }

  std::vector<std::string> tokenize(const std::string &s, const std::string &delim, const bool keep_empty_tokens=false)
  {
    std::vector<std::string> tokens;

    if(delim.empty())
      return std::vector<std::string>(1,s);

    std::string::const_iterator front,back;
    front = s.begin();
    const std::string::const_iterator end = s.end();

    while(true)
    {

      back = std::search(front,end,delim.begin(),delim.end());
      std::string token(front,back);
      if(keep_empty_tokens || !!!token.empty())
        tokens.push_back(token);

      if(back == end)
        break;

      front = back + delim.size();
    }

    return tokens; 
  }



  std::vector<std::string> getSubTableBosonKeys(void)
  {
    std::vector<std::string> bosons;
    std::vector<std::string> all_keys(getSubTableKeys());
    std::vector<std::string>::const_iterator it;
    std::string dont_match("o"); // all of the fermion labels feature something like J1o2 for spin half
    for(it = all_keys.begin(); it != all_keys.end(); ++it)
      if(!!!match_string(dont_match,*it))
        bosons.push_back(*it);
    return bosons;
  }

  std::vector<std::string> getBosonKeysPattern(const std::string &pattern)
  {
    std::vector<std::string> boson_keys = getSubTableBosonKeys();
    std::vector<std::string> oct_keys;
    std::vector<std::string>::const_iterator it;
    for(it = boson_keys.begin(); it != boson_keys.end(); ++it)
      if(match_string(pattern,*it))
        oct_keys.push_back(*it);
    return oct_keys;
  }

  std::vector<std::string> getBosonKeys(const std::string &group)
  {
    if(group == "oct")
      return getBosonKeysPattern(std::string("J"));
    if(group == "D4")
      return getBosonKeysPattern(group);
    if(group == "D3")
      return getBosonKeysPattern(group);
    if(group == "D2") 
      return getBosonKeysPattern(group);
    if(group == "C4nm0") 
      return getBosonKeysPattern(group);
    if(group == "C4nnm")
      return getBosonKeysPattern(group);

    std::cerr << __func__ <<": ERROR: the group " << group 
      << " is not a supported/available symmetry group" << std::endl;
    exit(1);
  }

  bool inFlight(const std::string &group)
  {
    if(group == "oct")
      return false;
    if(group == "D4")
      return true;
    if(group == "D3")
      return true;
    if(group == "D2") 
      return true;
    if(group == "C4nm0") 
      return true;
    if(group == "C4nnm")
      return true;

    std::cerr << __func__ << ": ERROR: the group " << group 
      << " is not a supported/available symmetry group" << std::endl;
    exit(1);
  }

  //! move from (-j..0..j) to (1..(2j + 1)) indexing 
  int remapHelicity_1based(const radmat::ContinuumBosonExprPrimitive &expr)
  {
    return expr.H + expr.J + 1;
  }

  // get all of the irreps that the cont operator went into dependent on whatever group it was in
  std::vector<std::string> getIrreps(const radmat::ContinuumBosonExprPrimitive &expr)
  {
    const std::string delimiter1("->");
    const std::string delimiter2(",1");
    const std::vector<std::string> group_keys(getBosonKeys(expr.group));
    std::stringstream parse_id_stream; 

    if(inFlight(expr.group))
      parse_id_stream << "H" << abs(expr.H);
    else
      parse_id_stream << "J" << expr.J; 

    std::string parse_id = parse_id_stream.str();

    std::vector<std::string> irreps;  
    std::vector<std::string>::const_iterator it;

    for(it = group_keys.begin(); it != group_keys.end(); ++it)
    {
      if(match_string(parse_id,*it))
      {
        //   std::cout << __func__ << ": matched " << *it << std::endl;
        std::vector<std::string> tokens = tokenize(*it,delimiter1);
        irreps.push_back( tokenize(tokens[1],delimiter2)[0] );
      }
    }
    return irreps;
  }


  radmat::LatticeIrrepExpr_t makeLatticeIrrepExpr(const ENSEM::Complex &coefficient,
      const std::string &group, 
      const std::string &irrep,
      const int row)
  {
    return radmat::LatticeIrrepExpr_t(coefficient, radmat::LatticeExprPrimitive(group,irrep,row));
  }




} // namespace anonymous







namespace radmat
{










  //! write this into a string for the factory
  std::string toString(const ContinuumBosonExprPrimitive &expr)
  {
    std::stringstream ss;

    ss << "_J_" << expr.J << "_H_" << expr.H << "_parity_";
    if(expr.parity)
      ss << "p";
    else
      ss << "m";

    ss << "_" << expr.group;

    return ss.str();
  }

  //! stream a ContinuumExprPrimitive
  std::ostream& operator<<(std::ostream &os, const ContinuumBosonExprPrimitive &expr)
  {
    os << toString(expr);
    return os;
  }

  //! write into a string
  std::string toString(const LatticeExprPrimitive &p)
  {
    std::stringstream ss;
    ss << "_" << p.group << "_" << p.irrep << "_row_" << p.row;
    return ss.str();
  }

  //! stream a LatticeExprPrimitive
  std::ostream& operator<<(std::ostream &os, const LatticeExprPrimitive &p)
  {
    os << toString(p);
    return os;
  }



  ListLatticeIrrepExpr_t invertSubduction(const ContinuumBosonExprPrimitive & expr)
  {
    std::vector<std::string> irreps = getIrreps(expr);
    std::vector<std::string>::const_iterator it;
    ListLatticeIrrepExpr_t lattice_expr;

    for(it = irreps .begin(); it != irreps.end(); ++it)
    {
      lattice_expr.push_back(makeLatticeIrrepExpr(SEMBLE::toScalar(std::complex<double>(1,0)), expr.group, *it ,1));
    }
    return lattice_expr; 
  }




} // namespace radmat
