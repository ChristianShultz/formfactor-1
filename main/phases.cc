/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : phases.cc

* Purpose :

* Creation Date : 18-01-2013

* Last Modified : Mon Jan 21 13:20:50 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "hadron/clebsch.h"
#include "hadron/irrep_util.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/obj_expr_t.h"
#include "radmat/utils/pow2assert.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "itpp/itbase.h"

typedef radmat::ObjExpr_t<std::complex<double> , std::string> J_t; 
typedef radmat::ListObjExpr_t<std::complex<double> , std::string> listJ_t; 
const std::complex<double> one(1.,0.);
const std::complex<double> cplx_i(0.,1.);
const double root2(sqrt(2.)); 
const double root2inv(1./root2); 


template<typename C, typename O>
  typename itpp::Vec< radmat::ListObjExpr_t<C,O> > 
operator*(const itpp::Mat<C> &m, 
    const typename itpp::Vec< radmat::ListObjExpr_t<C,O> > &v)
{
  typename itpp::Vec< radmat::ListObjExpr_t<C,O> > ret(v.size());

  POW2_ASSERT(v.size() == m.cols()); 
  radmat::ListObjExpr_t<C,O> tmp; 

  for(int row = 0; row < m.rows(); ++row)
    for(int col = 0; col < v.size(); ++col)
    {
      tmp = ret[row] + m(row,col)*v(col);
      ret[row] = tmp; 
    }

  return ret; 
}


itpp::Vec<std::complex<double> > eps3_z(const std::string &qn)
{
  itpp::Vec<std::complex<double> > ret(3);
  ret.zeros();

  if(qn == "p")
  {
    ret[0] = -root2inv;
    ret[1] = -cplx_i/root2; 
  }
  else if(qn == "m")
  {
    ret[0] = root2inv;
    ret[1] = -cplx_i/root2;
  }
  else if(qn == "0")
  {
    ret[2] = 1;
  }
  else
  {
    std::cerr << "qn(" << qn << ") is ill defined, use p,m,0" << std::endl;
    exit(1);
  }

  return ret;
}


// columns are indexed -m to m
// rows go -lambda to lambda

itpp::Mat<std::complex<double> > Wigner_D(const ADATXML::Array<int> mom, const int J)
{
  Hadron::CubicCanonicalRotation_t rot = Hadron::cubicCanonicalRotation(mom);
  itpp::Mat<std::complex<double> > ret(2*J + 1, 2*J+1); 
  
  const int tJ = 2*J;

  for(int lambda = -J; lambda < J +1; ++lambda)
    for(int m = -J; m < J +1; ++m)
    {
      ret(lambda+J,m+J) = SEMBLE::toScalar(Hadron::Wigner_D(tJ,2*m,2*lambda,rot.alpha,rot.beta,rot.gamma));
    }

  return ret; 
}




itpp::Mat<std::complex<double> > Wigner_D_jz_p(const ADATXML::Array<int> mom, const int J)
{
  Hadron::CubicCanonicalRotation_t rot = Hadron::cubicCanonicalRotation(mom);
  itpp::Mat<std::complex<double> > ret(2*J + 1, 2*J+1); 
  
  const int tJ = 2*J;

  for(int mp = -J; mp < J +1; ++mp)
    for(int m = -J; m < J +1; ++m)
    {
      ret(mp+J,m+J) = SEMBLE::toScalar(Hadron::Wigner_D(tJ,2*mp,2*m,rot.alpha,rot.beta,rot.gamma));
    }

  return ret; 
}






std::string doPrint(const listJ_t &l)
{
  std::map<std::string,std::complex<double> > mp;
  std::map<std::string,std::complex<double> >::iterator mpit;
  listJ_t::const_iterator it;

  for(it = l.begin(); it != l.end(); ++it)
  {
    mpit = mp.find(it->m_obj);
    if(mpit != mp.end())
    {
      mpit->second  += it->m_coeff; 
    }
    else
    {
      mp.insert(std::pair<std::string,std::complex<double> >(it->m_obj,it->m_coeff));
    }
  }

  std::stringstream ss; 
  for(mpit = mp.begin(); mpit != mp.end(); ++mpit)
    if(!!! (itpp::round_to_zero(mpit->second,0.0001) == 0.))
      ss << mpit->second  << " " << mpit->first << std::endl;
    else
      ss << "0+0i " << mpit->first << std::endl;

  return ss.str(); 
}



int main(void)
{

  ADATXML::Array<int> mom;
  mom.resize(3);
  mom[0] = -1;
  mom[1] = 0;
  mom[2] = 0;
  radmat::ListObjExpr_t<ENSEM::Complex,std::string> m_list; 
  Hadron::CubicCanonicalRotation_t rot = Hadron::cubicCanonicalRotation(mom);

  std::cout << "alpha " << rot.alpha << " beta " << rot.beta << " gamma " << rot.gamma << std::endl;


  listJ_t x,y,z;
  x.push_back(J_t(one,"x"));
  y.push_back(J_t(one,"y"));
  z.push_back(J_t(one,"z"));

  itpp::Vec<listJ_t> j_i(3), j_m(3),j_Deps(3),j_lambda,j_inv;
  j_i[0] = x;
  j_i[1] = y;
  j_i[2] = z;


  itpp::Mat<std::complex<double> > epsz, D,Deps; 
  epsz.set_size(3,3); 

  epsz.set_row(0,eps3_z("m"));
  epsz.set_row(1,eps3_z("0"));
  epsz.set_row(2,eps3_z("p"));

  j_m = epsz*j_i;


  for(int m = 0; m < j_m.size(); ++m)
  {
    std::cout << "m = " << m - 1 << std::endl;
    std::cout << j_m[m] << std::endl;
  }


  D = itpp::round_to_zero(Wigner_D(mom,1), 0.0001);


  j_lambda = D*j_m;

  Deps = D*epsz;


  for(int m = 0; m < j_m.size(); ++m)
  {
    std::cout << "lambda = " << m -1 << std::endl;
    std::cout << doPrint(j_lambda[m]) << std::endl;
  }


  j_Deps = Deps*j_i; 
  for(int m = 0; m < j_m.size(); ++m)
  {
    std::cout << "lambda = " << m -1 << std::endl;
    std::cout << doPrint(j_Deps[m]) << std::endl;
  }

  std::cout << Deps << std::endl; 
  
  itpp::Mat<std::complex<double> > invert = itpp::round_to_zero( itpp::inv(Deps), 0.00001); 

  j_inv = invert*j_Deps;

   for(int m = 0; m < j_m.size(); ++m)
  {
    std::cout << "lambda = " << m -1 << std::endl;
    std::cout << doPrint(j_inv[m]) << std::endl;
  }



  return 0;
}



#if 0


double rd(const double &inp)
{
  if(fabs(inp) < 0.0000001)
    return 0.;
  else return inp; 
}

std::complex<double> rd(const std::complex<double> &inp)
{
  return std::complex<double>(rd(inp.real()),rd(inp.imag()));
}


ADATXML::Array<int> mom;
mom.resize(3);
mom[0] = -1;
mom[1] = 0;
mom[2] = 0;
radmat::ListObjExpr_t<ENSEM::Complex,std::string> m_list; 
Hadron::CubicCanonicalRotation_t rot = Hadron::cubicCanonicalRotation(mom);

std::cout << "D_m_l" << std::endl;
for(int m = -1; m < 2; ++m)
{
  std::cout << "m = " << m << "  ";
  for(int lambda = -1; lambda < 2; ++lambda)
    std::cout << std::setw(20) << rd(SEMBLE::toScalar(Hadron::Wigner_D(2,2*m,2*lambda,rot.alpha,rot.beta,rot.gamma))) << "  ";
  std::cout << std::endl;
}

std::cout << std::endl; 
std::cout << "D_m_l(R^-1)" << std::endl;
for(int m = -1; m < 2; ++m)
{
  std::cout << "m = " << m << "  ";
  for(int lambda = -1; lambda < 2; ++lambda)
    std::cout << std::setw(20) << rd(SEMBLE::toScalar(Hadron::Wigner_D(2,2*m,2*lambda,-rot.gamma,-rot.beta,-rot.alpha))) << "  ";
  std::cout << std::endl;
}


std::cout << std::endl;
std::cout << "D_l_m" << std::endl;
for(int lambda = -1; lambda < 2; ++lambda)
{
  std::cout << "l = " << lambda << "  ";
  for(int m = -1; m < 2; ++m)
    std::cout << std::setw(20) << rd(SEMBLE::toScalar(Hadron::Wigner_D(2,2*m,2*lambda,rot.alpha,rot.beta,rot.gamma))) << "  ";
  std::cout << std::endl;
}




std::map<int,std::string> hel;
std::map<int,std::string>::const_iterator it;
int twoJ = 2; 
int target = 1; 

hel[2] = "+";
hel[0] = "0";
hel[-2] = "-";


for(it = hel.begin(); it != hel.end(); ++it)
{
  ENSEM::Complex ph = Hadron::Wigner_D(twoJ,2*target,it->first,rot.alpha,rot.beta,rot.gamma);
  m_list.push_back(radmat::ObjExpr_t<ENSEM::Complex,std::string>(ph,it->second));
}


radmat::ListObjExpr_t<ENSEM::Complex,std::string>::const_iterator list_it;
for(list_it = m_list.begin(); list_it != m_list.end(); ++list_it)
std::cout << SEMBLE::toScalar(list_it->m_coeff) << "  *  " << list_it->m_obj << std::endl;

#endif

