/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 15-10-2012

 * Last Modified : Tue Oct 16 09:21:57 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "../headers/tester.h"
#include "../headers/test_common.h"
#include "radmat/load_data/invert_subduction.h"

#include <string>



namespace radmat
{

  tester test_invert_subduction(void)
  {
    tester m_test(__func__);

    for(int i = 0; i < 5; ++i)
    {
      ContinuumBosonExprPrimitive cont_expr(i,false,i,std::string("oct"));

      ListLatticeIrrepExpr_t  inverse_subducer = invertSubduction(cont_expr);

      std::cout << "J = " << i << "\n" << inverse_subducer << std::endl;
    }


    return m_test;
  }


}

