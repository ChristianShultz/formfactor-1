#ifndef LLSQ_SOLUTION_H
#define LLSQ_SOLUTION_H 

#include "semble/semble_meta.h"
#include "radmat/construct_data/lattice_multi_data_tag.h"
#include "io/adat_xmlio.h"


namespace radmat
{
  
  template<typename T> 
    struct FormFacSolutions
    {
      FormFacSolutions(const SEMBLE::SembleMatrix<T> &ff, 
          const std::vector<LatticeMultiDataTag> &i)
        : FF_t(ff) , Ingredients(i) 
      { }

      void append_ingredients(const LatticeMultiDataTag &t)
      {
        Ingredients.push_back(t); 
      }

      SEMBLE::SembleMatrix<T> FF_t; 
      std::vector<LatticeMultiDataTag> Ingredients; 
    }; 

  template<typename T> 
    void 
    write(ADATIO::BinaryWriter &bin, const FormFacSolutions<T> &f)
   {
     int nr = f.FF_t.getN(); 
     write(bin,nr); 
      for(int ff = 0; ff < nr; ++ff)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type e; 
        SEMBLE::SembleVector<T> foo = f.FF_t.getRow(ff); 
        e.resize(foo.getB()); 
        e.resizeObs(foo.getN()); 

        for(int n = 0; n < foo.getN(); ++n)
          ENSEM::pokeObs(e,foo.getEnsemElement(n),n); 

        ENSEM::write(bin,e); 
      }
  
      int nt = f.Ingredients.size(); 
      write(bin,nt); 
      for(int i =0; i < nt; ++i)
        write(bin,f.Ingredients[i]); 
   } 

  template<typename T> 
    void 
    read(ADATIO::BinaryReader &bin, FormFacSolutions<T> &out)
    {
      FormFacSolutions<T> f; 
      int nr;
      read(bin,nr); 
      for(int ff = 0; ff < nr; ++ff)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type e; 
        SEMBLE::SembleVector<T> foo; 
        ENSEM::read(bin,e); 
        foo = e; 
        f.FF_t.append_row(foo); 
      }

      int nt; 
      read(bin,nt); 
      for(int i =0; i < nt; ++i)
        read(bin,f.Ingredients[i]); 

      out = f; 
    } 

}



#endif /* LLSQ_SOLUTION_H */
