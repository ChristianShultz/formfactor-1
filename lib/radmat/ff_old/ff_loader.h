#ifndef FF_LOADER_H_H_GUARD
#define FF_LOADER_H_H_GUARD

#include "ff_base.h"
#include "formfac/fake.h"
#include "radmat/utils/pow2assert.h"
#include <string>


/**
   @file ff_loader.h
   @brief load the set of functors to create the linear system based on the mat_elem type
 */

namespace radmat
{
    struct ff_loader
    {
        typedef std::list<FFBlockBase_t *> ff_list;

        ff_loader(void) ; // hidden
        ff_loader(const std::string mat_elem_id)
            : elem_id(mat_elem_id)
        {  }

        ff_list genList(void)
        {
            POW2_ASSERT(elem_id == "fake");
            ff_list list;
            list.push_back(new fake::ff_fake());
            return list;
        }

        const std::string elem_id;
    };

}
#endif
