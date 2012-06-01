#ifndef FF_LLSQ_ROW_H_H_GUARD
#define FF_LLSQ_ROW_H_H_GUARD

#include <list>
#include <iostream>
#include <string>
#include "ff_base.h"
#include "ff_loader.h"


/**
   @file ff_llsq_row.h
   @brief defines a class to generate a row of the linear system at the cfg level
 */

namespace radmat
{

    struct FFLLSqRowBase_t
    {
        typedef ff_loader::ff_list ff_list;
        typedef std::list<FFBlock_rt> ff_row;
        typedef FFBlockBase_t::pack_handle pack_handle;

        FFLLSqRowBase_t(void) ; // hidden

        FFLLSqRowBase_t(const pack_handle &_pack_f, const pack_handle &_pack_i, const std::string type)
            : pack_f(_pack_f) , pack_i(_pack_i)
        {
            load_ffs(type);
        }

        ~FFLLSqRowBase_t(void)
        {
            ff_list::iterator first = m_list.begin(), last = m_list.end();

            while(first != last)
                delete *first++;

            m_list.clear();
        }


        ff_row operator()(void)
        {
            ff_list::const_iterator first = m_list.begin();
            ff_list::const_iterator last = m_list.end();
            ff_row ret(m_list.size());
            ff_row::iterator result = ret.begin();

            while(first != last)
                {
                    *result = (**first)(pack_f, pack_i);
                    result ++;
                    first ++;
                }

            return ret;
        }

    private:

        void load_ffs(const std::string type)
        {
            std::cout << __func__ << std::endl;
            std::cout << "warning hardwired fake for dev" << std::endl;
            m_list = ff_loader("fake").genList();
        }

        pack_handle pack_f, pack_i;
        ff_list m_list;
    };

}

#endif
