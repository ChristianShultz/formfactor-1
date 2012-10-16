#ifndef MAKE_FAKE_DATABASE_H
#define MAKE_FAKE_DATABASE_H



#include <cstdio>
#include <cctype>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "adat/singleton.h"
#include "adat/funcmap.h"

#include "dbtype.h"
#include "formfac/hadron_1pt_corr.h"
#include "formfac/hadron_2pt_corr.h"
#include "formfac/hadron_3pt_corr.h"

#include "hadron/hadron_npart_npt_corr.h"

#include "io/key_val_db.h"

#include "AllConfStoreDB.h"
#include "ConfDataDBMerger.h"




using namespace std;
using namespace FILEDB;
using namespace ADATIO;
using namespace Util;
using namespace FF;










 
#endif /* MAKE_FAKE_DATABASE_H */
