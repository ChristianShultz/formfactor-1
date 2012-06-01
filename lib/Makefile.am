#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

#utils
nobase_include_HEADERS = radmat/utils/aux.h radmat/utils/breit_frame.h radmat/utils/euler_angles.h radmat/utils/lorentz.h radmat/utils/perThreadStorage.h radmat/utils/polarisation_tensors.h radmat/utils/pow2assert.h radmat/utils/stringify.h radmat/utils/tensor.h radmat/utils/wigner.h

#lorentz_ff
#nobase_include_HEADERS += radmat/lorentz_ff/ff_base.h radmat/lorentz_ff/ff_llsq_row.h radmat/lorentz_ff/ff_loader.h radmat/lorentz_ff/formfac/fake.h radmat/lorentz_ff/formfac/ff_0p0p.h

#ff
nobase_include_HEADERS += radmat/ff/common_cfg.h radmat/ff/lorentzff_PiPi.h

#
# Library
#

lib_LIBRARIES = libradmat.a

#utils
libradmat_a_SOURCES = radmat/utils/breit_frame.cc radmat/utils/euler_angles.cc radmat/utils/pow2assert.cc radmat/utils/wigner.cc radmat/utils/tensor.cc

#ff
libradmat_a_SOURCES += radmat/ff/common_cfg.cc radmat/ff/lorentzff_PiPi.cc
