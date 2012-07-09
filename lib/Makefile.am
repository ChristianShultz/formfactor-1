#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @SEMBLE_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

#utils
nobase_include_HEADERS = radmat/utils/aux.h radmat/utils/breit_frame.h radmat/utils/euler_angles.h radmat/utils/lorentz.h radmat/utils/perThreadStorage.h radmat/utils/polarisation_tensors.h radmat/utils/pow2assert.h radmat/utils/stringify.h radmat/utils/tensor.h radmat/utils/wigner.h radmat/utils/splash.h radmat/utils/harmonic.h

#ff
nobase_include_HEADERS += radmat/ff/ff_gen_llsq_row.h radmat/ff/formfactor_factory.h radmat/ff/llsq_gen_system.h radmat/ff/llsq_solver.h radmat/ff/llsq_solvers.h radmat/ff/lorentzff_PiPi.h

#fake_data
nobase_include_HEADERS += radmat/fake_data/covarrying_vectors.h radmat/fake_data/fake_data.h radmat/fake_data/fake_data_ini.h radmat/fake_data/fake_overlaps.h radmat/fake_data/fake_spectrum.h radmat/fake_data/minimal_fake_data.h



#
# Library
#

lib_LIBRARIES = libradmat.a

#utils
libradmat_a_SOURCES = radmat/utils/breit_frame.cc radmat/utils/euler_angles.cc radmat/utils/pow2assert.cc radmat/utils/wigner.cc radmat/utils/tensor.cc radmat/utils/splash.cc radmat/utils/harmonic.cc

#ff
libradmat_a_SOURCES += radmat/ff/ff_gen_llsq_row.cc radmat/ff/formfactor_factory.cc radmat/ff/llsq_solvers.cc radmat/ff/lorentzff_PiPi.cc

#fake_data 
libradmat_a_SOURCES +=radmat/fake_data/covarrying_vectors.cc radmat/fake_data/fake_data_ini.cc radmat/fake_data/fake_overlaps.cc radmat/fake_data/fake_spectrum.cc
