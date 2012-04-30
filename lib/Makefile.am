#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

#tensor
nobase_include_HEADERS = 	radmat/tensor/tensorbase.h 	\
				radmat/tensor/tensor_fwd.h	\
				radmat/tensor/narray.h		\
				radmat/tensor/tensor.h		\
				radmat/tensor/tensor_algebra.h	\
				radmat/tensor/test_NArray.h	\
				radmat/tensor/test_Tensor.h	\
				radmat/tensor/test_utils.h	


#utils
nobase_include_HEADERS += 	radmat/utils/stringify.h			\
				radmat/utils/pow2assert.h			\
				radmat/utils/breit_frame.h			\
				radmat/utils/wigner.h 				\
				radmat/utils/euler_angles.h			\
				radmat/utils/perThreadStorage.h			\
				radmat/utils/aux.h

#polarisation
nobase_include_HEADERS +=	radmat/polarisation/polarisation_factory_keys.h	\
				radmat/polarisation/polarisation_factory_inventory.h	\
				radmat/polarisation/polarisation_factory.h		\
				radmat/polarisation/polarisationbase.h



#
# Library
#

lib_LIBRARIES = libradmat.a

#tensor
libradmat_a_SOURCES = 		radmat/tensor/test_utils.cc			\
				radmat/tensor/tensor_algebra.cc

#utils
libradmat_a_SOURCES +=		radmat/utils/pow2assert.cc 			\
				radmat/utils/breit_frame.cc			\
				radmat/utils/wigner.cc				\
				radmat/utils/euler_angles.cc

#polarisation
libradmat_a_SOURCES +=		radmat/polarisation/polarisation_factory_keys.cc	\
				radmat/polarisation/polarisation_factory_inventory.cc	\
				radmat/polarisation/polarisation_factory.cc

