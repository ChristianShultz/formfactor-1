#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

#tensor
nobase_include_HEADDERS = 	tensor/tensorbase.h 	\
				tensor/tensor_fwd.h	\
				tensor/narray.h		\
				tensor/tensor.h		\
				tensor/tensor_algebra.h	\
				tensor/test_NArray.h	\
				tensor/test_Tensor.h	\
				tensor/test_utils.h	\
				tensor/tensor_aux.h


#utils
nobase_include_HEADDERS += 	utils/stringify.h			\
				utils/pow2assert.h			\
				utils/rodrigues_rotation_matrix.h	\
				utils/lorentz.h

#polarisation
nobase_include_HEADDERS +=	polarisation/polarisation_factory_keys.h	\
				polarisation/polarisation_factory_inventory.h	\
				polarisation/polarisation_factory.h		\
				polarisation/polarisationbase.h

#kinematic_factors
nobase_include_HEADDERS += 	kinematic_factors/kinfac_fwd.h		\
				kinematic_factors/kinfac_key.h		\
				kinematic_factors/kinematic_factors.h	\
				kinematic_factors/kinfacbase.h

#
# Library
#

lib_LIBRARIES = libradmat.a

#tensor
libradmat_a_SOURCES = 		tensor/test_utils.cc			\
				tensor/tensor_algebra.cc

#utils
libradmat_a_SOURCES +=		utils/pow2assert.cc 			\
				utils/rodrigues_rotation_matrix.cc	\
				utils/lorentz.cc

#polarisation
libradmat_a_SOURCES +=		polarisation/polarisation_factory_keys.cc	\
				polarisation/polarisation_factory_inventory.cc	\
				polarisation/polarisation_factory.cc

#kinematic_factors
libradmat_a_SOURCES += 		kinematic_factors/kinfac_key.cc