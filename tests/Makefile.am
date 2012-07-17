#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/lib \
               @CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ @SEMBLE_CXXFLAGS@

AM_LDFLAGS  = -L@top_builddir@/lib/ \
               @LDFLAGS@ @ADAT_LDFLAGS@ @SEMBLE_LDFLAGS@

LDADD = -lradmat @SEMBLE_LIBS@ @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ 


check_PROGRAMS = test test_fake_ini test_cov_vecs fake_3pt_play

test_SOURCES = src/test_polarisation_tensor.cc src/test_ff_pipi.cc src/test_LLSQ_solver.cc src/test_factories.cc src/test_minimal_fake_data.cc src/test_covarrying_vector.cc src/test_make_fake_overlaps.cc src/test_make_fake_spectrum.cc src/test_fake_3pt_aux.cc src/test_fake_3pt.cc src/test_make_fake_data_ini.cc src/test.cc 

test_fake_ini_SOURCES = src/test_fake_ini.cc

test_cov_vecs_SOURCES = src/test_cov_vecs.cc

fake_3pt_play_SOURCES = src/fake_3pt_play.cc

# no idea how/why (if?) this works, cut paste job
#
# The folow rule will always try to rebuild libadat.a when you 
# try to compile a main prog
#
DEPENDENCIES = build_lib
${bin_PROGRAMS}: ${DEPENDENCIES}   # uncomment this when there are binprograms
#${check_PROGRAMS}: ${DEPENDENCIES}
#${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
