#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/lib \
               @CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ 

AM_LDFLAGS  = -L@top_builddir@/lib \
               @LDFLAGS@ @ADAT_LDFLAGS@

LDADD = @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ -lradmat

check_PROGRAMS = test_CLEBSCH

test_CLEBSCH_SOURCES = test_clebsch.cc

#
# The following rule will always try to rebuild libadat.a when you 
# try to compile anything else
#
####DEPENDENCIES = build_lib ${lib_LIBRARIES}
DEPENDENCIES = ${lib_LIBRARIES}
${noinst_PROGRAMS}: ${DEPENDENCIES}
${bin_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)