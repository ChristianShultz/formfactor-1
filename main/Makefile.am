#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @SEMBLE_CXXFLAGS@ @JACKFITTER_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

AM_LDFLAGS  = -L@top_builddir@/lib/ \
               @LDFLAGS@ @SEMBLE_LDFLAGS@ @JACKFITTER_LDFLAGS@ @ADAT_LDFLAGS@

LDADD = -lradmat @SEMBLE_LIBS@ @JACKFITTER_LIBS@ @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ 

noinst_PROGRAMS = hello plot_fake_3pt

hello_SOURCES = hello.cc

plot_fake_3pt_SOURCES = plot_fake_3pt.cc

# no idea how/why this works, cut paste job
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