#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/lib \
               @CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ 

AM_LDFLAGS  = -L@top_builddir@/lib \
               @LDFLAGS@ @ADAT_LDFLAGS@

LDADD = @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ -lradmat

bin_PROGRAMS = hello

hello_SOURCES = hello.cc


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