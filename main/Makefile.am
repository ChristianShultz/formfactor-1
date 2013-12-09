#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @SEMBLE_CXXFLAGS@ @JACKFITTER_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@

AM_LDFLAGS  = -L@top_builddir@/lib/ \
               @LDFLAGS@ @SEMBLE_LDFLAGS@ @JACKFITTER_LDFLAGS@ @ADAT_LDFLAGS@ 

LDADD = -lradmat @SEMBLE_LIBS@ @JACKFITTER_LIBS@ @ADAT_LIBS@ @ITPP_LIBS@ @LIBS@ 

bin_PROGRAMS =   fake_radmat mass_overlap_db_create mass_overlap_db_combine single_plot phases mass_overlap_db_util grab_ff grab_pol grab_vector_mean_cmp invert_subduced phases2 check_rotations corr_norm radmat_util radmat test_tens test_redstar_ops


fake_radmat_SOURCES = fake_radmat.cc

mass_overlap_db_create_SOURCES = mass_overlap_db_create.cc

mass_overlap_db_combine_SOURCES = mass_overlap_db_combine.cc

mass_overlap_db_util_SOURCES = mass_overlap_db_util.cc

single_plot_SOURCES = single_plot.cc

phases_SOURCES = phases.cc

phases2_SOURCES = phases2.cc

grab_ff_SOURCES = grab_ff.cc

grab_pol_SOURCES = grab_pol.cc

grab_vector_mean_cmp_SOURCES = grab_vector_mean_cmp.cc

invert_subduced_SOURCES = invert_subduced.cc

check_rotations_SOURCES = check_rotations.cc

corr_norm_SOURCES = corr_norm.cc

radmat_util_SOURCES = radmat_util.cc

radmat_SOURCES = radmat.cc

test_tens_SOURCES = test_tens.cc

test_redstar_ops_SOURCES = test_redstar_ops.cc


# test_jfit_SOURCES = test_jfit.cc

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
