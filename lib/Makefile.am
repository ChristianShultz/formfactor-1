#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
		@CXXFLAGS@ @SEMBLE_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ @JACKFITTER_CXXFLAGS@

#utils
nobase_include_HEADERS = radmat/utils/aux.h radmat/utils/breit_frame.h radmat/utils/euler_angles.h radmat/utils/harmonic.h radmat/utils/lorentz.h radmat/utils/obj_expr_t.h radmat/utils/perThreadStorage.h radmat/utils/polarisation_tensors.h radmat/utils/pow2assert.h radmat/utils/splash.h radmat/utils/stringify.h radmat/utils/tensor.h radmat/utils/type_computations.h radmat/utils/wigner.h radmat/utils/levi_civita.h radmat/utils/mink_qsq.h

#ff
nobase_include_HEADERS += radmat/ff/ff_gen_llsq_row.h radmat/ff/formfactor_factory.h radmat/ff/lorentzff_PiPi.h radmat/ff/lorentzff_PiPiStar.h radmat/ff/lorentzff_polarization_embedding.h radmat/ff/lorentzff_PiRho.h radmat/ff/lorentzff_RhoPi.h

#fake_data
nobase_include_HEADERS += radmat/fake_data/covarrying_vectors.h radmat/fake_data/fake_3pt_function.h radmat/fake_data/fake_3pt_function_aux.h radmat/fake_data/fake_data_ini.h radmat/fake_data/fake_overlaps.h radmat/fake_data/fake_spectrum.h radmat/fake_data/gen_fake_dataset.h radmat/fake_data/minimal_fake_data.h radmat/fake_data/normalize_fake_3pt_data.h

#llsq
nobase_include_HEADERS += radmat/llsq/llsq_multi_data_serialize.h radmat/llsq/llsq_driver.h radmat/llsq/llsq_gen_system.h radmat/llsq/llsq_multi_data.h radmat/llsq/llsq_multi_driver.h radmat/llsq/llsq_q2_pack.h radmat/llsq/llsq_solver.h radmat/llsq/llsq_solvers.h radmat/llsq/llsq_zero_filter.h radmat/llsq/llsq_chisq.h


#fitting
nobase_include_HEADERS += radmat/fitting/axis_plotter.h radmat/fitting/fit_tins.h

#driver
nobase_include_HEADERS += radmat/driver/radmat_driver.h radmat/driver/radmat_driver_props.h radmat/driver/radmat_single_q2_driver.h radmat/driver/radmat_driver_aux.h


#jfit

# nobase_include_HEADERS += radmat/jfit/abstract_function_default_errors.h radmat/jfit/abstract_function_default_values.h radmat/jfit/abstract_function_exit.h radmat/jfit/abstract_function.h radmat/jfit/abstract_function_named_parameter_container.h radmat/jfit/abstract_function_parameter_container.h radmat/jfit/abstract_function_parameter_fixing.h radmat/jfit/abstract_function_parameter_limits.h radmat/jfit/abstract_function_parameter_names.h

# redstar_interface 
nobase_include_HEADERS += radmat/redstar_interface/abstract_named_object.h radmat/redstar_interface/redstar_abstract_block_base.h radmat/redstar_interface/redstar_abstract_merge_npoint_factory.h radmat/redstar_interface/redstar_abstract_merge_npoint.h radmat/redstar_interface/redstar_abstract_xml_factory.h radmat/redstar_interface/redstar_interface.h radmat/redstar_interface/redstar_merge_function_xml.h radmat/redstar_interface/redstar_merge_vector_current_npoint.h radmat/redstar_interface/redstar_npoint_function_xml.h radmat/redstar_interface/redstar_single_particle_meson_block.h radmat/redstar_interface/redstar_unimproved_vector_current.h

# construct_data
nobase_include_HEADERS += radmat/construct_data/construct_correlators_bad_data_repository.h radmat/construct_data/construct_correlators.h radmat/construct_data/construct_correlators_utils.h radmat/construct_data/construct_correlators_xml.h radmat/construct_data/construct_q2_packs.h radmat/construct_data/invert_subduction.h radmat/construct_data/lattice_multi_data_object.h radmat/construct_data/lattice_multi_data_tag.h radmat/construct_data/lattice_multi_data_tag_redstar_interface.h radmat/construct_data/load_fake_data.h radmat/construct_data/radmat_database_interface.h radmat/construct_data/radmat_overlap_key_val_db.h radmat/construct_data/three_point.h

#
# Library
#

lib_LIBRARIES = libradmat.a

#utils
libradmat_a_SOURCES = radmat/utils/aux.cc radmat/utils/breit_frame.cc radmat/utils/euler_angles.cc radmat/utils/harmonic.cc radmat/utils/pow2assert.cc radmat/utils/splash.cc radmat/utils/tensor.cc radmat/utils/wigner.cc radmat/utils/polarisation_tensors.cc radmat/utils/mink_qsq.cc

#ff
libradmat_a_SOURCES += radmat/ff/ff_gen_llsq_row.cc radmat/ff/formfactor_factory.cc radmat/ff/lorentzff_PiPi.cc radmat/ff/lorentzff_PiPiStar.cc

#fake_data 
libradmat_a_SOURCES += radmat/fake_data/covarrying_vectors.cc radmat/fake_data/fake_3pt_function_aux.cc radmat/fake_data/fake_data_ini.cc radmat/fake_data/fake_overlaps.cc radmat/fake_data/fake_spectrum.cc

#llsq
libradmat_a_SOURCES += radmat/llsq/llsq_gen_system.cc radmat/llsq/llsq_multi_driver.cc radmat/llsq/llsq_q2_pack.cc radmat/llsq/llsq_solvers.cc radmat/llsq/llsq_chisq.cc


#fitting
libradmat_a_SOURCES += radmat/fitting/axis_plotter.cc radmat/fitting/fit_tins.cc

#driver
libradmat_a_SOURCES += radmat/driver/radmat_driver.cc radmat/driver/radmat_driver_props.cc radmat/driver/radmat_single_q2_driver.cc radmat/driver/radmat_driver_aux.cc

#redstar_interface
libradmat_a_SOURCES +=radmat/redstar_interface/redstar_abstract_merge_npoint_factory.cc radmat/redstar_interface/redstar_abstract_xml_factory.cc radmat/redstar_interface/redstar_merge_function_xml.cc radmat/redstar_interface/redstar_merge_vector_current_npoint.cc radmat/redstar_interface/redstar_npoint_function_xml.cc radmat/redstar_interface/redstar_single_particle_meson_block.cc radmat/redstar_interface/redstar_unimproved_vector_current.cc

#construct_data
libradmat_a_SOURCES += radmat/construct_data/construct_correlators_bad_data_repository.cc radmat/construct_data/construct_correlators.cc radmat/construct_data/construct_correlators_utils.cc radmat/construct_data/construct_correlators_xml.cc radmat/construct_data/invert_subduction.cc radmat/construct_data/lattice_multi_data_tag.cc radmat/construct_data/lattice_multi_data_tag_redstar_interface.cc radmat/construct_data/radmat_database_interface.cc radmat/construct_data/radmat_overlap_key_val_db.cc

