#
# Path and flags for all files used in building the library
#

AM_CXXFLAGS = 	-I@top_srcdir@/lib		\
								@CXXFLAGS@ @SEMBLE_CXXFLAGS@ @ADAT_CXXFLAGS@ @ITPP_CXXFLAGS@ @JACKFITTER_CXXFLAGS@

#utils
nobase_include_HEADERS =   \
													 radmat/construct_data/construct_correlators_bad_data_repository.h		\
													 radmat/construct_data/construct_correlators.h		\
													 radmat/construct_data/construct_correlators_utils.h		\
													 radmat/construct_data/construct_correlators_xml.h		\
													 radmat/construct_data/construct_data.h		\
													 radmat/construct_data/construct_q2_packs.h		\
													 radmat/construct_data/lattice_multi_data_object.h		\
													 radmat/construct_data/lattice_multi_data_tag.h		\
													 radmat/construct_data/lattice_multi_data_tag_redstar_interface.h		\
													 radmat/construct_data/load_fake_data.h		\
													 radmat/construct_data/radmat_database_interface.h		\
													 radmat/construct_data/radmat_overlap_key_val_db.h		\
													 radmat/construct_data/three_point.h		\
													 radmat/driver/driver.h		\
													 radmat/driver/radmat_driver_aux.h		\
													 radmat/driver/radmat_driver.h		\
													 radmat/driver/radmat_driver_props.h		\
													 radmat/driver/radmat_single_q2_driver.h		\
													 radmat/fake_data/covarrying_vectors.h		\
													 radmat/fake_data/fake_3pt_function_aux.h		\
													 radmat/fake_data/fake_3pt_function.h		\
													 radmat/fake_data/fake_data.h		\
													 radmat/fake_data/fake_data_ini.h		\
													 radmat/fake_data/fake_overlaps.h		\
													 radmat/fake_data/fake_spectrum.h		\
													 radmat/fake_data/gen_fake_dataset.h		\
													 radmat/fake_data/minimal_fake_data.h		\
													 radmat/fake_data/normalize_fake_3pt_data.h		\
													 radmat/ff/ff_gen_llsq_row.h		\
													 radmat/ff/ff.h		\
													 radmat/ff/formfactor_factory.h		\
													 radmat/ff/lorentzff_canonical_rotations.h		\
													 radmat/ff/lorentzff_canonical_rotation_groups.h		\
													 radmat/ff/lorentzff_canonical_rotations_utils.h		\
													 radmat/ff/lorentzff_PiPi.h		\
													 radmat/ff/lorentzff_PiPiStar.h		\
													 radmat/ff/lorentzff_PiRho.h		\
													 radmat/ff/lorentzff_polarization_embedding.h		\
													 radmat/ff/lorentzff_polarization_tensors.h		\
													 radmat/ff/lorentzff_RhoPi.h		\
													 radmat/ff/lorentzff_RhoRho.h		\
													 radmat/ff/lorentzff_Zpolarization_tensors.h		\
													 radmat/ff/rest_polarization_convention.h		\
													 radmat/fitting/axis_plotter.h		\
													 radmat/fitting/fit_tins.h		\
													 radmat/llsq/llsq_chisq.h		\
													 radmat/llsq/llsq_driver.h		\
													 radmat/llsq/llsq_gen_system.h		\
													 radmat/llsq/llsq.h		\
													 radmat/llsq/llsq_multi_data.h		\
													 radmat/llsq/llsq_multi_data_serialize.h		\
													 radmat/llsq/llsq_multi_driver.h		\
													 radmat/llsq/llsq_q2_pack.h		\
													 radmat/llsq/llsq_solution.h		\
													 radmat/llsq/llsq_solver.h		\
													 radmat/llsq/llsq_solvers.h		\
													 radmat/llsq/llsq_zero_filter.h		\
													 radmat/redstar_interface/abstract_named_object.h		\
													 radmat/redstar_interface/redstar_abstract_block_base.h		\
													 radmat/redstar_interface/redstar_abstract_merge_npoint_factory.h		\
													 radmat/redstar_interface/redstar_abstract_merge_npoint.h		\
													 radmat/redstar_interface/redstar_abstract_xml_factory.h		\
													 radmat/redstar_interface/redstar_canonical_lattice_rotations.h		\
													 radmat/redstar_interface/redstar_canonical_rotations.h		\
													 radmat/redstar_interface/redstar_cartesian_interface.h		\
													 radmat/redstar_interface/redstar_improved_vector_current.h		\
													 radmat/redstar_interface/redstar_interface.h		\
													 radmat/redstar_interface/redstar_invert_subduction.h		\
													 radmat/redstar_interface/redstar_merge_function_xml.h		\
													 radmat/redstar_interface/redstar_merge_vector_current_npoint.h		\
													 radmat/redstar_interface/redstar_npoint_function_xml.h		\
													 radmat/redstar_interface/redstar_photon_polarization_tensor.h		\
													 radmat/redstar_interface/redstar_single_particle_meson_block.h		\
													 radmat/redstar_interface/redstar_unimproved_vector_current.h		\
													 radmat/register_all/register_all.h		\
													 radmat/utils/aux.h		\
													 radmat/utils/debug_handler.h		\
													 radmat/utils/handle.h		\
													 radmat/utils/harmonic.h		\
													 radmat/utils/levi_civita.h		\
													 radmat/utils/map_obj.h		\
													 radmat/utils/mink_qsq.h		\
													 radmat/utils/obj_expr_t.h		\
													 radmat/utils/perThreadStorage.h		\
													 radmat/utils/polarization_tensors.h		\
													 radmat/utils/pow2assert.h		\
													 radmat/utils/splash.h		\
													 radmat/utils/stringify.h		\
													 radmat/utils/tensor.h		\
													 radmat/utils/type_computations.h		\
													 radmat/utils/utils.h		

# the lib
lib_LIBRARIES = libradmat.a


libradmat_a_SOURCES =	\
											radmat/construct_data/construct_correlators_bad_data_repository.cc		\
											radmat/construct_data/construct_correlators.cc		\
											radmat/construct_data/construct_correlators_utils.cc		\
											radmat/construct_data/construct_correlators_xml.cc		\
											radmat/construct_data/lattice_multi_data_tag.cc		\
											radmat/construct_data/lattice_multi_data_tag_redstar_interface.cc		\
											radmat/construct_data/radmat_database_interface.cc		\
											radmat/construct_data/radmat_overlap_key_val_db.cc		\
											radmat/driver/radmat_driver_aux.cc		\
											radmat/driver/radmat_driver.cc		\
											radmat/driver/radmat_driver_props.cc		\
											radmat/driver/radmat_single_q2_driver.cc		\
											radmat/fake_data/covarrying_vectors.cc		\
											radmat/fake_data/fake_3pt_function_aux.cc		\
											radmat/fake_data/fake_data_ini.cc		\
											radmat/fake_data/fake_overlaps.cc		\
											radmat/fake_data/fake_spectrum.cc		\
											radmat/ff/ff_gen_llsq_row.cc		\
											radmat/ff/formfactor_factory.cc		\
											radmat/ff/lorentzff_canonical_rotations.cc		\
											radmat/ff/lorentzff_canonical_rotations_utils.cc		\
											radmat/ff/lorentzff_PiPi.cc		\
											radmat/ff/lorentzff_PiPiStar.cc		\
											radmat/fitting/axis_plotter.cc		\
											radmat/fitting/fit_tins.cc		\
											radmat/llsq/llsq_chisq.cc		\
											radmat/llsq/llsq_gen_system.cc		\
											radmat/llsq/llsq_multi_driver.cc		\
											radmat/llsq/llsq_q2_pack.cc		\
											radmat/llsq/llsq_solvers.cc		\
											radmat/redstar_interface/redstar_abstract_merge_npoint_factory.cc		\
											radmat/redstar_interface/redstar_abstract_xml_factory.cc		\
											radmat/redstar_interface/redstar_canonical_lattice_rotations.cc		\
											radmat/redstar_interface/redstar_canonical_rotations.cc		\
											radmat/redstar_interface/redstar_cartesian_interface.cc		\
											radmat/redstar_interface/redstar_improved_vector_current.cc		\
											radmat/redstar_interface/redstar_invert_subduction.cc		\
											radmat/redstar_interface/redstar_merge_function_xml.cc		\
											radmat/redstar_interface/redstar_merge_vector_current_npoint.cc		\
											radmat/redstar_interface/redstar_npoint_function_xml.cc		\
											radmat/redstar_interface/redstar_single_particle_meson_block.cc		\
											radmat/redstar_interface/redstar_unimproved_vector_current.cc		\
											radmat/register_all/register_all.cc		\
											radmat/utils/aux.cc		\
											radmat/utils/harmonic.cc		\
											radmat/utils/mink_qsq.cc		\
											radmat/utils/polarization_tensors.cc		\
											radmat/utils/pow2assert.cc		\
											radmat/utils/splash.cc		\
											radmat/utils/tensor.cc		
