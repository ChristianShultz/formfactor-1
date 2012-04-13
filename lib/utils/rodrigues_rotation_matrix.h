#ifndef RODRIGUES_ROTATION_MATRIX_H_H_GUARD
#define RODRIGUES_ROTATION_MATRIX_H_H_GUARD

#include "itpp/itbase.h"
#include "tensor/tensorbase.h"


/**
   @file rodrigues_rotation_matrix.h
   
   @brief determine the rotation matrix between two 3-vectors
 */

//! rotated = R*orig
itpp::Mat<double> rodRotMat(const itpp::Vec<double> &orig,
			    const itpp::Vec<double> &rotated);

// assumes a 4 vector and gets rotation for compontents orig[1-3]
//! for the spatial compontent rotated = R*orig
itpp::Mat<double> rodRotMat(const tensor::Tensor<double, 1> &orig, 
			    const tensor::Tensor<double, 1> &rotated);



#endif
