/*
 * potential2D.cpp
 *
 *  Created on: 8 nov. 2015
 *      Author: aschils
 */

#include "../model/LaplaceSolver.hpp"
#include "../model/Rect2DBoundaryValues.hpp"
#include "../model/ZeroRightHandSide.hpp"
#include "../model/Rect2DDetector.hpp"
#include "../model/SerratedRect2DDetector.hpp"

void test_2D_potential() {

	ZeroRightHandSide<2> rhs;
	double strip_potential = 1;
	unsigned strip_length = 100;
	unsigned pitch = 100;

	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 5; nbr_of_strips++) {
		Rect2DDetector rdd(nbr_of_strips, strip_length, pitch, strip_potential);
		rdd.compute_potential();
	}
}

void test_serrated_2D_potential() {

	ZeroRightHandSide<2> rhs;
	double strip_potential = 1;
	unsigned strip_length = 100;
	unsigned strip_width = 30;
	unsigned pitch = 100;

	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 0; nbr_of_strips++) {
		SerratedRect2DDetector srdd(nbr_of_strips, strip_length, strip_width, pitch, strip_potential);
		srdd.compute_potential();
	}
}

