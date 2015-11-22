/*
 * SerratedRect2DDetector.hpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#ifndef __SERRATED_RECT_2D_DETECTOR_HPP__
#define __SERRATED_RECT_2D_DETECTOR_HPP__

#include "SerratedLaplaceSolver.hpp"
#include "ZeroRightHandSide.hpp"
#include "SerratedRect2DBoundaryValues.hpp"

class SerratedRect2DDetector {

public:
	SerratedRect2DDetector(unsigned nbr_of_strips, unsigned strip_length,
			unsigned strip_width, unsigned pitch, double strip_potential);
	~SerratedRect2DDetector();
	void compute_potential();

private:
	unsigned nbr_of_strips, strip_length, strip_width, pitch, total_length = 1;
	double strip_potential = 1.0;
	const double RECT_WIDTH = 300.0; //i.e. in domain language (microm,..)
	const double RECT_LENGTH_FE = 10000.0;
	double rect_width_fe = 1.0;

	ZeroRightHandSide<2> *zero_right_hand_side;
	SerratedRect2DBoundaryValues<2> *boundary_val;
	SerratedLaplaceSolver<2> *rect_potential_solver;

	double compute_total_length();
	double compute_rect_width_fe();
};

#endif
