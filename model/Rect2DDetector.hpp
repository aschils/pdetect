/*
 * Rect2DDetector.hpp
 *
 *  Created on: 10 nov. 2015
 *      Author: aschils
 */

#ifndef __RECT_2D_DETECTOR_HPP__
#define __RECT_2D_DETECTOR_HPP__

#include "LaplaceSolver.hpp"
#include "ZeroRightHandSide.hpp"
#include "Rect2DBoundaryValues.hpp"

class Rect2DDetector {

public:
	Rect2DDetector(unsigned nbr_of_strips, unsigned strip_length,
			unsigned pitch, double strip_potential);
	~Rect2DDetector();
	void compute_potential();

private:
	unsigned nbr_of_strips, strip_length, pitch, total_length = 1;
	double strip_potential = 1.0;
	const double RECT_WIDTH = 300.0; //i.e. in domain language (microm,..)
	const double RECT_LENGTH_FE = 10000.0;
	double rect_width_fe = 1.0;

	ZeroRightHandSide<2> *zero_right_hand_side;
	Rect2DBoundaryValues<2> *boundary_val;
	LaplaceSolver<2> *rect_potential_solver;

	double compute_total_length();
	double compute_rect_width_fe();
};

#endif
