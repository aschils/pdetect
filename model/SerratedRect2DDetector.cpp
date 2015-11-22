/*
 * SerratedRect2DDetector.cpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#include "SerratedRect2DDetector.hpp"

/**
 * Compute the length of the rectangle (detector) in the domain language unit
 * (microm) depending on the number of strips, the strip length and the pitch.
 */
double SerratedRect2DDetector::compute_total_length(){
	if(nbr_of_strips == 0 || strip_length == 0)
		return RECT_WIDTH;
	return nbr_of_strips * (strip_length + pitch);
}

/**
 * Rectangle width should always be 300microm. Thus finite element rectangle
 * width should be adapted depending on the length of the rectangle in microm.
 */
double SerratedRect2DDetector::compute_rect_width_fe() {
	double one_micron_fe = RECT_LENGTH_FE/total_length;
	return RECT_WIDTH*one_micron_fe;
}

SerratedRect2DDetector::SerratedRect2DDetector(unsigned nbr_of_strips,
		unsigned strip_length, unsigned strip_width, unsigned pitch,
		double strip_potential) {

	this->nbr_of_strips = nbr_of_strips;
	this->strip_length = strip_length;
	this->pitch = pitch;
	this->strip_potential = strip_potential;
	total_length = compute_total_length();
	rect_width_fe = compute_rect_width_fe();

	zero_right_hand_side = new ZeroRightHandSide<2>();
	boundary_val = new SerratedRect2DBoundaryValues<2>(nbr_of_strips, strip_length,
			strip_width,
			pitch, RECT_LENGTH_FE, rect_width_fe, strip_potential);


	/*
	 * double rect_length_fe,
			double rect_width_fe, unsigned nbr_of_strips, unsigned strip_length,
			unsigned strip_width, unsigned pitch,
			const Function<dim> *right_hand_side,
			Function<dim> *boundary_values_fun, std::string result_file_path
	 */

	rect_potential_solver = new SerratedLaplaceSolver<2>(RECT_LENGTH_FE, rect_width_fe,
			nbr_of_strips, strip_length, strip_width, pitch,
			zero_right_hand_side, boundary_val,
			std::to_string(nbr_of_strips) + ".vtk");
}

SerratedRect2DDetector::~SerratedRect2DDetector() {
	delete zero_right_hand_side;
	delete boundary_val;
	delete rect_potential_solver;
}

void SerratedRect2DDetector::compute_potential() {
	rect_potential_solver->run();
}

