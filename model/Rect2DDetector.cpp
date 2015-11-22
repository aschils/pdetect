/*
 * Rect2DDetector.cpp
 *
 *  Created on: 10 nov. 2015
 *      Author: aschils
 */

#include "Rect2DDetector.hpp"

/**
 * Compute the length of the rectangle (detector) in the domain language unit
 * (microm) depending on the number of strips, the strip length and the pitch.
 */
double Rect2DDetector::compute_total_length() {
	if (nbr_of_strips == 0 || strip_length == 0)
		return RECT_WIDTH;
	return nbr_of_strips * (strip_length + pitch);
}

/**
 * Rectangle width should always be 300microm. Thus finite element rectangle
 * width should be adapted depending on the length of the rectangle in microm.
 */
double Rect2DDetector::compute_rect_width_fe() {
	double one_micron_fe = RECT_LENGTH_FE / total_length;
	return RECT_WIDTH * one_micron_fe;
}

Rect2DDetector::Rect2DDetector(unsigned nbr_of_strips, unsigned strip_length,
		unsigned pitch, double strip_potential, unsigned refine_level,
		unsigned max_iter, double stop_accuracy, std::string output_file) {

	this->nbr_of_strips = nbr_of_strips;
	this->strip_length = strip_length;
	this->pitch = pitch;
	this->strip_potential = strip_potential;
	this->triangulation = new Triangulation<2>();
	total_length = compute_total_length();
	rect_width_fe = compute_rect_width_fe();

	zero_right_hand_side = new ZeroRightHandSide<2>();
	boundary_val = new Rect2DBoundaryValues<2>(nbr_of_strips, strip_length,
			pitch, RECT_LENGTH_FE, rect_width_fe, strip_potential);

	MyGridGenerator<2>::hyper_rectangle(*triangulation, RECT_LENGTH_FE,
			rect_width_fe);

	rect_potential_solver = new LaplaceSolver<2>(triangulation, refine_level,
			max_iter, stop_accuracy, zero_right_hand_side, boundary_val,
			output_file, true);

}

Rect2DDetector::~Rect2DDetector() {
	delete zero_right_hand_side;
	delete boundary_val;
	delete rect_potential_solver;
	delete triangulation;
}

void Rect2DDetector::compute_potential() {
	rect_potential_solver->run();
}

