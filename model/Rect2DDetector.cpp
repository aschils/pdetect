/*
 * Rect2DDetector.cpp
 *
 *  Created on: 10 nov. 2015
 *      Author: aschils
 */

#include "Rect2DDetector.hpp"

/**
 * Rectangle width should always be 300microm. So finite element rectangle
 * width should be adapted depending on the length of the rectangle in microm.
 */
void compute_rect_width_fe(){

}

Rect2DDetector::Rect2DDetector(unsigned nbr_of_strips, unsigned strip_length,
		unsigned pitch, double strip_potential) {
	//TODO compute rect_length_fe depending on problem size
	double rect_length_fe = 10000;
	double rect_width_fe = 3000;
	zero_right_hand_side = new ZeroRightHandSide<2>();
	boundary_val = new Rect2DBoundaryValues<2>(nbr_of_strips, strip_length,
			pitch, rect_length_fe, rect_width_fe, strip_potential);
	rect_potential_solver = new LaplaceSolver<2>(rect_length_fe, rect_width_fe,
			zero_right_hand_side, boundary_val,
			std::to_string(nbr_of_strips) + ".vtk");

	this->nbr_of_strips = nbr_of_strips;
	this->strip_length = strip_length;
	this->pitch = pitch;
	this->strip_potential = strip_potential;
}

Rect2DDetector::~Rect2DDetector() {
	delete zero_right_hand_side;
	delete boundary_val;
	delete rect_potential_solver;
}

void Rect2DDetector::compute_potential() {
	rect_potential_solver->run();
}

