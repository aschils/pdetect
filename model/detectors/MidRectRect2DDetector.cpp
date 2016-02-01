/*
 * MidRectRect2DDetector.cpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#include "MidRectRect2DDetector.hpp"

MidRectRect2DDetector::MidRectRect2DDetector(unsigned width,
		unsigned strip_length, unsigned strip_width, unsigned inter_strip_dist,
		unsigned nbr_of_strips, double potential, unsigned refine_level,
		unsigned max_iter, double stop_accuracy) {
	this->width = width;
	this->strip_length = strip_length;
	this->strip_width = strip_width;
	this->inter_strip_dist = inter_strip_dist;
	this->nbr_of_strips = nbr_of_strips;
	this->strip_potential = potential;
	this->refine_level = refine_level;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;

	MyGridGenerator<2>::rectangle_with_circular_holes(*triangulation, width,
			potential, inter_strip_dist, nbr_of_strips);
//	boundary_conditions = new MidCircleRect2DBoundaryCond<2>(width, potential,
//			nbr_of_potential_src, inter_potential_srcs_dist);
//	potential_solver = new LaplaceSolver<2>(triangulation, refine_level,
//			max_iter, stop_accuracy, zero_right_hand_side, boundary_conditions,
//			true);
//
//	boundary_conditions_weight = new MidCircleRect2DBoundaryCond<2>(width,
//			potential, nbr_of_potential_src, inter_potential_srcs_dist);
//	MyGridGenerator<2>::rectangle_with_circular_holes(*triangulation_weight,
//			width, potential_src_radius, inter_potential_srcs_dist,
//			nbr_of_potential_src);
//	potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
//			refine_level, max_iter, stop_accuracy, zero_right_hand_side,
//			boundary_conditions_weight, true);
}

std::string MidRectRect2DDetector::params_to_string() {
	//TODO
	return "";
}
