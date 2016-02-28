/*
 * MidRectRect2DDetector.cpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#include "MidRectRect2DDetector.hpp"

MidRectRect2DDetector::MidRectRect2DDetector(unsigned half_width,
		unsigned strip_length, unsigned half_strip_width,
		unsigned half_inter_strip_dist, unsigned nbr_of_strips,
		double potential, unsigned refine_level, unsigned max_iter,
		double stop_accuracy) {
	this->half_width = half_width;
	this->strip_length = strip_length;
	this->half_strip_width = half_strip_width;
	this->half_inter_strip_dist = half_inter_strip_dist;
	this->nbr_of_strips = nbr_of_strips;
	this->strip_potential = potential;
	this->refine_level = refine_level;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;

	this->geo_info = new MidRectRectGeoInfo(half_width, half_strip_width,
			strip_length, nbr_of_strips, half_inter_strip_dist);

	MyGridGenerator<2>::rectangle_with_rectangular_holes(*triangulation,
			half_width, strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips);

	boundary_conditions = new MidRectRect2DBoundaryCond<2>(half_width,
			potential, nbr_of_strips, half_inter_strip_dist, strip_length,
			half_strip_width);
	potential_solver = new LaplaceSolver<2>(triangulation, refine_level,
			max_iter, stop_accuracy, zero_right_hand_side, boundary_conditions,
			true);

	boundary_conditions_weight = new MidRectRect2DBoundaryCondWeight<2>(
			half_width, half_strip_width, strip_length, half_inter_strip_dist,
			nbr_of_strips, weight_strip_potential);
	MyGridGenerator<2>::rectangle_with_rectangular_holes(*triangulation_weight,
			half_width, strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips);
	potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
			refine_level, max_iter, stop_accuracy, zero_right_hand_side,
			boundary_conditions_weight, true);
}

std::string MidRectRect2DDetector::params_to_string() {
	std::string str = "half_width" + std::to_string(half_width)
			+ "_nbr_of_strips_" + std::to_string(nbr_of_strips)
			+ "_strip_length_" + std::to_string(strip_length)
			+ "_half_strip_width_" + std::to_string(half_strip_width)
			+ "_half_inter_strip_dist_" + std::to_string(half_inter_strip_dist)
			+ "_potential_" + std::to_string(strip_potential) + "_refine_level_"
			+ std::to_string(refine_level) + "_max_iter_"
			+ std::to_string(max_iter) + "_stop_accuracy_"
			+ std::to_string(stop_accuracy);
	return str;
}
