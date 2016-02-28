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

	this->strip_potential = potential;
	this->refine_level = refine_level;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;

	this->mrr_geo_info = new MidRectRectGeoInfo(half_width, half_strip_width,
			strip_length, nbr_of_strips, half_inter_strip_dist);
	this->geo_info = mrr_geo_info;

	MyGridGenerator<2>::rectangle_with_rectangular_holes(*triangulation,
			half_width, strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips);

	boundary_conditions = new MidRectRect2DBoundaryCond<2>(mrr_geo_info,
			potential);
	potential_solver = new LaplaceSolver<2>(triangulation, refine_level,
			max_iter, stop_accuracy, zero_right_hand_side, boundary_conditions,
			true);

	boundary_conditions_weight = new MidRectRect2DBoundaryCondWeight<2>(mrr_geo_info,
			weight_strip_potential);
	MyGridGenerator<2>::rectangle_with_rectangular_holes(*triangulation_weight,
			half_width, strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips);
	potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
			refine_level, max_iter, stop_accuracy, zero_right_hand_side,
			boundary_conditions_weight, true);
}

std::string MidRectRect2DDetector::params_to_string() {



	std::string str = "width" + std::to_string(geo_info->get_width())
			+ "_nbr_of_strips_" + std::to_string(geo_info->get_nbr_of_strips())
			+ "_strip_length_" + std::to_string(mrr_geo_info->get_strip_length())
			+ "_strip_width_" + std::to_string(mrr_geo_info->get_strip_width())
			+ "_half_inter_strip_dist_" +
			std::to_string(mrr_geo_info->get_half_inter_strip_dist())
			+ "_potential_" + std::to_string(strip_potential) + "_refine_level_"
			+ std::to_string(refine_level) + "_max_iter_"
			+ std::to_string(max_iter) + "_stop_accuracy_"
			+ std::to_string(stop_accuracy);
	return str;
}
