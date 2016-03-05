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
		double potential, double refine_accuracy, unsigned max_iter,
		double stop_accuracy, unsigned material_id):
		Detector2D(max_iter, strip_potential, stop_accuracy, refine_accuracy,
				material_id) {

	this->mrr_geo_info = new MidRectRectGeoInfo(half_width, half_strip_width,
			strip_length, nbr_of_strips, half_inter_strip_dist);
	this->geo_info = mrr_geo_info;

	MyGridGenerator<2>::rectangle_with_rectangular_holes(*triangulation,
			half_width, strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips);

	boundary_conditions = new MidRectRect2DBoundaryCond<2>(mrr_geo_info,
			potential);
	potential_solver = new LaplaceSolver<2>(triangulation, refine_accuracy,
			max_iter, stop_accuracy, zero_right_hand_side, boundary_conditions,
			true);

	boundary_conditions_weight = new MidRectRect2DBoundaryCondWeight<2>(mrr_geo_info,
			WEIGHT_STRIP_POTENTIAL);
	MyGridGenerator<2>::rectangle_with_rectangular_holes(*triangulation_weight,
			half_width, strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips);
	potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
			refine_accuracy_weight, max_iter, stop_accuracy, zero_right_hand_side,
			boundary_conditions_weight, true);
}

std::string MidRectRect2DDetector::params_to_string() {



	std::string str = "width" + std::to_string(geo_info->get_width())
			+ "_nbr_of_strips_" + std::to_string(geo_info->get_nbr_of_strips())
			+ "_strip_length_" + std::to_string(mrr_geo_info->get_strip_length())
			+ "_strip_width_" + std::to_string(mrr_geo_info->get_strip_width())
			+ "_half_inter_strip_dist_" +
			std::to_string(mrr_geo_info->get_half_inter_strip_dist())
			+ "_potential_" + std::to_string(strip_potential) + "refine_accuracy_"
			+ std::to_string(refine_accuracy) + "_max_iter_"
			+ std::to_string(max_iter) + "_stop_accuracy_"
			+ std::to_string(stop_accuracy);
	return str;
}
