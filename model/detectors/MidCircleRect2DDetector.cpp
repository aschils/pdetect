/*
 *
 *
 *  Created on: 29 janv. 2016
 *      Author: aschils
 */

#include "MidCircleRect2DDetector.hpp"

MidCircleRect2DDetector::MidCircleRect2DDetector(unsigned half_width,
		unsigned nbr_of_potential_src, unsigned potential_src_radius,
		unsigned half_inter_potential_srcs_dist, double potential,
		unsigned material_id) :
		MidCircleRect2DDetector(half_width, nbr_of_potential_src,
				potential_src_radius, half_inter_potential_srcs_dist, potential,
				0.01, 10000, 10e-12, material_id) {
}

MidCircleRect2DDetector::MidCircleRect2DDetector(unsigned half_width,
		unsigned nbr_of_potential_src, unsigned potential_src_radius,
		unsigned half_inter_potential_srcs_dist, double potential,
		double refine_accuracy, unsigned max_iter, double stop_accuracy,
		unsigned material_id):

		Detector2D(max_iter, potential, stop_accuracy, refine_accuracy,
						material_id){

	this->half_width = half_width;
	this->nbr_of_potential_src = nbr_of_potential_src;
	this->potential_src_radius = potential_src_radius;
	this->half_inter_potential_srcs_dist = half_inter_potential_srcs_dist;
	this->geo_info = new MidCircleRectGeoInfo();

	MyGridGenerator<2>::rectangle_with_circular_holes(*triangulation,
			half_width, potential_src_radius, half_inter_potential_srcs_dist,
			nbr_of_potential_src);

	boundary_conditions = new MidCircleRect2DBoundaryCond<2>(half_width,
			potential, nbr_of_potential_src,
			half_inter_potential_srcs_dist);
	potential_solver = new LaplaceSolver<2>(triangulation, refine_accuracy,
			max_iter, stop_accuracy, zero_right_hand_side, boundary_conditions,
			true);

	MyGridGenerator<2>::rectangle_with_circular_holes(*triangulation_weight,
			half_width, potential_src_radius, half_inter_potential_srcs_dist,
			nbr_of_potential_src);
	boundary_conditions_weight = new MidCircleRect2DBoundaryCondWeight<2>(
			half_width, WEIGHT_STRIP_POTENTIAL, nbr_of_potential_src, potential_src_radius,
			half_inter_potential_srcs_dist);
	potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
			refine_accuracy_weight, max_iter, stop_accuracy, zero_right_hand_side,
			boundary_conditions_weight, true);
}

std::string MidCircleRect2DDetector::params_to_string() {

	std::string str = "half_width" + std::to_string(half_width)
			+ "_nbr_of_potential_src_" + std::to_string(nbr_of_potential_src)
			+ "_potential_src_radius_" + std::to_string(potential_src_radius)
			+ "_half_inter_potential_srcs_dist_"
			+ std::to_string(half_inter_potential_srcs_dist) + "_potential_"
			+ std::to_string(strip_potential) + "refine_accuracy_"
			+ std::to_string(refine_accuracy) + "_max_iter_"
			+ std::to_string(max_iter) + "_stop_accuracy_"
			+ std::to_string(stop_accuracy);
	return str;
}

