/*
 *
 *
 *  Created on: 29 janv. 2016
 *      Author: aschils
 */

#pragma once

#include "../boundary_conditions/MidCircleRect2DBoundaryCond.hpp"
#include "../boundary_conditions/MidCircleRect2DBoundaryCondWeight.hpp"
#include "../geometry_info/MidCircleRectGeoInfo.hpp"
#include "Detector2D.hpp"

class MidCircleRect2DDetector: public Detector2D {

public:

	MidCircleRect2DDetector(unsigned half_width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned half_inter_potential_srcs_dist, double potential,
			double refine_accuracy, unsigned max_iter,
			double stop_accuracy, unsigned material_id);

	MidCircleRect2DDetector(unsigned half_width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned half_inter_potential_srcs_dist, double potential,
			unsigned material_id);

	std::string params_to_string();

private:
	unsigned half_width;
	unsigned nbr_of_potential_src;
	unsigned potential_src_radius;
	unsigned half_inter_potential_srcs_dist;
};

