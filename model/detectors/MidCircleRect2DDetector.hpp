/*
 *
 *
 *  Created on: 29 janv. 2016
 *      Author: aschils
 */

#pragma once

#include "../boundary_conditions/MidCircleRect2DBoundaryCond.hpp"
#include "../boundary_conditions/MidCircleRect2DBoundaryCondWeight.hpp"
#include "Detector2D.hpp"

class MidCircleRect2DDetector: public Detector2D {

public:

	MidCircleRect2DDetector(unsigned width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned inter_potential_srcs_dist, double potential,
			unsigned refine_level, unsigned max_iter,
			double stop_accuracy);

	MidCircleRect2DDetector(unsigned width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned inter_potential_srcs_dist, double potential);

	std::string params_to_string();

private:
	unsigned width;
	unsigned nbr_of_potential_src;
	unsigned potential_src_radius;
	unsigned inter_potential_srcs_dist;
};

