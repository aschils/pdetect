/*
 *
 *
 *  Created on: 29 janv. 2016
 *      Author: aschils
 */

#pragma once

#include "Detector2D.hpp"
#include "../boundary_conditions/CirclePotential2DBoundaryCond.hpp"

class CirclePotential2DDet: public Detector2D {

public:

	CirclePotential2DDet(unsigned width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned inter_potential_srcs_dist, unsigned potential,
			unsigned refine_level, unsigned max_iter,
			double stop_accuracy);

	CirclePotential2DDet(unsigned width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned inter_potential_srcs_dist, unsigned potential);

	std::string params_to_string();

private:
	unsigned width;
	unsigned nbr_of_potential_src;
	unsigned potential_src_radius;
	unsigned inter_potential_srcs_dist;
	unsigned potential;

};

