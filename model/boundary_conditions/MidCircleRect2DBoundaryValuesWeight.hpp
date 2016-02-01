/*
 * MidCircleRect2DBoundaryValuesWeight.hpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "MidCircleRect2DBoundaryValues.hpp"

template<unsigned dim>
class MidCircleRect2DBoundaryValuesWeight: public Function<dim> {

public:
	MidCircleRect2DBoundaryValuesWeight(unsigned width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned inter_potential_srcs_dist, double potential) {
		this->half_width = ceil(width / 2.0);
		this->potential = potential;
		this->potential_src_radius = potential_src_radius;

		unsigned half_inter_potential_srcs_dist = ceil(
				inter_potential_srcs_dist/2.0);
		unsigned nbr_of_src_before = nbr_of_potential_src / 2;
		middle_potential_src_left_x = half_inter_potential_srcs_dist-potential_src_radius
				+ nbr_of_src_before * inter_potential_srcs_dist;
		middle_potential_src_right_x = middle_potential_src_left_x
				+ 2 * potential_src_radius;
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		double y = p[1];
		double epsilon = 0.000001;

		/*if (Utils::equals_double(y, -half_width, epsilon)
				|| Utils::equals_double(y, half_width, epsilon)) {
			std::cout << "coucou1" << std::endl;
			return 0.0;
		}*/
		if(is_middle_strip(p)){
			return potential;
		}
		else{
			return 0.0;
		}
	}

private:
	unsigned half_width;
	double potential;
	int potential_src_radius;

	//Range of x values where the middle potential source remains
	unsigned middle_potential_src_left_x;
	unsigned middle_potential_src_right_x;

	bool is_middle_strip(const Point<dim> &p) const {

		double x = p[0];
		double y = p[1];

		double epsilon = 0.000000001;

		return
				Utils::greater_than_or_equals_double(x,
						middle_potential_src_left_x,
						epsilon)
				&& Utils::less_than_or_equals_double(x,
						middle_potential_src_right_x,
						epsilon)
				&& Utils::greater_than_or_equals_double(y,
						-potential_src_radius,
						epsilon)
				&& Utils::less_than_or_equals_double(y,
						potential_src_radius,
						epsilon);
	}
};
