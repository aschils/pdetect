/*
 * MidRectRect2DBoundaryValuesWeight.hpp
 *
 *  Created on: 3 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <deal.II/base/function.h>

using namespace dealii;

template<unsigned dim>
class MidRectRect2DBoundaryValuesWeight: public Function<dim> {

public:
	MidRectRect2DBoundaryValuesWeight(unsigned half_width,
			unsigned half_strip_width,
			unsigned strip_length,
			unsigned half_inter_strip_dist,
			unsigned nbr_of_strips,
			double potential) {

		this->potential = potential;

		unsigned nbr_of_strips_before = nbr_of_strips/2;
		middle_strip_left_boundary_x = nbr_of_strips_before*(strip_length
				+2*half_inter_strip_dist)+half_inter_strip_dist;
		middle_strip_right_boundary_x = middle_strip_left_boundary_x+
				strip_length;
		middle_strip_bottom_boundary_y = half_width-half_strip_width;
		middle_strip_top_boundary_y = half_width+half_strip_width;
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		if (is_middle_strip(p))
			return potential;
		else
			return 0.0;
	}

private:

	double potential;

	int middle_strip_left_boundary_x, middle_strip_right_boundary_x;
	int middle_strip_bottom_boundary_y, middle_strip_top_boundary_y;

	bool is_middle_strip(const Point<dim> &p) const {
		double x = p[0];
		double y = p[1];
		return x >= middle_strip_left_boundary_x
				&& x <= middle_strip_right_boundary_x
				&& y <= middle_strip_top_boundary_y
				&& y >= middle_strip_bottom_boundary_y;
	}

};

