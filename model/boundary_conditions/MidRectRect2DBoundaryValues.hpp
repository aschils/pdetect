/*
 * MidRectRect2DBoundaryValues.hpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <deal.II/base/function.h>
#include "../Utils.hpp"

using namespace dealii;

template<unsigned dim>
class MidRectRect2DBoundaryValues: public Function<dim> {

public:
	MidRectRect2DBoundaryValues(double potential, unsigned strip_width,
			unsigned width) {
		this->potential = potential;
		unsigned half_wo_hole_width = ceil((width-strip_width)/2.0);
		this->real_width = 2*half_wo_hole_width+strip_width;
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		double y = p[1];
		double epsilon = 0.000001;

		if (Utils::equals_double(y, 0.0, epsilon)
				|| Utils::equals_double(y, real_width, epsilon)) {
			return 0.0;
		} else
			return potential;
	}

private:
	double potential;
	unsigned real_width;
};
