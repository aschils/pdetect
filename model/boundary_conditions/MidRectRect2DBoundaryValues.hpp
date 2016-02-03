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
	MidRectRect2DBoundaryValues(double potential, unsigned half_strip_width,
			unsigned half_width) {
		this->potential = potential;
		this->width = 2*half_width;
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		double y = p[1];
		double epsilon = 0.000001;

		if (Utils::equals_double(y, 0.0, epsilon)
				|| Utils::equals_double(y, width, epsilon)) {
			return 0.0;
		} else
			return potential;
	}

private:
	double potential;
	unsigned width;
};
