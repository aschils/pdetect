/*
 * CirclePotential2DBoundaryValues.hpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */

#pragma once

#include <deal.II/grid/tria.h>
#include <deal.II/base/function.h>
#include "../Utils.hpp"

using namespace dealii;

template<unsigned dim>
class MidCircleRect2DBoundaryValues: public Function<dim> {

public:
	MidCircleRect2DBoundaryValues(unsigned half_width, double potential){
		this->half_width = half_width;
		this->potential = potential;
	}

	double value(const Point<dim> &p,
			const unsigned int /*component*/) const {

		double y = p[1];
		double epsilon = 0.000001;

		if(Utils::equals_double(y, -half_width, epsilon) ||
				Utils::equals_double(y, half_width, epsilon))
			return 0.0;
		else
			return potential;
	}

private:
	int half_width;
	double potential;
};
