/*
 * MidRectRect2DBoundaryValuesWeight.hpp
 *
 *  Created on: 3 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <deal.II/base/function.h>
#include "../Utils.hpp"

using namespace dealii;

template<unsigned dim>
class MidRectRect2DBoundaryValuesWeight: public Function<dim> {

public:

	MidRectRect2DBoundaryValuesWeight(MidRectRectGeoInfo *geo_info,
			double potential) {

		this->potential = potential;
		this->geo_info = geo_info;

	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		if (geo_info->is_middle_strip(p))
			return potential;
		else
			return 0.0;
	}

private:
	MidRectRectGeoInfo *geo_info;
	double potential;

};

