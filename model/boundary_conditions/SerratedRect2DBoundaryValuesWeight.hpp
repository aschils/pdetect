/*
 * SerratedRect2DBoundaryValuesWeight.hpp
 *
 *  Created on: 28 nov. 2015
 *      Author: aschils
 */

#pragma once

#include "SerratedRect2DBoundaryValues.hpp"

template<unsigned dim>
class SerratedRect2DBoundaryValuesWeight: public SerratedRect2DBoundaryValues<
		dim> {

public:
	SerratedRect2DBoundaryValuesWeight(unsigned nbr_of_strips,
			double strip_potential,
			SerratedRectGeoInfo *geo_info) :
			SerratedRect2DBoundaryValues<dim>(nbr_of_strips, strip_potential,
					geo_info) {
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		bool (SerratedRectGeoInfo::*func)(const Point<dim> &p);
		func = &SerratedRectGeoInfo::is_middle_strip<dim>;

		if ((this->geo_info->*func)(p)) {
			return this->strip_potential;
		} else
			return 0.0;
	}
};
