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

	typedef bg::model::point<double, 2, bg::cs::cartesian> bpoint;

	SerratedRect2DBoundaryValuesWeight(unsigned nbr_of_strips,
			double strip_potential,
			SerratedRectGeoInfo *geo_info) :
			SerratedRect2DBoundaryValues<dim>(nbr_of_strips, strip_potential,
					geo_info) {
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		bpoint bp = Utils::dealii_point_to_bpoint<dim>(p);

		bool (SerratedRectGeoInfo::*func)(const bpoint &bp);
		func = &SerratedRectGeoInfo::is_middle_strip<dim>;

		if ((this->geo_info->*func)(bp)) {
			return this->strip_potential;
		} else
			return 0.0;
	}
};
