/*
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#pragma once

#include <deal.II/grid/tria.h>
#include <deal.II/base/function.h>

#include <math.h>

#include "../errors.hpp"
#include "../geometry_info/SerratedRectGeoInfo.hpp"
#include "../Utils.hpp"

using namespace dealii;

template<unsigned dim>
class SerratedRect2DBoundaryValues: public Function<dim> {
public:

	SerratedRect2DBoundaryValues(unsigned nbr_of_strips,
			double strip_potential,
			SerratedRectGeoInfo *geo_info);

	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;

protected:
	unsigned nbr_of_strips;
	double strip_potential;
	SerratedRectGeoInfo *geo_info;
};

template<unsigned dim>
SerratedRect2DBoundaryValues<dim>::SerratedRect2DBoundaryValues(
		unsigned nbr_of_strips,
		double strip_potential,
		SerratedRectGeoInfo *geo_info) {
	this->nbr_of_strips = nbr_of_strips;
	this->strip_potential = strip_potential;
	this->geo_info = geo_info;
}

template<unsigned dim>
double SerratedRect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {
	if (nbr_of_strips > 0
			&& geo_info->is_strip<dim>(p))
		return strip_potential;
	else
		return 0.0;
}

