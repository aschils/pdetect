/*
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#ifndef __SERRATED_BOUNDARY_VALUES_HPP__
#define __SERRATED_BOUNDARY_VALUES_HPP__

#include <deal.II/grid/tria.h>
#include <deal.II/base/function.h>

#include <math.h>

#include "errors.hpp"
#include "Utils.hpp"

using namespace dealii;

template<int dim>
class SerratedRect2DBoundaryValues: public Function<dim> {
public:

	SerratedRect2DBoundaryValues(unsigned nbr_of_strips, double rect_length_fe,
			double rect_width_fe, double strip_potential, double pitch_fe,
			double strip_length_fe, double strip_width_fe);

	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;

private:
	unsigned nbr_of_strips = 1;
	double rect_length_fe = 1.0;
	double rect_width_fe = 1.0;
	double strip_potential = 1.0;

	double strip_width_fe, strip_length_fe = 1.0;
	double periodic_str_length_fe, half_pitch_fe = 1.0;

	bool is_strip(const Point<dim> &p) const;
};

#include "SerratedRect2DBoundaryValues.cpp"

#endif

