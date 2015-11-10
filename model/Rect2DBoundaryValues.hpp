/*
 * BoundaryValues.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#ifndef __BOUNDARY_VALUES__HPP__
#define __BOUNDARY_VALUES__HPP__

#include <deal.II/grid/tria.h>
#include <deal.II/base/function.h>

#include <math.h>

#include "errors.hpp"

using namespace dealii;

template<int dim>
class Rect2DBoundaryValues: public Function<dim> {
public:

	/**
	 * @pre:
	 * - rect_length_fe is the length of the rectangular finite element
	 * mathematical domain.
	 * - rect_length_fe > 0.0
	 * - strip_width and pitch are values in the domain language
	 * (i.e. microm,...).
	 * TODO complete this @pre
	 */
	Rect2DBoundaryValues(unsigned nbr_of_strip, unsigned strip_length,
			unsigned pitch, double rect_length_fe);
	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;

private:
	double strip_length_fe = 1.0;
	double pitch_length_fe = 1.0;
	double strip_pitch_pair_length_fe = 1.0;
	double rect_length_fe = 1.0;

	bool is_strip(double length_component) const;
};

#include "Rect2DBoundaryValues.cpp"

#endif
