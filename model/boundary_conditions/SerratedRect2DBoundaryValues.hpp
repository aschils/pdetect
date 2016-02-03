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
#include "../Utils.hpp"

using namespace dealii;

template<unsigned dim>
class SerratedRect2DBoundaryValues: public Function<dim> {
public:

	SerratedRect2DBoundaryValues(unsigned nbr_of_strips, unsigned rect_length,
			unsigned rect_width, double strip_potential, unsigned half_pitch,
			unsigned strip_length, unsigned strip_width);

	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;

	static bool is_strip(const Point<dim> &p, unsigned rect_width,
			unsigned strip_width, unsigned strip_length, unsigned half_pitch);
protected:
	unsigned nbr_of_strips = 1;
	unsigned rect_length = 1;
	unsigned rect_width = 1;
	double strip_potential = 1;

	unsigned strip_width, strip_length = 1;
	unsigned periodic_str_length, half_pitch = 1;
};

template<unsigned dim>
SerratedRect2DBoundaryValues<dim>::SerratedRect2DBoundaryValues(
		unsigned nbr_of_strips, unsigned rect_length, unsigned rect_width,
		double strip_potential, unsigned half_pitch, unsigned strip_length,
		unsigned strip_width) {
	this->rect_width = rect_width;
	this->rect_length = rect_length;
	this->strip_potential = strip_potential;
	this->strip_length = strip_length;
	this->strip_width = strip_width;
	this->nbr_of_strips = nbr_of_strips;
	this->periodic_str_length = 2*half_pitch + strip_length;
	this->half_pitch = half_pitch;
}

/**
 *
 * We define a periodic structure as a strip surrounded by half a pitch on
 * each side. It is the structure repeated on the detector upper surface to
 * implement a grid of strips/pixels. Examples:
 *
 * Detector with one periodic structure:
 * |  ---  |
 *
 * Detector with two periodic structures:
 * |  ---    ---  |
 */
template<unsigned dim>
bool SerratedRect2DBoundaryValues<dim>::is_strip(const Point<dim> &p,
		unsigned rect_width, unsigned strip_width, unsigned strip_length,
		unsigned half_pitch) {

	unsigned pitch = 2*half_pitch;
	unsigned periodic_str_length = pitch + strip_length;

	double epsilon = 0.00000001;

	double x = p[0];
	double y = p[1];

	double bottom_of_strip = rect_width - strip_width;

	if (!Utils::greater_than_or_equals_double(y, bottom_of_strip, epsilon)
			|| periodic_str_length == 0.0 || strip_length == 0.0)
		return false;

	unsigned nbr_of_prev_periodic_str = x / periodic_str_length;
	double delta_from_prev_periodic_str = x
			- nbr_of_prev_periodic_str * periodic_str_length;

	double strip_border_right = half_pitch + strip_length;

	return Utils::greater_than_or_equals_double(delta_from_prev_periodic_str,
			half_pitch, epsilon)
			&& Utils::less_than_or_equals_double(delta_from_prev_periodic_str,
					strip_border_right, epsilon);
}

template<unsigned dim>
double SerratedRect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {
	if (nbr_of_strips > 0
			&& is_strip(p, rect_width, strip_width, strip_length, half_pitch))
		return strip_potential;
	else
		return 0.0;
}

