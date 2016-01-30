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

	SerratedRect2DBoundaryValues(unsigned nbr_of_strips, double rect_length_fe,
			double rect_width_fe, double strip_potential, double pitch_fe,
			double strip_length_fe, double strip_width_fe);

	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;

protected:
	unsigned nbr_of_strips = 1;
	double rect_length_fe = 1.0;
	double rect_width_fe = 1.0;
	double strip_potential = 1.0;

	double strip_width_fe, strip_length_fe = 1.0;
	double periodic_str_length_fe, half_pitch_fe = 1.0;

	bool is_strip(const Point<dim> &p) const;
};

template<unsigned dim>
SerratedRect2DBoundaryValues<dim>::SerratedRect2DBoundaryValues(
		unsigned nbr_of_strips, double rect_length_fe, double rect_width_fe,
		double strip_potential, double pitch_fe, double strip_length_fe,
		double strip_width_fe) {
	this->rect_width_fe = rect_width_fe;
	this->rect_length_fe = rect_length_fe;
	this->strip_potential = strip_potential;
	this->strip_length_fe = strip_length_fe;
	this->strip_width_fe = strip_width_fe;
	this->nbr_of_strips = nbr_of_strips;
	this->periodic_str_length_fe = pitch_fe + strip_length_fe;
	this->half_pitch_fe = pitch_fe / 2;
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
bool SerratedRect2DBoundaryValues<dim>::is_strip(const Point<dim> &p) const {

	double epsilon = 0.00001;

	double x = p[0];
	double y = p[1];

	double bottom_of_strip = rect_width_fe - strip_width_fe;

	if (!Utils::greater_than_or_equals_double(y, bottom_of_strip, epsilon)
			|| periodic_str_length_fe == 0.0 || strip_length_fe == 0.0)
		return false;

	unsigned nbr_of_prev_periodic_str = x / periodic_str_length_fe;
	double delta_from_prev_periodic_str = x
			- nbr_of_prev_periodic_str * periodic_str_length_fe;

	double strip_border_right = half_pitch_fe + strip_length_fe;

	return Utils::greater_than_or_equals_double(delta_from_prev_periodic_str,
			half_pitch_fe, epsilon)
			&& Utils::less_than_or_equals_double(delta_from_prev_periodic_str,
					strip_border_right, epsilon);
}

template<unsigned dim>
double SerratedRect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {
	if (nbr_of_strips > 0 && is_strip(p))
		return strip_potential;
	else
		return 0.0;
}

