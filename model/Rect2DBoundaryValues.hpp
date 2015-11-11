/*
 * BoundaryValues.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#ifndef __BOUNDARY_VALUES_HPP__
#define __BOUNDARY_VALUES_HPP__

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
	 * mathematical domain (i.e. the longest edge of the rectangle).
	 * - rect_length_fe > 0
	 * - strip_width and pitch are values in the domain language
	 * (i.e. micrometers,...).
	 *
	 * @throw:
	 * - ZERO_OR_NEGATIVE_DOMAIN_LENGTH_ERROR
	 *
	 */
	Rect2DBoundaryValues(unsigned nbr_of_strip, unsigned strip_length,
			unsigned pitch, double rect_length_fe, double rect_width_fe,
			double strip_potential);
	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;

private:
	double strip_length_fe = 1.0;
	double pitch_length_fe = 1.0;
	double strip_pitch_pair_length_fe = 1.0;
	double rect_length_fe = 1.0;
	double half_rect_width_fe = 1.0;
	double half_pitch_length_fe = 1.0;
	double strip_potential = 1.0;

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
	 *
	 * @pre:
	 * - if no strip or strip length is 0, half_pitch_length_fe must have been set
	 * to 0 in the constructor.
	 *
	 * @return: true if the there is a strip at position p, false otherwise.
	 */
	bool is_strip(const Point<dim> &p) const;
};

#include "Rect2DBoundaryValues.cpp"

#endif
