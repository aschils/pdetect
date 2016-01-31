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
			unsigned rect_length, unsigned rect_width, double strip_potential,
			unsigned pitch, unsigned strip_length, unsigned strip_width) :
			SerratedRect2DBoundaryValues<dim>(nbr_of_strips, rect_length,
					rect_width, strip_potential, pitch, strip_length,
					strip_width) {
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		unsigned nbr_of_prev_periodic_str = p[0] / this->periodic_str_length;

		if (this->nbr_of_strips > 0 && this->is_strip(p)
				&& nbr_of_prev_periodic_str == this->nbr_of_strips / 2) {
			return this->strip_potential;
		} else
			return 0.0;
	}
};
