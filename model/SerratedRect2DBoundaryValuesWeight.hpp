/*
 * SerratedRect2DBoundaryValuesWeight.hpp
 *
 *  Created on: 28 nov. 2015
 *      Author: aschils
 */

#include "SerratedRect2DBoundaryValues.hpp"

#ifndef __SERRATED_BOUNDARY_VALUES_WEIGHT_HPP__
#define __SERRATED_BOUNDARY_VALUES_WEIGHT_HPP__

template<int dim>
class SerratedRect2DBoundaryValuesWeight: public SerratedRect2DBoundaryValues<
		dim> {

public:
	SerratedRect2DBoundaryValuesWeight(unsigned nbr_of_strips,
			double rect_length_fe, double rect_width_fe, double strip_potential,
			double pitch_fe, double strip_length_fe, double strip_width_fe) :
			SerratedRect2DBoundaryValues<dim>(nbr_of_strips, rect_length_fe,
					rect_width_fe, strip_potential, pitch_fe, strip_length_fe,
					strip_width_fe) {
	}

	double value(const Point<dim> &p, const unsigned int /*component*/) const {

		unsigned nbr_of_prev_periodic_str = p[0] / this->periodic_str_length_fe;

		if (this->nbr_of_strips > 0 && this->is_strip(p)
				&& nbr_of_prev_periodic_str == this->nbr_of_strips / 2) {
			return this->strip_potential;
		} else
			return 0.0;
	}
};

#endif
