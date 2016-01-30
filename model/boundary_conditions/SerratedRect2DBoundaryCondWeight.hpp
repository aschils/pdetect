/*
 * SerratedRect2DBoundaryCondWeight.hpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */

#pragma once

#include "SerratedRect2DBoundaryCond.hpp"
#include "SerratedRect2DBoundaryValuesWeight.hpp"

template <unsigned dim>
class SerratedRect2DBoundaryCondWeight : public SerratedRect2DBoundaryCond<dim> {

public:
	SerratedRect2DBoundaryCondWeight(unsigned nbr_of_strips, double rect_length_fe,
			double rect_width_fe, double strip_potential, double pitch_fe,
			double strip_length_fe, double strip_width_fe) {
		this->rect_length_fe = rect_length_fe;
		this->rect_width_fe = rect_width_fe;
		this->values = new SerratedRect2DBoundaryValuesWeight<dim>(nbr_of_strips,
				rect_length_fe, rect_width_fe, strip_potential, pitch_fe,
				strip_length_fe, strip_width_fe);
	}
};
