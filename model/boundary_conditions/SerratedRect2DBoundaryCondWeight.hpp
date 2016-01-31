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
	SerratedRect2DBoundaryCondWeight(unsigned nbr_of_strips, unsigned rect_length,
			unsigned rect_width, double strip_potential, unsigned pitch,
			unsigned strip_length, unsigned strip_width) {
		this->rect_length = rect_length;
		this->rect_width = rect_width;
		this->values = new SerratedRect2DBoundaryValuesWeight<dim>(nbr_of_strips,
				rect_length, rect_width, strip_potential, pitch,
				strip_length, strip_width);
	}
};
