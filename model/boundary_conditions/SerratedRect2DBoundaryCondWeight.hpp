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
	SerratedRect2DBoundaryCondWeight(SerratedRectGeoInfo *geo_info,
			unsigned nbr_of_strips, double strip_potential, unsigned half_pitch,
			unsigned strip_length, unsigned strip_width){
		this->set_class_var(geo_info, half_pitch, strip_length, strip_width);
		this->values = new SerratedRect2DBoundaryValuesWeight<dim>(nbr_of_strips,
				strip_potential, geo_info);
	}
};
