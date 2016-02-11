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
			 double strip_potential){
		this->set_class_var(geo_info);
		this->values = new SerratedRect2DBoundaryValuesWeight<dim>(
				geo_info->get_nbr_of_strips(),
				strip_potential, geo_info);
	}
};
