/*
 * MidRectRect2DBoundaryCondWeight.hpp
 *
 *  Created on: 3 févr. 2016
 *      Author: aschils
 */



#pragma once

#include "MidRectRect2DBoundaryCond.hpp"
#include "MidRectRect2DBoundaryValuesWeight.hpp"
#include "../geometry_info/MidRectRectGeoInfo.hpp"

template <unsigned dim>
class MidRectRect2DBoundaryCondWeight: public MidRectRect2DBoundaryCond<dim> {

public:
	MidRectRect2DBoundaryCondWeight(MidRectRectGeoInfo *geo_info,
			double potential){
		this->set_class_var(geo_info->get_nbr_of_strips(),
				geo_info->get_half_inter_strip_dist(),
				geo_info->get_strip_length());
		this->values = new MidRectRect2DBoundaryValuesWeight<dim>(geo_info,
				potential);
	}
};
