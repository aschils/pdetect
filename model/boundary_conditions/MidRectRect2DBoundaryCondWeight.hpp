/*
 * MidRectRect2DBoundaryCondWeight.hpp
 *
 *  Created on: 3 févr. 2016
 *      Author: aschils
 */



#pragma once

#include "MidRectRect2DBoundaryCond.hpp"
#include "MidRectRect2DBoundaryValuesWeight.hpp"

template <unsigned dim>
class MidRectRect2DBoundaryCondWeight: public MidRectRect2DBoundaryCond<dim> {

public:
	MidRectRect2DBoundaryCondWeight(unsigned half_width,
			unsigned half_strip_width,
			unsigned strip_length,
			unsigned half_inter_strip_dist,
			unsigned nbr_of_strips,
			double potential){
		this->set_class_var(nbr_of_strips, half_inter_strip_dist,
				strip_length);
		this->values = new MidRectRect2DBoundaryValuesWeight<dim>(half_width,
				half_strip_width,
				strip_length,
				half_inter_strip_dist,
				nbr_of_strips,
				potential);
	}
};
