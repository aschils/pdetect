/*
 * MidCircleRect2DBoundaryCondWeight.hpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "MidCircleRect2DBoundaryCond.hpp"
#include "MidCircleRect2DBoundaryValuesWeight.hpp"

template <unsigned dim>
class MidCircleRect2DBoundaryCondWeight: public MidCircleRect2DBoundaryCond<dim> {

public:

	MidCircleRect2DBoundaryCondWeight(unsigned width,
			unsigned nbr_of_potential_src, unsigned potential_src_radius,
			unsigned inter_potential_srcs_dist, double potential):
				MidCircleRect2DBoundaryCond<dim>(width, potential,
							nbr_of_potential_src, inter_potential_srcs_dist){
		this->values = new MidCircleRect2DBoundaryValuesWeight<dim>(width,
				nbr_of_potential_src, potential_src_radius,
				inter_potential_srcs_dist, potential);
	}

};
