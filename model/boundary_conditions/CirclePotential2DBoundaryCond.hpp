/*
 * CirclePotential2DBoundaryCond.hpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */
#pragma once

#include "CirclePotential2DBoundaryValues.hpp"

template<unsigned dim>
class CirclePotential2DBoundaryCond: public BoundaryConditions<dim> {

public:

	CirclePotential2DBoundaryCond(double width, unsigned potential,
			unsigned nbr_of_potential_src, unsigned inter_potential_srcs_dist) {
		this->values = new CirclePotential2DBoundaryValues<dim>(width,
				potential);
		unsigned half_periodic_struct_length = ceil(inter_potential_srcs_dist / 2.0);
		left_boundary_x = -half_periodic_struct_length;
		right_boundary_x = (nbr_of_potential_src == 1)? half_periodic_struct_length:
				half_periodic_struct_length*(2*nbr_of_potential_src-1);
	}

	void set_periodicity_constraints(
			typename Triangulation<dim>::active_cell_iterator &cell) {

		double epsilon = 0.000001;

		for (unsigned f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (cell->face(f)->at_boundary()) {
				if (Utils::equals_double(cell->face(f)->center()[0], left_boundary_x,
						epsilon)
						|| Utils::equals_double(cell->face(f)->center()[0],
								right_boundary_x, epsilon)) {
					cell->face(f)->set_boundary_id(1);
				}
			}
		}
	}

private:
	int left_boundary_x, right_boundary_x;

};
