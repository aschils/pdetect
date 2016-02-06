/*
 * CirclePotential2DBoundaryCond.hpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */
#pragma once

#include "BoundaryConditions.hpp"
#include "MidCircleRect2DBoundaryValues.hpp"

template<unsigned dim>
class MidCircleRect2DBoundaryCond: public BoundaryConditions<dim> {

public:

	MidCircleRect2DBoundaryCond(unsigned half_width, double potential,
			unsigned nbr_of_potential_src,
			unsigned half_inter_potential_srcs_dist) {
		set_class_var(nbr_of_potential_src, half_inter_potential_srcs_dist);
		this->values = new MidCircleRect2DBoundaryValues<dim>(half_width,
				potential);
	}

	void set_periodicity_constraints(
			typename Triangulation<dim>::active_cell_iterator &cell) {

		double epsilon = 0.000001;

		for (unsigned f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (cell->face(f)->at_boundary()) {
				if (Utils::equals_double(cell->face(f)->center()[0],
						left_boundary_x, epsilon)
						|| Utils::equals_double(cell->face(f)->center()[0],
								right_boundary_x, epsilon)) {
					cell->face(f)->set_boundary_id(1);
				}
			}
		}
	}

protected:

	MidCircleRect2DBoundaryCond(){}

	void set_class_var(unsigned nbr_of_potential_src,
			unsigned half_inter_potential_srcs_dist) {

		//Compute domain left and right boundaries x coordinate
		unsigned periodic_struct_length = 2 * half_inter_potential_srcs_dist;
		left_boundary_x = 0.0;
		right_boundary_x = periodic_struct_length * nbr_of_potential_src;
	}

private:
	int left_boundary_x, right_boundary_x;
};
