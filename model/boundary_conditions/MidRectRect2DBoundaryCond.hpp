/*
 * MidRectRect2DBoundaryCond.hpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "BoundaryConditions.hpp"
#include "MidRectRect2DBoundaryValues.hpp"
#include "../geometry_info/MidRectRectGeoInfo.hpp"

template<unsigned dim>
class MidRectRect2DBoundaryCond: public BoundaryConditions<dim> {

public:
	MidRectRect2DBoundaryCond(MidRectRectGeoInfo *geo_info,
			double potential) {
		set_class_var(geo_info->get_nbr_of_strips(),
					  geo_info->get_half_inter_strip_dist(),
					  geo_info->get_strip_length());
		this->values = new MidRectRect2DBoundaryValues<dim>(potential,
						geo_info->get_half_strip_width(),
						geo_info->get_half_width());
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

	MidRectRect2DBoundaryCond(){}

	void set_class_var(unsigned nbr_of_strips,
			unsigned half_inter_potential_srcs_dist, unsigned strip_length) {
		left_boundary_x = 0.0;
		right_boundary_x = nbr_of_strips
				* (strip_length + 2 * half_inter_potential_srcs_dist);
	}

private:
	unsigned left_boundary_x, right_boundary_x;
};
