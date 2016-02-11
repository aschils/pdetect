/*
 * SerratedRect2DBoundaryCond.hpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */

#pragma once

#include "BoundaryConditions.hpp"
#include "SerratedRect2DBoundaryValues.hpp"
#include "../Utils.hpp"

using namespace dealii;

template<unsigned dim>
class SerratedRect2DBoundaryCond: public BoundaryConditions<dim> {

public:

	SerratedRect2DBoundaryCond() {
	}

	SerratedRect2DBoundaryCond(SerratedRectGeoInfo *geo_info,
			double strip_potential) {

		set_class_var(geo_info);
		this->values = new SerratedRect2DBoundaryValues<dim>(
				geo_info->get_nbr_of_strips(),
				strip_potential, geo_info);
	}

	void set_periodicity_constraints(
			typename Triangulation<dim>::active_cell_iterator &cell) {

		for (unsigned f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (cell->face(f)->at_boundary()) {
				if (Utils::equals_double(cell->face(f)->center()[0], 0.0,
						0.000001)
						|| Utils::equals_double(cell->face(f)->center()[0],
								geo_info->get_length(), 0.000001)) {
					cell->face(f)->set_boundary_id(1);
				} else if (Utils::equals_double(cell->face(f)->center()[1],
						geo_info->get_width(), 0.000001)) {

					if (!(geo_info->get_strip_width() == 0
							&& geo_info->is_strip<dim>(
									cell->face(f)->center())))
						cell->face(f)->set_boundary_id(1);
				}
			}
		}
	}

protected:

	void set_class_var(SerratedRectGeoInfo *geo_info) {
		this->geo_info = geo_info;
	}

private:
	SerratedRectGeoInfo *geo_info;
};
