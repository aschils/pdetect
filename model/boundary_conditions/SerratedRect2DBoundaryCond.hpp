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

	SerratedRect2DBoundaryCond(){}

	SerratedRect2DBoundaryCond(unsigned nbr_of_strips, double rect_length_fe,
			double rect_width_fe, double strip_potential, double pitch_fe,
			double strip_length_fe, double strip_width_fe) {
		this->rect_length_fe = rect_length_fe;
		this->rect_width_fe = rect_width_fe;
		this->values = new SerratedRect2DBoundaryValues<dim>(nbr_of_strips,
				rect_length_fe, rect_width_fe, strip_potential, pitch_fe,
				strip_length_fe, strip_width_fe);
	}

	void set_periodicity_constraints(
			typename Triangulation<dim>::active_cell_iterator &cell) {

		for (unsigned f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (cell->face(f)->at_boundary()) {
				if (Utils::equals_double(cell->face(f)->center()[0], 0.0,
						0.000001)
						|| Utils::equals_double(cell->face(f)->center()[0],
								rect_length_fe, 0.000001)
						|| Utils::equals_double(cell->face(f)->center()[1],
								rect_width_fe, 0.000001)) {
					cell->face(f)->set_boundary_id(1);
				}
			}
		}
	}

protected:
	double rect_length_fe, rect_width_fe;
}
;
