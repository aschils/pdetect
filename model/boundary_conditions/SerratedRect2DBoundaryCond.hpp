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

	SerratedRect2DBoundaryCond(unsigned nbr_of_strips, unsigned rect_length,
			unsigned rect_width, double strip_potential, unsigned half_pitch,
			unsigned strip_length, unsigned strip_width) {
		this->rect_length = rect_length;
		this->rect_width = rect_width;
		this->strip_width = strip_width;
		this->strip_length = strip_length;
		this->periodic_str_length = strip_length + 2*half_pitch;
		this->half_pitch = half_pitch;
		this->values = new SerratedRect2DBoundaryValues<dim>(nbr_of_strips,
				rect_length, rect_width, strip_potential, half_pitch, strip_length,
				strip_width);
	}

	void set_periodicity_constraints(
			typename Triangulation<dim>::active_cell_iterator &cell) {

		for (unsigned f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {
			if (cell->face(f)->at_boundary()) {
				if (Utils::equals_double(cell->face(f)->center()[0], 0.0,
						0.000001)
						|| Utils::equals_double(cell->face(f)->center()[0],
								rect_length, 0.000001)) {
					cell->face(f)->set_boundary_id(1);
				} else if (Utils::equals_double(cell->face(f)->center()[1],
						rect_width, 0.000001)) {

					if (!(strip_width == 0
							&& SerratedRect2DBoundaryValues<dim>::is_strip(
									cell->face(f)->center(), rect_width,
									strip_width, strip_length, half_pitch)))
						cell->face(f)->set_boundary_id(1);
				}
			}
		}
	}

protected:
	unsigned rect_length, rect_width;
	unsigned strip_length, strip_width, periodic_str_length, half_pitch;
}
;
