/*
 * BoundaryConditions.hpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */

#pragma once

#include <deal.II/base/function.h>

using namespace dealii;

template <unsigned dim>
class BoundaryConditions {

public:

	virtual void set_periodicity_constraints(
			typename Triangulation<dim>::active_cell_iterator &cell) = 0;

	Function<dim>* get_values(){
		return this->values;
	}

	virtual ~BoundaryConditions(){
		delete values;
	}

protected:
	Function<dim> *values;
};