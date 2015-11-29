/*
 * Solution.hpp
 *
 *  Created on: 29 nov. 2015
 *      Author: aschils
 */

#ifndef __SOLUTION_HPP__
#define __SOLUTION_HPP__

#include <deal.II/numerics/data_out_dof_data.h>

using namespace dealii;

template <int dim>
class Solution{

public:

	Solution(Vector<double> data, 	DoFHandler<dim> *dof_handler){
		this->data = data;
		this->dof_handler = dof_handler;
	}

	Vector<double> data;
	DoFHandler<dim> *dof_handler;
};

#endif
