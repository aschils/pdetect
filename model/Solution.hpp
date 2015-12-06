/*
 * Solution.hpp
 *
 *  Created on: 29 nov. 2015
 *      Author: aschils
 */

#ifndef __SOLUTION_HPP__
#define __SOLUTION_HPP__

#include <fstream>
#include <iostream>
#include <tuple>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>

#include "Utils.hpp"
#include "VectorUtils.hpp"

using namespace dealii;

template<unsigned dim>
class SolutionData {
public:
	double fun;
	Tensor<1, dim> gradient;
	Tensor<2, dim> hessian;

	SolutionData(double fun, Tensor<1, dim> gradient, Tensor<2, dim> hessian) {
		this->fun = fun;
		this->gradient = gradient;
		this->hessian = hessian;
	}

	SolutionData() {
	}
};

template<int dim>
class Solution {

public:

	std::vector<std::pair<std::vector<double>, SolutionData<dim> > > coord_and_data;

	Solution() {
	}

	void set_fun_drawer(DataOut<dim> fun_drawer){
		this->fun_drawer = fun_drawer;
	}

	void set_derivatives_drawer(DataOut<dim> derivatives_drawer){
		this->derivatives_drawer = derivatives_drawer;
	}

	void draw_vtk_graph_fun(std::string output_file) {
		std::ofstream output(output_file);
		fun_drawer.write_vtk(output);
	}

	void draw_vtk_graph_derivatives(std::string output_file) {
		std::ofstream output(output_file);
		derivatives_drawer.write_vtk(output);
	}

	void sort_by_coord() {
		VectorUtils::sort_by_coord<dim, SolutionData<dim> >(coord_and_data);
	}

private:
	//These two structures contain data already available in coord_and_data,
	//but it is useful to keep them as such to easily output vtk graph file
	//using deal.ii DataOut class.
	DataOut<dim> fun_drawer;
	DataOut<dim> derivatives_drawer;
};

#endif
