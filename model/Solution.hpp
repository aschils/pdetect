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

	//These two structures contain data already available in coord_and_data,
	//but it is useful to keep them as such to easily output vtk graph file
	//using deal.ii DataOut class.
	DataOut<dim> fun_drawer;
	DataOut<dim> derivatives_drawer;
	std::vector<std::pair<std::vector<double>, SolutionData<dim> > > coord_and_data;

	Solution() {
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

	void print() {

		for (unsigned i = 0; i < coord_and_data.size(); i++) {
			std::cout << "[coord: (";
			VectorUtils::print_vec_components(coord_and_data[i].first);
			std::cout << "), function: " << coord_and_data[i].second.fun;
			std::cout << " gradient: (";
			VectorUtils::print_tensor_components<dim>(
					coord_and_data[i].second.gradient);
			std::cout << ")] ";
		}
		std::cout << std::endl;
	}
};

#endif
