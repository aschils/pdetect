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
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>

#include "Utils.hpp"

using namespace dealii;

class Solution {

public:

	virtual void draw_vtk_graph(std::string output_file) = 0;

	virtual ~Solution() {

	}
};

template<int dim>
class SolutionScalar: public Solution {

public:

	std::vector<std::vector<double> > coord;

	DoFHandler<dim> *dof_handler;
	Vector<double> data;

	SolutionScalar(Vector<double> data, DoFHandler<dim> *dof_handler) {
		this->dof_handler = dof_handler;
		this->data = data;
	}

	void draw_vtk_graph(std::string output_file) {
		DataOut<2> data_out;
		data_out.attach_dof_handler(*(this->dof_handler));
		data_out.add_data_vector(data, "solution");
		data_out.build_patches();
		std::ofstream output(output_file);
		data_out.write_vtk(output);
	}
};

template<int dim>
class SolutionVector: public Solution {

public:

	DoFHandler<dim> *dof_handler;
	std::vector<std::pair<std::vector<double>, std::vector<double> > >
	coord_and_data_sorted;
	DataOut<dim> data_out; //Used to plot vtk file with deal.ii

	SolutionVector(std::vector<std::pair<std::vector<double>,
			std::vector<double> > > coord_and_data,
			DoFHandler<dim> *dof_handler, DataOut<dim> data_out) {
		this->dof_handler = dof_handler;
		this->data_out = data_out;
		this->coord_and_data_sorted = coord_and_data;
		sort_data_and_point();
	}

	void draw_vtk_graph(std::string output_file) {
		std::ofstream output(output_file);
		data_out.write_vtk(output);
	}

	void sort_data_and_point() {

		auto cmp =
				[](std::pair<std::vector<double>,
						std::vector<double>> const & a,
						std::pair<std::vector<double>,
						std::vector<double>> const & b){
			double epsilon = 0.00001;
			for(int i=dim-1; i>=0; i--){
				if(!Utils::greater_than_or_equals_double((a.first)[i], (b.first)[i], epsilon))
					return true;
				else if(!Utils::less_than_or_equals_double(a.first[i], b.first[i], epsilon))
					return false;
			}
			return false;
				};

		std::sort(coord_and_data_sorted.begin(), coord_and_data_sorted.end(), cmp);
		//Utils::print_vec_of_pair_of_vec(coord_and_data_sorted);
	}
};

#endif
