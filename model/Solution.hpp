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

class SolutionData {
public:
	double solution;
	std::vector<double> gradient;
	std::vector<double> seconde_derivatives;

	SolutionData(double solution, std::vector<double> gradient,
			std::vector<double> seconde_derivatives) {
		this->solution = solution;
		this->gradient = gradient;
		this->seconde_derivatives = seconde_derivatives;
	}

	SolutionData(){}
};

template<int dim>
class Solution {

public:

	DoFHandler<dim> *dof_handler;
	std::vector<std::pair<std::vector<double>, SolutionData> > coord_and_data;
	//These two structures contain data already available in coord_and_data,
	//but it is useful to keep them as such to easily output vtk graph file
	//using deal.ii DataOut class.
	Vector<double> solution;
	DataOut<dim> derivatives_data_container;

	Solution() {
	}

	Solution(Vector<double> solution_vec, Derivatives<dim> &derivatives,
			DoFHandler<dim> *dof_handler) {

		this->dof_handler = dof_handler;
		this->solution = solution_vec;

		derivatives_data_container.attach_dof_handler(*dof_handler);
		derivatives_data_container.add_data_vector(solution_vec, derivatives);
		derivatives_data_container.build_patches();

		std::stringstream derivatives_stream;
		//TODO find better way to retrieve data from DataPostProcessor
		derivatives_data_container.write_gnuplot(derivatives_stream);
		std::vector<std::pair<std::vector<double>, std::vector<double> > >
		coord_and_derivatives;
		Utils::parse_gnuplot<dim>(derivatives_stream, dim + dim * dim,
				coord_and_derivatives);

		std::cout << derivatives_stream.str() << std::endl;

		coord_and_data.resize(coord_and_derivatives.size());

		for(unsigned i=0; i< coord_and_derivatives.size(); i++){
			std::pair<std::vector<double>, SolutionData>
			coord_and_data_one_point;
			coord_and_data_one_point.first = coord_and_derivatives[i].first;

			std::vector<double> gradient(dim);
			std::vector<double> second_derivatives(dim*dim);

			for(unsigned j=0; j<dim; j++)
				gradient[j] = coord_and_derivatives[i].second[j];
			for(unsigned j=0; j<(dim*dim); j++)
				second_derivatives[j] = coord_and_derivatives[i].second[dim+j];

			SolutionData data_one_point(solution[i], gradient, second_derivatives);
			coord_and_data_one_point.second = data_one_point;
			coord_and_data[i] = coord_and_data_one_point;
		}
	}

	void draw_vtk_graph_solution(std::string output_file) {
		DataOut<2> data_out;
		data_out.attach_dof_handler(*(this->dof_handler));
		data_out.add_data_vector(solution, "solution");
		data_out.build_patches();
		std::ofstream output(output_file);
		data_out.write_vtk(output);
	}

	void draw_vtk_graph_derivatives(std::string output_file) {
		std::ofstream output(output_file);
		derivatives_data_container.write_vtk(output);
	}

	void sort_by_coord() {
		VectorUtils::sort_by_coord<dim, SolutionData>(coord_and_data);
	}
};

#endif
