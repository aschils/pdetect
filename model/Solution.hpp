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
#include "VectorUtils.hpp"

using namespace dealii;

class Solution {

public:

	virtual void draw_vtk_graph(std::string output_file) = 0;

	virtual void sort_by_coord() = 0;

	virtual ~Solution() {

	}
};

template<int dim>
class SolutionScalar: public Solution {

public:

	std::vector<std::vector<double> > coord;

	DoFHandler<dim> *dof_handler;
	Vector<double> data;

	SolutionScalar() {
	}

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

	void sort_by_coord(){
		//TODO
	}
};

template<int dim>
class SolutionVector: public Solution {

public:

	DoFHandler<dim> *dof_handler;
	std::vector<std::pair<std::vector<double>, std::vector<double> > > coord_and_data;
	DataOut<dim> data_out; //Used to plot vtk file with deal.ii

	SolutionVector() {
	}

	SolutionVector(
			std::vector<std::pair<std::vector<double>, std::vector<double> > > coord_and_data,
			DoFHandler<dim> *dof_handler, DataOut<dim> data_out) {
		this->dof_handler = dof_handler;
		this->data_out = data_out;
		this->coord_and_data = coord_and_data;
	}

	void draw_vtk_graph(std::string output_file) {
		std::ofstream output(output_file);
		data_out.write_vtk(output);
	}

	void sort_by_coord(){
		VectorUtils::sort_by_coord<dim, std::vector<double>>(coord_and_data);
	}

};

#endif
