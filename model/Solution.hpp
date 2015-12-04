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

	DoFHandler<dim> *dof_handler;
	Vector<double> data;

	SolutionScalar(Vector<double> data, DoFHandler<dim> *dof_handler){
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
	Vector<double> data;
	DataOut<dim> data_out; //Used to plot vtk file with deal.ii

	SolutionVector(Vector<double> data, DoFHandler<dim> *dof_handler,
			DataOut<dim> data_out){
		this->dof_handler = dof_handler;
		this->data = data;
		this->data_out = data_out;
	}

	void draw_vtk_graph(std::string output_file) {
		std::ofstream output(output_file);
		data_out.write_vtk(output);
	}
};

#endif
