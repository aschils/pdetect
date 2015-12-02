/*
 * GraphDrawer.hpp
 *
 *  Created on: 29 nov. 2015
 *      Author: aschils
 */

#ifndef __GRAPH_DRAWER_HPP__
#define __GRAPH_DRAWER_HPP__

#include <fstream>
#include <iostream>
#include <deal.II/numerics/data_out.h>
#include "../model/Solution.hpp"

template<int dim>
class GraphDrawer {

public:

	static void draw_vtk_graph(Solution<dim> sol, std::string output_file) {
		DataOut<dim> data_out;
		data_out.attach_dof_handler(*sol.dof_handler);
		data_out.add_data_vector(sol.data, "solution");
		data_out.build_patches();
		std::ofstream output(output_file);
		data_out.write_vtk(output);
	}

	static void draw_vtk_graph(DataOut<dim> data_out, std::string output_file) {
		std::ofstream output(output_file);
		data_out.write_vtk(output);
	}

};

#endif
