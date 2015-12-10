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
#include <algorithm>

#include "TensorUtils.hpp"
#include "Utils.hpp"

using namespace dealii;

/*template<int dim>
class ValuesAtPoint {
public:
	double fun;
	Tensor<1, dim> gradient;
	Tensor<2, dim> hessian;

	ValuesAtPoint(double fun, 
			Tensor<1, dim> gradient,
			Tensor<2, dim> hessian) {
		this->fun = fun;
		this->gradient = gradient;
		this->hessian = hessian;
	}

	ValuesAtPoint() {
	}
};*/

template<unsigned dim>
class ValuesAtCell {
public:
	std::vector<double> fun;
	std::vector<Tensor<1, dim> >  gradient;
	std::vector<Tensor<2, dim> > hessian;

	ValuesAtCell(std::vector<double> fun,
			std::vector<Tensor<1, dim> > gradient,
			std::vector<Tensor<2, dim> > hessian) {
		this->fun = fun;
		this->gradient = gradient;
		this->hessian = hessian;
	}

	ValuesAtCell() {
	}
};

template<int dim>
class Solution {

public:

	std::vector<std::pair<typename DoFHandler<dim>::active_cell_iterator,
	ValuesAtCell<dim> > > values_at_cells;

	Solution() {
	}

	void set_fun_drawer(DataOut<dim> fun_drawer){
		//this->fun_drawer = fun_drawer;
	}

	void set_derivatives_drawer(DataOut<dim> derivatives_drawer){
		//this->derivatives_drawer = derivatives_drawer;
	}

	void draw_vtk_graph_fun(std::string output_file) {
		std::ofstream output(output_file);
		//fun_drawer.write_vtk(output);
	}

	void draw_vtk_graph_derivatives(std::string output_file) {
		std::ofstream output(output_file);
		//derivatives_drawer.write_vtk(output);
	}

	void sort_cells_by_coord() {
		Utils::sort_cells_by_coord<dim, ValuesAtCell<dim> >(&values_at_cells);
	}

	/**
	 * Take the coordinates of a point, find the cell in which the point lies
	 * and extrapole the values of fun, gradient, and hessian at this point.
	 */
	void get_values(Point<dim> &point) {
		auto cmp =
				[](std::pair<typename DoFHandler<dim>::active_cell_iterator,
						ValuesAtCell<dim>> const &values_in_cell,
						Point<dim> const &point) {

					int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
					Point<dim> bottom_left = values_in_cell.first->vertex(0);
					Point<dim> top_right = values_in_cell.first->vertex(vertices_per_cell-1);

					double epsilon = 0.0000001;

					if(Utils::less_than_or_equals_double(top_right[dim-1], point[dim-1],
							epsilon)) //Under the cell
						return true;
					else if(Utils::less_than_or_equals_double(top_right[dim-2], point[dim-2],
							epsilon) && Utils::less_than_or_equals_double(bottom_left[dim-1], 
							point[dim-1], epsilon)) //Before the cell
						return true;
					else
						return false;
				};

		typename std::vector<std::pair<typename DoFHandler<dim>::active_cell_iterator,
						ValuesAtCell<dim>>>::iterator low;
		low = std::lower_bound(values_at_cells.begin(), values_at_cells.end(), point, cmp);

		int pos = low - values_at_cells.begin();

		/*std::cout << "x in [" << values_at_cells[pos].first->vertex(0)[0] << "," 
				<< values_at_cells[pos].first->vertex(3)[0] << "]" << std::endl 
				<< "y in [" << values_at_cells[pos].first->vertex(0)[1] << ","
				<< values_at_cells[pos].first->vertex(3)[1] << "]" << std::endl;*/
	}

	void print(){

		for(unsigned i=0; i<values_at_cells.size(); i++){
			std::cout << "[coord: (";
			TensorUtils::print_vec_components(values_at_cells[i].first);
			std::cout << "), function: " << values_at_cells[i].second.fun;
			std::cout << " gradient: (";
			TensorUtils::print_tensor_components<dim>(
					values_at_cells[i].second.gradient);
			std::cout << ") hessian: (";
			TensorUtils::print_tensor_components<dim>(
								values_at_cells[i].second.hessian);
			std::cout << ")] ";
		}
		std::cout << std::endl;
	}

private:
	//These two structures contain data already available in coord_and_data,
	//but it is useful to keep them as such to easily output vtk graph file
	//using deal.ii DataOut class.
	//DataOut<dim> fun_drawer;
	//DataOut<dim> derivatives_drawer;
};

#endif
