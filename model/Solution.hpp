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

template<unsigned dim>
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
};

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

template<unsigned dim>
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

	ValuesAtPoint<dim> extrapolate_values(std::pair<typename DoFHandler<dim>::active_cell_iterator,
						ValuesAtCell<dim>> const &values_in_cell,
						Point<dim> const &point) {
		double x_left = values_in_cell.first->vertex(0)[0];
		double x_right = values_in_cell.first->vertex(3)[0];
		double y_bottom = values_in_cell.first->vertex(0)[1];
		double y_top = values_in_cell.first->vertex(3)[1];

		double epsilon = 0.00001;

		bool left = false, bottom = false;

		if(Utils::less_than_or_equals_double(abs(x_left-point[0]), abs(x_right-point[0]),
				epsilon))
			left = true;
		if(Utils::less_than_or_equals_double(abs(y_bottom-point[1]), abs(y_top-point[1]),
				epsilon))
			bottom = true;

		ValuesAtPoint<dim> extrapol;
		if(left) {
			if(bottom) {
				extrapol.fun = values_in_cell.second->fun[0] 
						+ (values_in_cell.second->fun[1]-values_in_cell.second->fun[0])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->fun[2]-values_in_cell.second->fun[0])
						* (point[1]-y_bottom)/(y_top-y_bottom);
				extrapol.gradient = values_in_cell.second->gradient[0] 
						+ (values_in_cell.second->gradient[1]-values_in_cell.second->gradient[0])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->gradient[2]-values_in_cell.second->gradient[0])
						* (point[1]-y_bottom)/(y_top-y_bottom);
				extrapol.hessian = values_in_cell.second->hessian[0] 
						+ (values_in_cell.second->hessian[1]-values_in_cell.second->hessian[0])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->hessian[2]-values_in_cell.second->hessian[0])
						* (point[1]-y_bottom)/(y_top-y_bottom);
			}
			else {
				extrapol.fun = values_in_cell.second->fun[2] 
						+ (values_in_cell.second->fun[3]-values_in_cell.second->fun[2])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->fun[0]-values_in_cell.second->fun[2])
						* (point[1]-y_top)/(y_top-y_bottom);
				extrapol.gradient = values_in_cell.second->gradient[2] 
						+ (values_in_cell.second->gradient[3]-values_in_cell.second->gradient[2])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->gradient[0]-values_in_cell.second->gradient[2])
						* (point[1]-y_bottom)/(y_top-y_bottom);
				extrapol.hessian = values_in_cell.second->hessian[2] 
						+ (values_in_cell.second->hessian[3]-values_in_cell.second->hessian[2])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->hessian[0]-values_in_cell.second->hessian[2])
						* (point[1]-y_top)/(y_top-y_bottom);
			}
		}
		else {
			if(bottom) {
				extrapol.fun = values_in_cell.second->fun[1] 
						+ (values_in_cell.second->fun[0]-values_in_cell.second->fun[1])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->fun[3]-values_in_cell.second->fun[1])
						* (point[1]-y_top)/(y_top-y_bottom);
				extrapol.gradient = values_in_cell.second->gradient[1] 
						+ (values_in_cell.second->gradient[0]-values_in_cell.second->gradient[1])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->gradient[3]-values_in_cell.second->gradient[1])
						* (point[1]-y_top)/(y_top-y_bottom);
				extrapol.hessian = values_in_cell.second->hessian[1] 
						+ (values_in_cell.second->hessian[0]-values_in_cell.second->hessian[1])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->hessian[3]-values_in_cell.second->hessian[1])
						* (point[1]-y_top)/(y_top-y_bottom);
			}
			else {
				extrapol.fun = values_in_cell.second->fun[3] 
						+ (values_in_cell.second->fun[2]-values_in_cell.second->fun[3])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->fun[1]-values_in_cell.second->fun[3])
						* (point[1]-y_top)/(y_top-y_bottom);
				extrapol.gradient = values_in_cell.second->gradient[3] 
						+ (values_in_cell.second->gradient[2]-values_in_cell.second->gradient[3])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->gradient[1]-values_in_cell.second->gradient[3])
						* (point[1]-y_top)/(y_top-y_bottom);
				extrapol.hessian = values_in_cell.second->hessian[3] 
						+ (values_in_cell.second->hessian[2]-values_in_cell.second->hessian[3])
						* (point[0]-x_left)/(x_right-x_left)
						+ (values_in_cell.second->hessian[1]-values_in_cell.second->hessian[3])
						* (point[1]-y_top)/(y_top-y_bottom);
			}
		}
		return extrapol;
	}

	/**
	 * Take the coordinates of a point, find the cell in which the point lies
	 * and extrapole the values of fun, gradient, and hessian at this point.
	 */
	ValuesAtPoint<dim> get_values(Point<dim> &point) {
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

		unsigned pos = low - values_at_cells.begin();

		ValuesAtPoint<dim> extrapol;
		extrapol = extrapolate_values(values_at_cells[pos], point);

		return extrapol;
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
