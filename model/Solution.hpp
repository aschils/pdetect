/*
 * Solution.hpp
 *
 *  Created on: 29 nov. 2015
 *      Author: aschils
 */

#pragma once

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
class PhysicalValues {
public:
	double potential;
	Tensor<1, dim> electric_field;

	PhysicalValues(double potential, Tensor<1, dim> electric_field) {
		this->potential = potential;
		this->electric_field = electric_field;
	}

	PhysicalValues() {}
};

template<unsigned dim>
class ValuesAtCell {
public:
	std::vector<double> fun;
	std::vector<Tensor<1, dim> > gradient;
	std::vector<Tensor<2, dim> > hessian;

	ValuesAtCell(std::vector<double> fun, std::vector<Tensor<1, dim> > gradient,
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

	std::vector<
			std::pair<typename DoFHandler<dim>::active_cell_iterator,
					ValuesAtCell<dim> > > values_at_cells;

	Solution() {
	}

	void set_fun_drawer(DataOut<dim> fun_drawer) {
		this->fun_drawer = fun_drawer;
	}

	void set_derivatives_drawer(DataOut<dim> derivatives_drawer) {
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

	void sort_cells_by_coord() {
		Utils::sort_cells_by_coord<dim, ValuesAtCell<dim> >(&values_at_cells);
	}

	unsigned get_cell(Point<dim> const &point) {
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
									point[dim-1], epsilon))//Before the cell
					return true;
					else
					return false;
				};

		typename std::vector<
				std::pair<typename DoFHandler<dim>::active_cell_iterator,
						ValuesAtCell<dim>>>::iterator low;

		//iterator that will find the cell in which lies point, using the comparison 
		//algorithm cmp
		low = std::lower_bound(values_at_cells.begin(), values_at_cells.end(),
				point, cmp);

		unsigned pos = low - values_at_cells.begin();
		if (pos == values_at_cells.size())
			pos--;

		return pos;
	}

	std::pair<unsigned, unsigned> get_closest_point(Point<dim> const &point) {

		unsigned pos = get_cell(point);
		std::pair<unsigned, unsigned> closest_point;

		double x_left = values_at_cells[pos].first->vertex(0)[0];
		double x_right = values_at_cells[pos].first->vertex(3)[0];
		double y_bottom = values_at_cells[pos].first->vertex(0)[1];
		double y_top = values_at_cells[pos].first->vertex(3)[1];

		double epsilon = 0.00001;

		bool left = false, bottom = false;

		if (Utils::less_than_or_equals_double(fabs(x_left - point[0]),
				fabs(x_right - point[0]), epsilon))
			left = true;
		if (Utils::less_than_or_equals_double(fabs(y_bottom - point[1]),
				fabs(y_top - point[1]), epsilon))
			bottom = true;

		unsigned closest;

		if (left) {
			if (bottom)
				closest = 0;
			else
				closest = 2;
		} else {
			if (bottom)
				closest = 1;
			else
				closest = 3;
		}

		closest_point.first = pos;
		closest_point.second = closest;

		return closest_point;
	}

	PhysicalValues<dim> extrapolate_values(Point<dim> const &point) {

		std::pair<unsigned, unsigned> closest_point = get_closest_point(point);
		unsigned pos = closest_point.first;
		unsigned closest = closest_point.second;

		PhysicalValues<dim> extrapol;

		double x_left = values_at_cells[pos].first->vertex(0)[0];
		double x_right = values_at_cells[pos].first->vertex(3)[0];
		double y_bottom = values_at_cells[pos].first->vertex(0)[1];
		double y_top = values_at_cells[pos].first->vertex(3)[1];

		double delta_x;
		double delta_y;
		if (closest == 0 || closest == 2)
			delta_x = point[0] - x_left;
		else
			delta_x = point[0] - x_right;
		if (closest == 0 || closest == 1)
			delta_y = point[1] - y_bottom;
		else
			delta_y = point[1] - y_top;

		//We use the Taylor expansion to calculate the values of the function and the gradient
		extrapol.potential = values_at_cells[pos].second.fun[closest]
				+ values_at_cells[pos].second.gradient[closest][0] * delta_x
				+ values_at_cells[pos].second.gradient[closest][1] * delta_y
				+ 0.5*(values_at_cells[pos].second.hessian[closest][0][0] * delta_x * delta_x
				+ values_at_cells[pos].second.hessian[closest][1][1] * delta_y * delta_y)
				+ values_at_cells[pos].second.hessian[closest][0][1] * delta_y * delta_x;

		for(unsigned i = 0; i < dim; i++) {
			extrapol.electric_field[i] = -(values_at_cells[pos].second.gradient[closest][i]
					+ values_at_cells[pos].second.hessian[closest][i][0] * delta_x
					+ values_at_cells[pos].second.hessian[closest][i][1] * delta_y);
		}

		return extrapol;
	}

	/**
	 * Take the coordinates of a point, find the cell in which the point lies
	 * and extrapole the values of fun, gradient, and hessian at this point.
	 */
	PhysicalValues<dim> get_values(Point<dim> const &point) {

		PhysicalValues<dim> extrapol;
		extrapol = extrapolate_values(point);

		return extrapol;
	}

	void print() {

		for (unsigned i = 0; i < values_at_cells.size(); i++) {
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

	DataOut<dim> fun_drawer;
	DataOut<dim> derivatives_drawer;
};

