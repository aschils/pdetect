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
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_tools.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

#include <boost/foreach.hpp>

#include "TensorUtils.hpp"
#include "Utils.hpp"

using namespace dealii;

template<unsigned dim>
class PhysicalValues {
public:
	double potential;
	double uncertainty;
	Tensor<1, dim> electric_field;

	PhysicalValues(double potential, double uncertainty,
			Tensor<1, dim> electric_field) {
		this->potential = potential;
		this->uncertainty = uncertainty;
		this->electric_field = electric_field;
	}

	PhysicalValues() {
	}
};

template<unsigned dim>
class ValuesAtCell {
public:
	std::vector<std::pair<double, double>> fun;
	std::vector<Tensor<1, dim> > gradient;
	std::vector<Tensor<2, dim> > hessian;

	ValuesAtCell(std::vector<std::pair<double, double>> fun,
			std::vector<Tensor<1, dim> > gradient,
			std::vector<Tensor<2, dim> > hessian) {
		this->fun = fun;
		this->gradient = gradient;
		this->hessian = hessian;
	}

	ValuesAtCell() {
	}
};

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<unsigned dim>
class Solution {

public:

	typedef bg::model::point<double, 2, bg::cs::cartesian> bpoint;
	typedef bg::model::box<bpoint> box;
	typedef std::pair<box, ValuesAtCell<dim>> value;

	bgi::rtree<value, bgi::quadratic<16> > values_at_cells;

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

	PhysicalValues<dim> extrapolate_values(Point<dim> const &point) {

		bpoint query_point(point[0], point[1]);
		std::vector<value> result_s;
		values_at_cells.query(bgi::intersects(query_point), std::back_inserter(result_s));
		value values_at_cell = result_s[0];
		bpoint bottom_left = values_at_cell.first.min_corner();

		double x_left = bottom_left.get<0>();
		double y_bottom = bottom_left.get<1>();
		double delta_x = point[0] - x_left;
		double delta_y = point[1] - y_bottom;

		PhysicalValues<dim> extrapol;

		//We use the Taylor expansion to calculate the values of the function and the gradient
		extrapol.potential = values_at_cell.second.fun[0].first
				+ values_at_cell.second.gradient[0][0] * delta_x
				+ values_at_cell.second.gradient[0][1] * delta_y
				+ 0.5
						* (values_at_cell.second.hessian[0][0][0] * delta_x
								* delta_x
								+ values_at_cell.second.hessian[0][1][1]
										* delta_y * delta_y)
				+ values_at_cell.second.hessian[0][0][1] * delta_y * delta_x;

		extrapol.uncertainty = values_at_cell.second.fun[0].second;

		for (unsigned i = 0; i < dim; i++) {
			extrapol.electric_field[i] = -(values_at_cell.second.gradient[0][i]
					+ values_at_cell.second.hessian[0][i][0] * delta_x
					+ values_at_cell.second.hessian[0][i][1] * delta_y);
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

private:
	//These two structures contain data already available in coord_and_data,
	//but it is useful to keep them as such to easily output vtk graph file
	//using deal.ii DataOut class.
	DataOut<dim> fun_drawer;
	DataOut<dim> derivatives_drawer;
};

