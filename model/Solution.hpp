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
	//typedef std::pair<box, ValuesAtCell<dim>> cell_coord_pair;
	typedef std::pair<box, std::pair<
				typename DoFHandler<dim>::active_cell_iterator, float>>
		cell_coord_pair;
	bgi::rtree<cell_coord_pair, bgi::quadratic<16> > values_at_cells;

	DoFHandler<dim> *dof_handler;
	Vector<double> solution_vec;
	Vector<float> uncertainty_per_cell;

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
		std::vector<cell_coord_pair> result_s;
		values_at_cells.query(bgi::intersects(query_point),
				std::back_inserter(result_s));

		cell_coord_pair cell_containing_point = result_s[0];
		//value values_at_cell = get_closest_cell(result_s, point);
		//bpoint bottom_left = cell_containing_point.first.min_corner();


		/*double x_left = bottom_left.get<0>();
		double y_bottom = bottom_left.get<1>();
		 double delta_x = point[0] - x_left;
		 double delta_y = point[1] - y_bottom;*/

		/*
		 double x = point[0];
		 double y = point[1];
		 bpoint top_right = values_at_cell.first.max_corner();
		 double x_right = top_right.get<0>();
		 double y_top = top_right.get<1>();

		 if (x < x_left || x > x_right || y < y_bottom || y > y_top) {

		 std::cout << "Point: (" << point[0] << "," << point[1] << ")"
		 << std::endl;

		 std::cout << "bot_left: (" << x_left << "," << y_bottom
		 << ") bot_right: (" << x_right << "," << y_bottom
		 << ") top_left: (" << x_left << "," << y_top
		 << ") top_right: (" << x_right << "," << y_top << ")"
		 << std::endl;
		 }*/

		PhysicalValues<dim> extrapol;

		//We use the Taylor expansion to calculate the values of the function and the gradient
//		extrapol.potential = values_at_cell.second.fun[0].first
//				+ values_at_cell.second.gradient[0][0] * delta_x
//				+ values_at_cell.second.gradient[0][1] * delta_y
//				+ 0.5
//						* (values_at_cell.second.hessian[0][0][0] * delta_x
//								* delta_x
//								+ values_at_cell.second.hessian[0][1][1]
//										* delta_y * delta_y)
//				+ values_at_cell.second.hessian[0][0][1] * delta_y * delta_x;
//
//		extrapol.uncertainty = values_at_cell.second.fun[0].second;
//
//		for (unsigned i = 0; i < dim; i++) {
//			extrapol.electric_field[i] = -(values_at_cell.second.gradient[0][i]
//					+ values_at_cell.second.hessian[0][i][0] * delta_x
//					+ values_at_cell.second.hessian[0][i][1] * delta_y);
//		}

		FE_Q<dim> fe(1);
		MappingQ1<dim> mapping;
		//const std::pair<typename DoFHandler<dim, dim>::active_cell_iterator,
		//		Point<dim> > cell_point = cell_containing_point.second;

		const Point<dim> p_cell = mapping.transform_real_to_unit_cell(
				cell_containing_point.second.first, point);

		const Quadrature<dim> quadrature(
				GeometryInfo<dim>::project_to_unit_cell(p_cell));
		FEValues<dim> fe_values(mapping, fe, quadrature,
				update_values | update_gradients);
		fe_values.reinit(cell_containing_point.second.first);

		std::vector<double> potential_v(1);
		std::vector<Tensor<1, dim> > gradient_v(1);
		fe_values.get_function_values(solution_vec, potential_v);
		fe_values.get_function_gradients(solution_vec, gradient_v);
		extrapol.potential = potential_v[0];
		extrapol.uncertainty = cell_containing_point.second.second;
		extrapol.electric_field = -gradient_v[0];

		/*std::vector<std::pair<double, double>> fun_at_cell;
				for (int i = 0; i < vertices_per_cell; i++) {
					fun_at_cell.push_back(
							std::pair<double, double>(fun_at_pts_of_one_cell[i],
									uncertainty_per_cell[pos]));
				}*/

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

	/*
	 value get_closest_cell(std::vector<value> &result_s,
	 Point<dim> const &point) {

	 double min_dist = std::numeric_limits<double>::max();
	 value closest;

	 for (unsigned i = 0; i < result_s.size(); i++) {

	 bpoint bottom_left = result_s[i].first.min_corner();
	 double x = bottom_left.get<0>();
	 double y = bottom_left.get<1>();
	 double dist = std::pow(x - point[0], 2) + std::pow(y - point[1], 2);

	 if (dist <= min_dist) {
	 min_dist = dist;
	 closest = result_s[i];
	 }
	 }

	 return closest;

	 }*/

};

