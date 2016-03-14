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
	typedef std::pair<box,
			std::pair<typename DoFHandler<dim>::active_cell_iterator, float>> cell_coord_pair;
	bgi::rtree<cell_coord_pair, bgi::quadratic<16> > values_at_cells;

	Vector<double> solution_vec;
	Vector<float> uncertainty_per_cell;

	Solution(): potential_v(1), gradient_v(1) {
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
	//These two structures contain data already available in values_at_cells,
	//but it is useful to keep them as such to easily output vtk graph file
	//using deal.ii DataOut class.
	DataOut<dim> fun_drawer;
	DataOut<dim> derivatives_drawer;

	MappingQ1<dim> mapping;
	std::vector<double> potential_v;
	std::vector<Tensor<1, dim> > gradient_v;
	PhysicalValues<dim> extrapol;

	PhysicalValues<dim> extrapolate_values(Point<dim> const &point) {

		FE_Q<dim> fe(1);

		bpoint query_point(point[0], point[1]);
		std::vector<cell_coord_pair> result_s;
		values_at_cells.query(bgi::intersects(query_point),
				std::back_inserter(result_s));
		cell_coord_pair cell_containing_point = result_s[0];

		const Point<dim> p_cell = mapping.transform_real_to_unit_cell(
				cell_containing_point.second.first, point);
		const Quadrature<dim> quadrature(
				GeometryInfo<dim>::project_to_unit_cell(p_cell));
		FEValues<dim> fe_values(mapping, fe, quadrature,
				update_values | update_gradients);
		fe_values.reinit(cell_containing_point.second.first);
		fe_values.get_function_values(solution_vec, potential_v);
		fe_values.get_function_gradients(solution_vec, gradient_v);

		extrapol.potential = potential_v[0];
		extrapol.uncertainty = cell_containing_point.second.second;
		extrapol.electric_field = -gradient_v[0];

		return extrapol;
	}

};

