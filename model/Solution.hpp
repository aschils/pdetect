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

template<unsigned dim>
class Solution {

public:

	std::vector<
			std::pair<typename DoFHandler<dim>::active_cell_iterator,
					ValuesAtCell<dim> > > values_at_cells;
	DoFHandler<dim> *dof_handler;
	Vector<double> solution_vec;

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
		//Utils::sort_cells_by_coord<dim, ValuesAtCell<dim> >(&values_at_cells);
		//std::cout << "AFTER THE BIG SORT" << std::endl;
		/*unsigned x_idx = 0;
		int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
		unsigned top_right_vertex_idx = vertices_per_cell-1;
		Utils::sort_cells_by_one_coord<dim, ValuesAtCell<dim> >(
				&values_at_cells, top_right_vertex_idx, x_idx);*/
	}

	/*
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

	 typename std::vector<std::pair<typename DoFHandler<dim>::active_cell_iterator,
	 ValuesAtCell<dim>>>::iterator low;

	 //iterator that will find the cell in which lies point, using the comparison
	 //algorithm cmp
	 low = std::lower_bound(values_at_cells.begin(), values_at_cells.end(),
	 point, cmp);

	 unsigned pos = low - values_at_cells.begin();
	 if (pos == values_at_cells.size())
	 pos--;

	 return pos;
	 }*/

//	void get_cells_in_x_range(Point<dim> const &point,
//			std::vector<
//					std::pair<typename DoFHandler<dim>::active_cell_iterator,
//							ValuesAtCell<dim>>>&cells_in_x_range) {
//
//				int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
//
//				auto cmp_x_coord =
//				[&vertices_per_cell](std::pair<typename DoFHandler<dim>::active_cell_iterator,
//				ValuesAtCell<dim>> const &values_in_cell,
//				Point<dim> const &point) {
//					//Point<dim> bottom_left = values_in_cell.first->vertex(0);
//					Point<dim> top_right = values_in_cell.first->vertex(vertices_per_cell-1);
//					//return bottom_left[0] < point[0];
//					return top_right[0] < point[0];
//				};
//
//				typename std::vector<
//				std::pair<typename DoFHandler<dim>::active_cell_iterator,
//				ValuesAtCell<dim>>>::iterator low =
//				std::lower_bound(values_at_cells.begin(), values_at_cells.end(),
//				point, cmp_x_coord);
//
//				std::vector<std::pair<typename DoFHandler<dim>::active_cell_iterator,
//							ValuesAtCell<dim>>> cells_top_right_x_ge(
//									low, values_at_cells.end());
//
///*
//
//
//				typename DoFHandler<dim>::active_cell_iterator cell =
//						values_at_cells[pos].first;
//
//				Point<dim> top_right = cell->vertex(vertices_per_cell-1);
//				double first_top_right_x = top_right[0];
//				double epsilon = 0.0000001;
//
//				while(Utils::equals_double(first_top_right_x, top_right[0], epsilon)) {
//					cells_in_x_range.push_back(values_at_cells[pos]);
//					pos++;
//
//					if(pos == values_at_cells.size())
//						break;
//
//					cell = values_at_cells[pos].first;
//					top_right = cell->vertex(vertices_per_cell-1);
//				}*/
//			}
//
//			std::pair<typename DoFHandler<dim>::active_cell_iterator,
//			ValuesAtCell<dim>> get_cell(Point<dim> const &point) {
//
//				//std::cout << "in get_cell" << std::endl;
//				std::vector<std::pair<typename DoFHandler<dim>::active_cell_iterator,
//				ValuesAtCell<dim>>> cells_in_x_range;
//
//				get_cells_in_x_range(point, cells_in_x_range);
//
//				//std::cout << "after get_cells_in_x_range" << std::endl;
//				//std::cout << "nbr_of_cells_in_x_range: " << cells_in_x_range.size() << std::endl;
//
//				unsigned y_idx = 1;
//				//Utils::sort_cells_by_one_coord<dim, ValuesAtCell<dim>>(
//				//&cells_in_x_range, y_idx);
//
//				/*for(unsigned i=0; i<cells_in_x_range.size(); i++){
//				 double y = cells_in_x_range[i].first->vertex(0)[1];
//				 std::cout << y << " ";
//				 }*/
//
//				//std::cout << "after sort_cells_by_one_coord" << std::endl;
//				auto cmp_y_coord =
//				[](std::pair<typename DoFHandler<dim>::active_cell_iterator,
//				ValuesAtCell<dim>> const &values_in_cell,
//				Point<dim> const &point) {
//					int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
//					Point<dim> top_right = values_in_cell.first->vertex(vertices_per_cell-1);
//					//std::cout << "in comp bazarre, " << (top_right[1] < point[1]) << std::endl;
//					return top_right[1] < point[1];
//				};
//
//				typename std::vector<
//				std::pair<typename DoFHandler<dim>::active_cell_iterator,
//				ValuesAtCell<dim>>>::iterator low;
//
//				//iterator that will find the cell in which lies point, using the comparison
//				//algorithm cmp
//				low = std::lower_bound(cells_in_x_range.begin(), cells_in_x_range.end(),
//				point, cmp_y_coord);
//
//				//std::cout << "after second lower_bound" << std::endl;
//
//				std::pair<typename DoFHandler<dim>::active_cell_iterator,
//				ValuesAtCell<dim>> found_cell;
//
//				if(low != cells_in_x_range.end()) {
//					found_cell = *low;
//				}
//				else {
//					found_cell = *(--low);
//				}
//
//				//std::cout << found_cell.first->vertex(0) << std::endl;
//
//				return found_cell;
//			}


			std::pair<typename DoFHandler<dim>::active_cell_iterator,
				ValuesAtCell<dim>> get_cell(Point<dim> const &point) {

				for(unsigned i=0; i<values_at_cells.size(); i++){

					int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
					Point<dim> top_right = values_at_cells[i].first->vertex(vertices_per_cell-1);
					Point<dim> bottom_left = values_at_cells[i].first->vertex(0);

					if(bottom_left[0] <= point[0] && top_right[0] >= point[0]
						&& bottom_left[1] <= point[1] && top_right[1] >= point[1]){
						return values_at_cells[i];
					}

				}

			}

			PhysicalValues<dim> extrapolate_values(Point<dim> const &point) {

				//std::pair<unsigned, unsigned> closest_point = get_closest_point(point);
				//unsigned pos = closest_point.first;
				//unsigned closest = closest_point.second;

				//unsigned pos = get_cell(point);

				std::pair<typename DoFHandler<dim>::active_cell_iterator,
				ValuesAtCell<dim>> values_at_cell = get_cell(point);

				//std::cout << "after get_cell" << std::endl;

				PhysicalValues<dim> extrapol;

				double x_left = values_at_cell.first->vertex(0)[0];
				double x_right = values_at_cell.first->vertex(3)[0];
				double y_bottom = values_at_cell.first->vertex(0)[1];
				double y_top = values_at_cell.first->vertex(3)[1];

				double x = point[0];
				double y = point[1];

				//std::cout << "before if error" << std::endl;

//				if (x < x_left || x > x_right || y < y_bottom || y > y_top) {
//
//					std::cout << "Point: (" << point[0] << "," << point[1] << ")"
//					<< std::endl;
//
//					std::cout << "bot_left: (" << x_left << "," << y_bottom
//					<< ") bot_right: (" << x_right << "," << y_bottom
//					<< ") top_left: (" << x_left << "," << y_top
//					<< ") top_right: (" << x_right << "," << y_top << ")"
//					<< std::endl;
//				}
				double delta_x;
				double delta_y;
				//if (closest == 0 || closest == 2)
				delta_x = point[0] - x_left;
				/*else
				 delta_x = point[0] - x_right;*/
				//if (closest == 0 || closest == 1)
				delta_y = point[1] - y_bottom;
				//else
				//delta_y = point[1] - y_top;

				//std::cout << "before first taylor" << std::endl;

				//We use the Taylor expansion to calculate the values of the function and the gradient
				extrapol.potential = values_at_cell.second.fun[0].first
				+ values_at_cell.second.gradient[0][0] * delta_x
				+ values_at_cell.second.gradient[0][1] * delta_y
				+ 0.5
				* (values_at_cell.second.hessian[0][0][0]
				* delta_x * delta_x
				+ values_at_cell.second.hessian[0][1][1]
				* delta_y * delta_y)
				+ values_at_cell.second.hessian[0][0][1] * delta_y
				* delta_x;

				extrapol.uncertainty = values_at_cell.second.fun[0].second;

				for (unsigned i = 0; i < dim; i++) {
					extrapol.electric_field[i] =
					-(values_at_cell.second.gradient[0][i]
					+ values_at_cell.second.hessian[0][i][0]
					* delta_x
					+ values_at_cell.second.hessian[0][i][1]
					* delta_y);
				}
				return extrapol;
			}

			/**
			 * Take the coordinates of a point, find the cell in which the point lies
			 * and extrapole the values of fun, gradient, and hessian at this point.
			 */
			PhysicalValues<dim> get_values(Point<dim> const &point) {

				/*
				 * 1. Create a quadrature that only contains your vertices
				 (QIterated<dim> of (QTrapez<1>(),1) should work)
				 2. On each cell (or face), init an FEValues (or FEFaceValues) as if
				 you were assembling.
				 3. fe_values.get_function_values(solution_vector, values) will give
				 you the values at the points specified in 1)
				 */

				PhysicalValues<dim> extrapol;
				//extrapol.potential = VectorTools::point_value(*dof_handler,
				//		solution_vec, point);
				//extrapol.electric_field = -VectorTools::point_gradient(*dof_handler,
				//		solution_vec, point);

				/*
				 QIterated<dim> quad(QTrapez<1>(),1);
				 FE_Q<dim> fe(1);
				 FEValues<dim> fe_values(fe, quad,
				 update_values | update_gradients | update_quadrature_points
				 | update_JxW_values | update_second_derivatives);*/

				extrapol = extrapolate_values(point);
				/*FE_Q<dim> fe(1);
				 MappingQ1<dim> mapping;
				 const std::pair<typename DoFHandler<dim, dim>::active_cell_iterator,
				 Point<dim> > cell_point =
				 GridTools::find_active_cell_around_point(mapping, *dof_handler,
				 point);

				 const Quadrature<dim> quadrature(
				 GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
				 FEValues<dim> fe_values(mapping, fe, quadrature,
				 update_values | update_gradients);
				 fe_values.reinit(cell_point.first);

				 std::vector<double> potential_v(1);
				 std::vector<Tensor<1, dim> > gradient_v(1);
				 fe_values.get_function_values(solution_vec, potential_v);
				 fe_values.get_function_gradients(solution_vec, gradient_v);
				 extrapol.potential = potential_v[0];
				 extrapol.electric_field = -gradient_v[0];*/

				return extrapol;
			}

			/*void print() {

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
			 }*/

		private:
			//These two structures contain data already available in coord_and_data,
			//but it is useful to keep them as such to easily output vtk graph file
			//using deal.ii DataOut class.
			DataOut<dim> fun_drawer;
			DataOut<dim> derivatives_drawer;
		};

