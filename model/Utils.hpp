/*
 * utils.hpp

 *
 *  Created on: 23 nov. 2015
 *      Author: aschils
 */

#pragma once

#include <stdlib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_handler.h>
//#include "Solution.hpp"

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <functional>
#include <math.h>

using namespace dealii;

class Utils {

public:
	static bool equals_double(double a, double b, double epsilon) {
		return fabs(a - b) <= epsilon;
	}

	static bool less_than_or_equals_double(double a, double b, double epsilon) {
		return a < b || equals_double(a, b, epsilon);
	}

	static bool greater_than_or_equals_double(double a, double b,
			double epsilon) {
		return a > b || equals_double(a, b, epsilon);
	}

	static void create_directory_if_not_exists(std::string path) {

		const char *cpath = path.c_str();

		struct stat st;
		if (!(stat(cpath, &st) == 0 && (st.st_mode & S_IFDIR) != 0))
			mkdir(cpath, 0777);
	}

	/**
	 * Sort cells in ascending order i.e. the bottom left point is used
	 * for comparison with the following ordering:
	 *
	 *  coordx = (x_1,...,x_n), coordy = (y_1,...,y_n)
	 *
	 * coordx < coordy if x_n < y_n
	 *  or x_n = y_n and x_(n-1) < y_(n-1)
	 *  or x_n = y_n and x_(n-1) = y_(n-1) and x_(n-2) < y_(n-2)
	 *  and so on.
	 *
	 */

	/*
	 template<unsigned dim, typename T>
	 static void sort_cells_by_coord(
	 std::vector<
	 std::pair<typename DoFHandler<dim>::active_cell_iterator, T> > *v) {
	 auto cmp =
	 [](std::pair<typename DoFHandler<dim>::active_cell_iterator,
	 T> const & a,
	 std::pair<typename DoFHandler<dim>::active_cell_iterator,
	 T> const & b) {

	 Point<dim> bottom_left_point_a = a.first->vertex(0);
	 Point<dim> bottom_left_point_b = b.first->vertex(0);

	 double epsilon = 0.00001;
	 for(int i=dim-1; i>=0; i--) {
	 if(bottom_left_point_a[i] - bottom_left_point_b[i] < -epsilon)
	 return true;
	 else if(bottom_left_point_a[i] - bottom_left_point_b[i] > epsilon)
	 return false;
	 }
	 return false;
	 };
	 std::sort(v->begin(), v->end(), cmp);
	 }
	 */

	template<unsigned dim, typename T>
	static void sort_cells_by_coord(
			std::vector<
					std::pair<typename DoFHandler<dim>::active_cell_iterator, T> > *v) {

		int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;

		auto cmp =
				[&vertices_per_cell](
						std::pair<typename DoFHandler<dim>::active_cell_iterator,
						T> const & a,
						std::pair<typename DoFHandler<dim>::active_cell_iterator,
						T> const & b) {

					Point<dim> bot_left_a = a.first->vertex(0);
					Point<dim> top_right_a = a.first->vertex(vertices_per_cell-1);
					Point<dim> bot_left_b = b.first->vertex(0);
					Point<dim> top_right_b = b.first->vertex(vertices_per_cell-1);

					std::cout << "before bool" << std::endl;

					bool A =  bot_left_a[0] <= bot_left_b[0];
					bool B = top_right_b[0] <= top_right_a[0];
					bool C = bot_left_a[1] <= bot_left_b[1];
					bool D = top_right_b[1] <= top_right_a[1];
					return !(A^B^C^D);
				};
		std::sort(v->begin(), v->end(), cmp);
		std::cout << "AFTER THE BIG SORT" << std::endl;
	}

	template<unsigned dim, typename T>
	static void sort_cells_by_one_coord(
			std::vector<
					std::pair<typename DoFHandler<dim>::active_cell_iterator, T> > *v,
			unsigned vertex_idx, unsigned coord_idx) {

		auto cmp = [&coord_idx, &vertex_idx](
				std::pair<typename DoFHandler<dim>::active_cell_iterator,
				T> const & a,
				std::pair<typename DoFHandler<dim>::active_cell_iterator,
				T> const & b) {
			Point<dim> top_right_a = a.first->vertex(vertex_idx);
			Point<dim> top_right_b = b.first->vertex(vertex_idx);
			return top_right_a[coord_idx] < top_right_b[coord_idx];
		};
		std::sort(v->begin(), v->end(), cmp);
	}

	//Increasing order, smaller if x1 < x2 or x1==x1 and y1<y2...
	template<unsigned dim>
	static void sort_points_by_coord(std::vector<Point<dim> > *v) {
		auto cmp = [](Point<dim> & a, Point<dim> & b) {
			for(unsigned i = 0; i < dim; i++) {
				if(a[i] < b[i])
				return true;
				else if(a[i] > b[i])
				return false;
			}
			return false;
		};
		std::sort(v->begin(), v->end(), cmp);
	}

	/*********** WARNING ***********/
	/*
	 * GAUSSIAN related function not yet tested / fully implemented
	 */

	static double gaussian(double peak, double standard_deviation,
			double center, double x) {
		double exp_arg = -std::pow((x - center), 2)
				/ (2 * std::pow(standard_deviation, 2));
		return peak * std::exp(exp_arg);
	}

	static double gaussian_near_zero_x(double peak, double standard_deviation,
			double center, double max_error) {

		//Init x to value for which f(x) = max/2 (half width at half maximum)
		double x = standard_deviation * sqrt(2 * log(2));
		double y = peak / 2;

		while (!equals_double(y, 0.0, max_error) || x == 0.0) {
			x = x + x / 2.0;
			y = gaussian(peak, standard_deviation, 0.0, x);
		}
		return x + center;
	}

	/*********** END WARNING ***********/

	template<unsigned dim>
	static std::vector<double> point_to_std_vec(Point<dim> p) {
		std::vector<double> v(dim);
		for (unsigned i = 0; i < dim; i++) {
			v[i] = p[i];
		}
		return v;
	}

	template<unsigned dim>
	static void print_point(Point<dim> p) {
		std::cout << "(";
		for (unsigned i = 0; i < dim; i++) {
			std::cout << p[i];
			if (i != dim - 1)
				std::cout << ",";
		}
		std::cout << ")";
	}

	/*
	 * This function write a file to open with gnuplot.
	 * The absciss of each point is the first element of the pair of each element of data
	 * The ordinate are the second element of the pair.
	 */
	template<unsigned dim>
	static void write_gnu_data_file(std::string file,
			std::vector<std::pair<double, double>> data) {

		std::ofstream gnu_graph;
		gnu_graph.open(file,
				std::fstream::in | std::fstream::out | std::fstream::trunc);

		if (gnu_graph.is_open()) {

			gnu_graph << "# X\tY" << std::endl;
			for (unsigned i = 0; i < data.size(); i++) {

				gnu_graph << data[i].first << "\t" << data[i].second
						<< std::endl;
				;
			}
		}
		gnu_graph.close();
	}

	template<unsigned dim>
	static void write_gnu_error_data_file(std::string file,
			std::vector<std::pair<double, std::pair<double, double>>>data) {

				std::ofstream gnu_graph;
				gnu_graph.open(file, std::fstream::in | std::fstream::out | std::fstream::trunc);

				if(gnu_graph.is_open()) {

					gnu_graph << "# X\tY\tZ" << std::endl;
					for(unsigned i = 0; i < data.size(); i++) {

						gnu_graph << data[i].first << "\t"
						<< data[i].second.first << "\t"
						<< data[i].second.second << std::endl;;
					}
				}
				gnu_graph.close();
			}

			static bool is_even(int nbr) {
				return nbr%2 == 0;
			}

			static bool is_odd(int nbr) {
				return !is_even(nbr);
			}

			template <typename T1, typename T2>
			static void write_vector_of_pair(std::string out_file_path,
			std::vector<std::pair<T1, T2> > &vec_of_pair,
			bool header, std::string var_name1, std::string var_name2) {

				std::ofstream out_file;
				out_file.open(out_file_path, std::fstream::trunc);

				if (out_file.is_open()) {

					if(header)
					out_file << var_name1 << "    " << var_name2 << std::endl;

					for(unsigned i=0; i<vec_of_pair.size(); i++) {
						out_file << vec_of_pair[i].first << "    ";
						out_file << vec_of_pair[i].second << std::endl;
					}
				}
				else {
					std::cout << "Error, unable to open file: " << out_file_path << std::endl;
				}

				out_file.close();
			}
		};
