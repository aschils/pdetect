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

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <functional>

using namespace dealii;

class Utils {

public:
	static bool equals_double(double a, double b, double epsilon) {
		return abs(a - b) <= epsilon;
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

	template<unsigned dim>
	static void parse_gnuplot(std::stringstream &gnuplot_stream,
			unsigned nbr_of_components,
			std::vector<std::pair<std::vector<double>, std::vector<double> > > &coord_and_data) {

		std::string line;
		unsigned numbers_per_line = dim + nbr_of_components;

		while (getline(gnuplot_stream, line)) {

			if (line.length() == 0 || line[0] == '#' || line[0] == '\n')
				continue;

			std::vector<double> numbers = parse_gnuplot_line(line);
			if (numbers.size() != numbers_per_line)
				continue;

			std::vector<double> point_coord(dim);
			for (unsigned i = 0; i < dim; i++)
				point_coord[i] = numbers[i];

			std::vector<double> point_data(nbr_of_components);
			for (unsigned i = dim; i < numbers_per_line; i++)
				point_data[i - dim] = numbers[i];

			std::pair<std::vector<double>, std::vector<double> > coord_and_data_pair(
					point_coord, point_data);
			coord_and_data.push_back(coord_and_data_pair);
		}
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

	template<unsigned dim, typename T>
	static void sort_cells_by_coord(std::vector<
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

private:

	static std::vector<double> parse_gnuplot_line(std::string line) {

		std::vector<double> numbers;
		std::stringstream stream(line);
		std::string token;

		while (getline(stream, token, ' ')) {
			double number = std::stod(token);
			numbers.push_back(number);
		}
		return numbers;
	}
};
