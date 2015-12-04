/*
 * utils.hpp

 *
 *  Created on: 23 nov. 2015
 *      Author: aschils
 */

#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <stdlib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unordered_map>
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

	static void print_vec_components(std::vector<double> vec) {
		for (unsigned j = 0; j < vec.size(); j++) {
			std::cout << vec[j];

			if (j != vec.size() - 1)
				std::cout << ",";
		}
	}

	static void print_vec_of_pair_of_vec(
			std::vector<std::pair<std::vector<double>, std::vector<double> > > vec) {
		for (unsigned i = 0; i < vec.size(); i++) {

			std::cout << "[(";
			print_vec_components(vec[i].first);
			std::cout << ") (";
			print_vec_components(vec[i].second);
			std::cout << ")] ";
		}
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

#endif
