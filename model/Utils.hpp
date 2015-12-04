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

	template <unsigned dim>
	static Vector<double> parse_gnuplot(std::stringstream &gnuplot_stream) {


		std::vector<double> vec;
		std::string line;
		unsigned numbers_per_line = 2*dim;

		unsigned count = 0;

		while (getline(gnuplot_stream, line)) {

			if (line.length() == 0 || line[0] == '#' || line[0] == '\n')
				continue;

			std::vector<double> numbers = parse_gnuplot_line(line);

			if (numbers.size() != numbers_per_line)
				continue;

			for(unsigned i=dim; i<numbers_per_line; i++)
				vec.push_back(numbers[i]);
		}

		Vector<double> data(vec.size());
		for(unsigned i=0; i<vec.size(); i++)
			data[i] = vec[i];
		return data;
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
