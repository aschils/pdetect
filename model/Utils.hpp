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

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unordered_map>
#include <functional>

template<>
struct std::hash<std::pair<double, double> > {
public:
	size_t operator()(std::pair<double, double> x) const throw () {
		return std::hash<double>()(x.first) ^ std::hash<double>()(x.second);
	}
};

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

	/**
	 * Parse gnuplot file with 4 decimal numbers per line x,y,v1,v2.
	 * File format:
	 * #
	 * # Gnu plot comments are ignored
	 * .
	 * .
	 * .
	 * #
	 *
	 * x y v1 v2
	 * .
	 * .
	 * .
	 *
	 * The returned hashmap
	 *
	 */
	static std::unordered_map<std::pair<double, double>,
			std::pair<double, double> > parse_gnuplot_2D(
			std::stringstream &gnuplot_stream) {

		std::unordered_map<std::pair<double, double>, std::pair<double, double> > data;
		std::string line;

		while (getline(gnuplot_stream, line)) {

			if (line.length() == 0 || line[0] == '#' || line[0] == '\n')
				continue;

			std::vector<double> numbers = parse_gnuplot_line(line);

			if (numbers.size() != 4)
				continue;

			std::pair<double, double> coord(numbers[0], numbers[1]);
			std::pair<double, double> vector_values(numbers[2], numbers[3]);

			data[coord] = vector_values;
		}

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
