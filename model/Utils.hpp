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

class Utils {

public:
	static bool equals_double(double a, double b, double epsilon) {
		return abs(a - b) <= epsilon;
	}

	static bool less_than_or_equals_double(double a, double b, double epsilon){
		return a < b || equals_double(a, b, epsilon);
	}

	static bool greater_than_or_equals_double(double a, double b, double epsilon){
		return a > b || equals_double(a, b, epsilon);
	}
};

#endif
