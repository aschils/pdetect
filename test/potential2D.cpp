/*
 * potential2D.cpp
 *
 *  Created on: 8 nov. 2015
 *      Author: aschils
 */

#include "../model/LaplaceSolver.hpp"
#include "../model/Rect2DBoundaryValues.hpp"
#include "../model/ZeroRightHandSide.hpp"

void test_2D_potential(){

	double rect_length_fe = 4;

	ZeroRightHandSide<2> rhs;
	double strip_potential = 1;
	unsigned strip_width = 1;
	unsigned pitch = 1;

	for(unsigned nbr_of_strips=0; nbr_of_strips<=5; nbr_of_strips++){
		Rect2DBoundaryValues<2> bv(nbr_of_strips, strip_width, pitch,
				rect_length_fe, strip_potential);
		LaplaceSolver<2> rect_potential_solver(rect_length_fe, &rhs, &bv, 
				std::to_string(nbr_of_strips)+".vtk");
		rect_potential_solver.run();
	}
}

