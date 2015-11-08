/*
 * potential2D.cpp
 *
 *  Created on: 8 nov. 2015
 *      Author: aschils
 */

#include "../model/LaplaceSolver.hpp"
#include "../model/Rect2DBoundaryValues.hpp"
#include "../model/ZeroRightHandSide.hpp"

void test_2D_potentiel(){

	double rect_length_fe = 10;

	ZeroRightHandSide<2> rhs;

	for(unsigned nbr_of_strips=0; nbr_of_strips<=10; nbr_of_strips++){
		Rect2DBoundaryValues<2> bv(nbr_of_strips,1,1,rect_length_fe);
		LaplaceSolver<2> rect_potential_solver(&rhs,&bv,
				std::to_string(nbr_of_strips)+".vtk");
		rect_potential_solver.run();
	}

	/*
	Rect2DBoundaryValues<2> bv(1,1,1,rect_length_fe);
	LaplaceSolver<2> rect_potential_solver(&rhs,&bv, "one_strip.vtk");
	rect_potential_solver.run();

	//unsigned nbr_of_strip, unsigned strip_length,	unsigned pitch, double rect_length_fe
	Rect2DBoundaryValues<2> bv2(2,1,1,rect_length_fe);
	LaplaceSolver<2> rect_potential_solver2(&rhs,&bv2, "two_strip.vtk");
	rect_potential_solver2.run();*/
}

