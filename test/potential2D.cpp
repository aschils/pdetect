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
	ZeroRightHandSide<2> rhs;
	Rect2DBoundaryValues<2> bv;
	LaplaceSolver<2> rect_potential_solver(&rhs,&bv);
	rect_potential_solver.run();
}

