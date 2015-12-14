/*
 * StraightLine.cpp
 *
 *  Created on: 14 déc. 2015
 *      Author: slardinois
 */

#ifndef __STRAIGHT_LINE_HPP__
#define __STRAIGHT_LINE_HPP__

#include "Solution.hpp"

using namespace dealii;

/**
 * This class will evaluate the values of the potential and the elecric field
 * along a straight line.
 * It will then plot the results in a graph to compare it with the known
 * theorical results to compare them.
 */

template<unsigned dim>
class StraightLine {
public:
	/**
	 * @param: alpha is the incline of the straight line
	 *		   pass is a point where our line pass
	 *		   sol is simply the object where lies the needed informations
	 */
	StraightLine(double alpha, Point<dim> pass, Solution<dim> sol);
	void construct_line(double alpha, Point<dim>);

private:
	Solution<dim> sol;
	std::vector<std::pair<ValuesAtPoint<dim>, Point<dim>>> values_on_line;
};

#endif