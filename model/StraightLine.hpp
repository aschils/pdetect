/*
 * StraightLine.hpp
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
	 * @param: -alpha is the incline of the straight line
	 *		   -pass is a crossing point of our line
	 *		   -sol is simply the object where lies the needed informations
	 *		   -precision is the space between two consecutives points on the
	 *			line
	 */
	StraightLine(double alpha, Point<dim> const &pass, 
				Solution<dim> *sol, double precision);

	/**
	 * This functiun simply gives us the two extremeties of our line for a given
	 * incline and a given crossing point.
	 *
	 * /!\ For now it only works for alpha = 0 or alpha = 90°
	 */
	std::pair<Point<dim>, Point<dim>> get_extremeties(double alpha, 
											Point<dim> const &pass);
	void construct_line(double alpha, Point<dim> const &pass);

private:
	double precision;
	double rect_length_fe;
	double rect_width_fe;
	Solution<dim> *sol;
	std::vector<std::pair<ValuesAtPoint<dim>, Point<dim>>> values_on_line;
};

#include "StraightLine.cpp"

#endif