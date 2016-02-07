/*
 * StraightLine.hpp
 *
 *  Created on: 14 déc. 2015
 *      Author: slardinois
 */

#pragma once

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
	 * @param: -alpha is the tilt of the straight line
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
	 */
	std::pair<Point<dim>, Point<dim>> get_beginning(double alpha, 
											Point<dim> const &pass);
	void construct_line(double alpha, Point<dim> const &pass);

private:
	double precision;
	double rect_length_fe;
	double rect_width_fe;
	Solution<dim> *sol;
	std::vector<std::pair<ValuesAtPoint<dim>, Point<dim>>> values_on_line;
};



template<unsigned dim>
StraightLine<dim>::StraightLine(double alpha, Point<dim> const &pass,
								Solution<dim> *sol,
								double precision) {

	this->precision = precision;
	this->sol = sol;

	//Here we get the the size of the detector by looking the coordinates of the 
	//last point in our finite element domain
	int l = sol->values_at_cells.size() - 1;
	rect_length_fe = sol->values_at_cells[l].first->vertex(3)[0];
	rect_width_fe = sol->values_at_cells[l].first->vertex(3)[1];

	construct_line(alpha, pass);
}

template<unsigned dim>
std::pair<Point<dim>, Point<dim>> StraightLine<dim>::get_beginning(double alpha,
														Point<dim> const &pass) {
	Point<dim> begin;

	begin[1] = 0;

	//We check the extrem points by simple trigonometry
	begin[0] = pass[0] - pass[1] / tan(alpha);
	if(begin[0] < 0){
		begin[0] = 0;
		begin[1] = pass[1] - pass[0]*tan(alpha);
	}

	return begin;
}

template <unsigned dim>
void StraightLine<dim>::construct_line(double alpha, Point<dim> const &pass) {

	Point<dim> point = get_beginning(alpha, pass);
	while(point[0] <= rect_length_fe && point[1] <= rect_width_fe) {

		std::cout << point[0] << std::endl
				<< point[1] << std::endl;
		ValuesAtPoint<dim> value = sol->get_values(point);

		std::pair<ValuesAtPoint<dim>, Point<dim>> values_at_point;

		values_at_point.first = value;
		values_at_point.second = point;

		values_on_line.push_back(values_at_point);

		point[0] = point[0] + precision*cos(alpha);
		point[1] = point[1] + precision*sin(alpha);
	}
}

