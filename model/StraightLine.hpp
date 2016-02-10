/*
 * StraightLine.hpp
 *
 *  Created on: 14 déc. 2015
 *      Author: slardinois
 */

#pragma once

#include "Solution.hpp"
 #include "boundary_conditions/SerratedRect2DBoundaryValues.hpp"
#include <fstream>
#include <iostream>

 #define PI 3.14159265

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
	Point<dim> get_beginning(double alpha, Point<dim> const &pass);
	//bool is_strip(const Point<dim> &p, unsigned strip_width, 
	//			unsigned strip_length, unsigned half_pitch);
	ValuesAtPoint<dim> exact_solution(Point<dim> const &point);
	void construct_line(double alpha, Point<dim> const &pass);
	void write_data_file();

private:
	double precision;
	double rect_length_fe;
	double rect_width;
	Solution<dim> *sol;
	std::vector<std::pair<ValuesAtPoint<dim>, Point<dim>>> values_on_line;
	std::vector<std::pair<ValuesAtPoint<dim>, Point<dim>>> exact_values_on_line;
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
	rect_width = sol->values_at_cells[l].first->vertex(3)[1];

	construct_line(alpha, pass);
}

template<unsigned dim>
Point<dim> StraightLine<dim>::get_beginning(double alpha,
														Point<dim> const &pass) {
	Point<dim> begin;

	begin[1] = 0;

	if(alpha != 90 && alpha != 0){
	//We check the extrem points by simple trigonometry
		begin[0] = pass[0] - pass[1] / tan(alpha);
		if(begin[0] < 0){
			begin[0] = 0;
			begin[1] = pass[1] - pass[0]*tan(alpha);
		}
	}
	else if (alpha == 0){
		begin[0] = 0;
		begin[1] = pass[1];
	}
	else {
		begin[0] = pass[0];
		begin[1] = 0;
	}

	return begin;
}

template<unsigned dim>
ValuesAtPoint<dim> StraightLine<dim>::exact_solution(Point<dim> const &point){
	ValuesAtPoint<dim> exact_value;

	double epsilon = 0.000000001;
	double x = point[0]-rect_length_fe/2;
	double y = -point[1]/rect_width + 1;
	double pot;
	pot = sin(PI*y)*sinh(PI*50);
	pot = pot/ (cosh(PI*x)-cos(PI*y)*cosh(PI*50));
	if(point[1] >= rect_width/2 - epsilon)
		pot = (atan(pot)+PI)/PI;
	else
		pot = (atan(pot))/PI;

	exact_value.fun = pot;

	return exact_value;
}

template <unsigned dim>
void StraightLine<dim>::construct_line(double alpha, Point<dim> const &pass) {

	Point<dim> point = get_beginning(alpha, pass);

	while(point[0] <= rect_length_fe && point[1] <= rect_width) {

		bool strip = SerratedRect2DBoundaryValues<dim>::is_strip(point, rect_width, 50, 100, 50);

		if(!strip){
			ValuesAtPoint<dim> value = sol->get_values(point);

			std::pair<ValuesAtPoint<dim>, Point<dim>> values_at_point;

			values_at_point.first = value;
			values_at_point.second = point;

			values_on_line.push_back(values_at_point);

			ValuesAtPoint<dim> exact_value = exact_solution(point);

			std::pair<ValuesAtPoint<dim>, Point<dim>> exact_values_at_point;

			exact_values_at_point.first = exact_value;
			exact_values_at_point.second = point;

			exact_values_on_line.push_back(exact_values_at_point);
		}

		point[0] = point[0] + precision*cos(alpha);
		point[1] = point[1] + precision*sin(alpha);
	}

	write_data_file();
}

template<unsigned dim>
void StraightLine<dim>::write_data_file() {

	std::ofstream gnu_graph;
	gnu_graph.open("sol_on_line", std::fstream::in | std::fstream::out | std::fstream::trunc);

	if(gnu_graph.is_open()){

		gnu_graph << "# X\tY" << std::endl;
		for(unsigned i = 0; i < values_on_line.size(); i++){
			
			gnu_graph << values_on_line[i].second[1] << "\t"
					  << values_on_line[i].first.fun << std::endl;;

		}
	}
	gnu_graph.close();

	std::ofstream agnu_graph;
	agnu_graph.open("exact_sol_on_line", std::fstream::in | std::fstream::out | std::fstream::trunc);

	if(agnu_graph.is_open()){

		agnu_graph << "# X\tY" << std::endl;
		for(unsigned i = 0; i < exact_values_on_line.size(); i++){
			
			agnu_graph << exact_values_on_line[i].second[1] << "\t"
					  << exact_values_on_line[i].first.fun << std::endl;;

		}
	}
	agnu_graph.close();
}