/*
 * StraightLine.hpp
 *
 *  Created on: 14 déc. 2015
 *      Author: slardinois
 */

#pragma once

#include "Solution.hpp"
#include "boundary_conditions/SerratedRect2DBoundaryValues.hpp"
#include "Utils.hpp"
#include "geometry_info/MyGeometryInfo.hpp"
#include <iostream>
#include <math.h>

//#define M_PI 3.14159265

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
				Solution<dim> *sol, 
				SerratedRectGeoInfo *geo_info,
				double precision);

	/**
	 * This functiun simply gives us the two extremeties of our line for a given
	 * incline and a given crossing point.
	 */
	Point<dim> get_beginning(double alpha, Point<dim> const &pass);
	//bool is_strip(const Point<dim> &p, unsigned strip_width, 
	//			unsigned strip_length, unsigned half_M_PItch);
	PhysicalValues<dim> exact_solution(Point<dim> const &point);
	void construct_line(double alpha, Point<dim> const &pass);
	std::vector<std::pair<double, double>> get_data();
	std::vector<std::pair<double, double>> get_exact_data();
	std::vector<std::pair<double, double>> get_ratio();

private:
	unsigned strip_length, strip_width, detector_width, detector_length;
	double precision;
	Solution<dim> *sol;
	SerratedRectGeoInfo *geo_info;
	std::vector<std::pair<PhysicalValues<dim>, Point<dim>>> values_on_line;
	std::vector<std::pair<double, double>> solution_data;
	std::vector<std::pair<PhysicalValues<dim>, Point<dim>>> exact_values_on_line;
	std::vector<std::pair<double, double>> exact_solution_data;

};



template<unsigned dim>
StraightLine<dim>::StraightLine(double alpha, Point<dim> const &pass,
								Solution<dim> *sol,
								SerratedRectGeoInfo *geo_info,
								double precision) {

	this->precision = precision;
	this->sol = sol;
	this->geo_info = geo_info;

	detector_length = geo_info->get_length();
	detector_width = geo_info->get_width();
	strip_length = geo_info->get_strip_length();
	strip_width = geo_info->get_strip_width();

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
PhysicalValues<dim> StraightLine<dim>::exact_solution(Point<dim> const &point){
	PhysicalValues<dim> exact_value;

	double x = point[0]-detector_length/2;
	double y = -point[1]/(detector_width-strip_width) + 1;
	double pot;

	pot = sin(M_PI*y)*sinh(M_PI*strip_length/(2*detector_width));
	pot = pot/ (cosh(M_PI*x)-cos(M_PI*y)*cosh(M_PI*strip_length/(2*detector_width)));

	if(atan(pot) <= 0)
		pot = (atan(pot)+M_PI)/M_PI;
	else
		pot = (atan(pot))/M_PI;

	exact_value.potential = pot;

	return exact_value;
}

template <unsigned dim>
void StraightLine<dim>::construct_line(double alpha, Point<dim> const &pass) {

	Point<dim> point = get_beginning(alpha, pass);

	double epsilon = 0.0001;
	while(Utils::less_than_or_equals_double(point[0],detector_length, epsilon) && 
			Utils::less_than_or_equals_double(point[1],detector_width, epsilon)) {

		bool strip = geo_info->is_strip<dim>(point);

		if(!strip){
			PhysicalValues<dim> value = sol->get_values(point);
			std::pair<PhysicalValues<dim>, Point<dim>> values_at_point(value, point);

			values_on_line.push_back(values_at_point);
			solution_data.push_back(std::pair<double, double>(point[1], value.potential));


			PhysicalValues<dim> exact_value = exact_solution(point);
			std::pair<PhysicalValues<dim>, Point<dim>> exact_values_at_point(exact_value, point);

			exact_values_on_line.push_back(exact_values_at_point);
			exact_solution_data.push_back(std::pair<double, double>(point[1], exact_value.potential));
		}
		point[0] = point[0] + precision*cos(alpha);
		point[1] = point[1] + precision*sin(alpha);
	}
}

template<unsigned dim>
std::vector<std::pair<double, double>> StraightLine<dim>::get_data() {

	return solution_data;
}

template<unsigned dim>
std::vector<std::pair<double, double>> StraightLine<dim>::get_exact_data() {

	return exact_solution_data;
}

template<unsigned dim>
std::vector<std::pair<double, double>> StraightLine<dim>::get_ratio() {

	std::vector<std::pair<double, double>> ratio;
	std::pair<double, double> point_ratio;

	for(unsigned i = 0; i < solution_data.size(); i++) {

		if(!Utils::equals_double(exact_solution_data[i].second, 0, 0.000001))
			point_ratio.second = solution_data[i].second /
								exact_solution_data[i].second;
		else
			point_ratio.second = 1;
		point_ratio.first = solution_data[i].first;
		ratio.push_back(point_ratio);
	}

	return ratio;
}