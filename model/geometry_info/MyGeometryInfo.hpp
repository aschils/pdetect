/*
 * MyGeometryInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "../geometry_info/BasicGeometry.hpp"

class MyGeometryInfo {

public:

	virtual bool is_point_inside_geometry(Point<2> p) = 0;

	//Returns the intersection points between the particle trajectory
	//and the domain boundaries
	virtual std::vector<Point<2>> boundaries_intersections(
			Line line) = 0;

	//Return the maximum width of the domain
	virtual unsigned get_width() = 0;

	//Return the maximum length of the domain
	virtual unsigned get_length() = 0;

	virtual unsigned get_nbr_of_strips() = 0;

	virtual unsigned get_dimension() = 0;

	virtual Line get_mid_length_vertical_line() = 0;

	~MyGeometryInfo() {
	}
};
