/*
 * MyGeometryInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "../Utils.hpp"

class MyGeometryInfo {

public:

	virtual bool is_point_inside_geometry(std::vector<double> point_coord) = 0;

	//Returns the segments part of the line which are inside the detector
	//domain
	virtual std::vector<Utils::Segment<2>> segments_at_intersec(
			Utils::Line<2> line) = 0;

	//Return the maximum width of the domain
	virtual unsigned get_width() = 0;

	//Return the maximum length of the domain
	virtual unsigned get_length() = 0;

	virtual unsigned get_nbr_of_strips() = 0;

	virtual unsigned get_dimension() = 0;

	~MyGeometryInfo() {
	}
};
