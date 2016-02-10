/*
 * MyGeometryInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

class MyGeometryInfo {

public:

	virtual bool is_point_inside_geometry(unsigned dim,
			std::vector<double> point_coord) = 0;

	//Return the maximum width of the domain
	virtual unsigned get_width() = 0;

	//Return the maximum length of the domain
	virtual unsigned get_length() = 0;

	~MyGeometryInfo() {
	}
};
