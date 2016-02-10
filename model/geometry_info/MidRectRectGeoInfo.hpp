/*
 * MidRectRectGeoInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "MyGeometryInfo.hpp"

class MidRectRectGeoInfo : public MyGeometryInfo {
	//TODO

	bool is_point_inside_geometry(unsigned dim,
			std::vector<double> point_coord) {
		return false;
	}

	//Return the maximum width of the domain
	unsigned get_width() {
		return 0;
	}

	//Return the maximum length of the domain
	unsigned get_length() {
		return 0;
	}

};
