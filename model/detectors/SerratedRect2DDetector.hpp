/*
 * SerratedRect2DDetector.hpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#pragma once

#include <unordered_map>

#include "Detector2D.hpp"
#include <functional>
#include "../TensorUtils.hpp"
#include "../boundary_conditions/SerratedRect2DBoundaryCond.hpp"
#include "../boundary_conditions/SerratedRect2DBoundaryCondWeight.hpp"
#include "../geometry_info/SerratedRectGeoInfo.hpp"

#define DEFAULT_RECT_WIDTH 300.0 //i.e. in domain language (microm,..)

using namespace std;

class SerratedRect2DDetector : public Detector2D {

public:
	SerratedRect2DDetector(unsigned nbr_of_strips,
			unsigned strip_length, unsigned strip_width, unsigned half_pitch,
			double strip_potential, unsigned refine_level, unsigned max_iter,
			double stop_accuracy);

	SerratedRect2DDetector(unsigned nbr_of_strips, unsigned width,
				unsigned strip_length, unsigned strip_width, unsigned half_pitch,
				double strip_potential, unsigned refine_level, unsigned max_iter,
				double stop_accuracy);

	std::vector<double> get_electric_field(Point<2> p);

	std::string params_to_string();

private:

	unsigned nbr_of_pts_along_x = 0, nbr_of_pts_along_y = 0;

	unsigned compute_total_length();
	SerratedRectGeoInfo *serr_geo_info;
	//void nbr_of_points_along_axes();
};

