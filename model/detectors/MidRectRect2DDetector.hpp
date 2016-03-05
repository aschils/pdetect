/*
 * MidRectRect2DDetector.hpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "Detector2D.hpp"
#include "../MyGridGenerator.hpp"
#include "../boundary_conditions/MidRectRect2DBoundaryCond.hpp"
#include "../boundary_conditions/MidRectRect2DBoundaryCondWeight.hpp"
#include "../geometry_info/MidRectRectGeoInfo.hpp"

class MidRectRect2DDetector: public Detector2D {

public:
	MidRectRect2DDetector(unsigned half_width, unsigned strip_length,
			unsigned half_strip_width, unsigned half_inter_strip_dist,
			unsigned nbr_of_strips, double potential, double refine_accuracy,
			unsigned max_iter, double stop_accuracy, unsigned material_id);

	std::string params_to_string();

private:
	MidRectRectGeoInfo *mrr_geo_info;
};
