/*
 * MidRectRect2DDetector.hpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "Detector2D.hpp"
#include "../MyGridGenerator.hpp"

class MidRectRect2DDetector: public Detector2D {


	MidRectRect2DDetector(unsigned width, unsigned strip_length,
			unsigned strip_width, unsigned inter_strip_dist,
			unsigned nbr_of_strips, double potential, unsigned refine_level,
			unsigned max_iter, double stop_accuracy);

	std::string params_to_string();

private:
	unsigned width, strip_length, strip_width;
	unsigned inter_strip_dist, nbr_of_strips;
};
