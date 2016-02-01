/*
 * MidRectRect2DDetector.cpp
 *
 *  Created on: 1 févr. 2016
 *      Author: aschils
 */

#include "MidRectRect2DDetector.hpp"


MidRectRect2DDetector::MidRectRect2DDetector(unsigned width, unsigned strip_length,
			unsigned strip_width, unsigned inter_strip_dist,
			unsigned nbr_of_strips, double potential, unsigned refine_level,
			unsigned max_iter, double stop_accuracy){
		this->width = width;
		this->strip_length = strip_length;
		this->strip_width = strip_width;
		this->inter_strip_dist = inter_strip_dist;
		this->nbr_of_strips = nbr_of_strips;
		this->strip_potential = potential;
		this->refine_level = refine_level;
		this->max_iter = max_iter;
		this->stop_accuracy = stop_accuracy;


}
