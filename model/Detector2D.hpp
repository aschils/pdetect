/*
 * Detector.hpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

#ifndef __DETECTOR_HPP__
#define __DETECTOR_HPP__

#include "Solution.hpp"

class Detector2D {

public:
	virtual Solution<2> compute_potential() = 0;
	virtual Solution<2> compute_weighting_potential() = 0;
	virtual Solution<2> compute_electric_field() = 0;
	virtual ~Detector2D(){}
};

#endif
