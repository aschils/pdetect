/*
 * Detector2D.hpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

#ifndef __DETECTOR2D_HPP__
#define __DETECTOR2D_HPP__

#include "Solution.hpp"
#include "Gradient.hpp"

class Detector2D {

public:
	virtual Solution<2> compute_potential() = 0;
	virtual Solution<2> compute_weighting_potential() = 0;
	virtual void compute_electric_field(std::string gradient_file_path) = 0;
	virtual ~Detector2D(){}
};

#endif
