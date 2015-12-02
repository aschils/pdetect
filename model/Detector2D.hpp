/*
 * Detector2D.hpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

#ifndef __DETECTOR2D_HPP__
#define __DETECTOR2D_HPP__

#include <deal.II/numerics/data_out.h>

#include "Solution.hpp"
#include "Gradient.hpp"

class Detector2D {

public:
	virtual Solution<2> compute_potential() = 0;
	virtual Solution<2> compute_weighting_potential() = 0;
	virtual DataOut<2> compute_electric_field() = 0;
	virtual ~Detector2D(){}
};

#endif
