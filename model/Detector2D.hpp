/*
 * Detector2D.hpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

#ifndef __DETECTOR2D_HPP__
#define __DETECTOR2D_HPP__

#include <deal.II/numerics/data_out.h>

#include "Derivatives.hpp"
#include "Solution.hpp"

class Detector2D {

public:
	virtual void compute() = 0;
	//virtual SolutionScalar<2> compute_weighting_potential() = 0;
	//virtual SolutionVector<2> compute_gradient_of_potential() = 0;
	virtual ~Detector2D(){}
};

#endif
