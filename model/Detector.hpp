/*
 * Detector.hpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

#ifndef __DETECTOR_HPP__
#define __DETECTOR_HPP__

class Detector {

public:
	virtual void compute_potential() = 0;
	virtual ~Detector(){}
};

#endif