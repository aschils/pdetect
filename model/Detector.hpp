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
	virtual void compute_potential(std::string output_file) = 0;
	virtual void compute_electric_field(std::string output_file) = 0;
	virtual ~Detector(){}
};

#endif
