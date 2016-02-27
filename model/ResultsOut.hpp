/*
 * ResultsOut.hpp
 *
 *  Created on: 27 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "Utils.hpp"

class ResultsOut {

public:

	static void write_current_vs_time(std::string file,
			std::vector<std::pair<double, double> > &current_vs_time){
		Utils::write_vector_of_pair<double, double>(file, current_vs_time);
	}

};

