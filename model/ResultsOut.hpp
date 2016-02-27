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

		//Convert seconds in nanoseconds

		std::vector<std::pair<double, double> >
		current_vs_nanosec(current_vs_time.size());

		for(unsigned i=0; i<current_vs_time.size(); i++){
			std::pair<double, double> p(current_vs_time[i].first*1e9,
					current_vs_time[i].second);
			current_vs_nanosec[i] = p;
		}

		Utils::write_vector_of_pair<double, double>(file, current_vs_nanosec,
				true, "TIME[ns],", "Itot[uA],");
	}

};

