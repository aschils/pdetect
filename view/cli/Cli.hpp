/*
 * Cli.hpp
 *
 *  Created on: 28 janv. 2016
 *      Author: aschils
 */

#pragma once

#include <boost/program_options.hpp>
#include <unordered_map>
#include <sstream>

#include "../../model/errors.hpp"

namespace po = boost::program_options;

class Cli {

public:

	int argc;
	char **argv;
	unsigned nbr_of_detector_types = 1;

	std::unordered_map<std::string, std::string> detector_types = { { "sr",
			"serrated rectangular" } };

	Cli(int argc, char** argv) {
		this->argc = argc;
		this->argv = argv;
	}

	void parse_cmd() {

		po::options_description desc("Allowed options");
		add_options_to_desc(desc);

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << desc << "\n";
			return;
		}

		validate_input(vm); //This throws Exception
	}

private:

	std::string detector_types_to_str() {

		std::stringstream ss;

		for (auto it = detector_types.begin();
				it != detector_types.end();
				it++)
			ss << " " << it->first << " (" << it->second << ")";

		ss << std::endl;
		return ss.str();
	}

	void add_options_to_desc(po::options_description &desc){

		std::string detector_types_str = detector_types_to_str();

		desc.add_options()("help,h", "produce help message")
						("detector-type,d", po::value<std::string>()->default_value("sr"),
								std::string("detector type: ").append(detector_types_str).c_str())
						("number-of-strips,n", po::value<unsigned>()->default_value(3),
								"number of strips in the detector")
						("strip-potential,v", po::value<unsigned>()->default_value(1),
								"strip potential (Volt)")
						("strip-length,l", po::value<unsigned>()->default_value(100),
								"strip length (nm)")
						("strip-width,w", po::value<unsigned>()->default_value(30),
								"strip width (nm)")
						("pitch,p", po::value<unsigned>()->default_value(100),
								"pitch: space between two consecutive strips (nm)")
						("refine-level,r", po::value<unsigned>()->default_value(1),
								"refine level: mesh precision increases exponentially with the refine level as 2^refine_level")
						("error,e", po::value<double>()->default_value(10e-12),
										"tolerance determining the success of the iteration: computation stops when error less or equal than specified error")
						("max-iteration,i", po::value<unsigned>()->default_value(10000),
										"maximum number of iterations allowed to reduce the error until user-specified error is reached (exception thrown if required precision not reached)")
						("output-file,o", po::value<std::string>()->default_value("out.vtk"),
										"results are outputed as a vtk file at specified path");
	}

	/**
	 * @throw INVALID_INPUT_EXCEPTION  if at least one user-provided input is invalid,
	 * does nothing otherwise.
	 */
	void validate_input(po::variables_map &vm){
		std::string detec_type = vm["detector-type"].as<std::string>();

		//Check if specified detector type exists
		try{
			detector_types.at(detec_type);
		}
		catch (out_of_range &ofr)  {
			std::cout << "unknown detector type" << std::endl;
		    throw INVALID_INPUT;
		}

		//Check if accuracy is a positive number
		double accuracy = vm["error"].as<double>();

		if(accuracy < 0){
			std::cout << "error is not allowed to be negative" << std::endl;
			throw INVALID_INPUT;
		}
	}
};

