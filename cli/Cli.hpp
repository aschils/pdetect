/*
 * Cli.hpp
 *
 *  Created on: 28 janv. 2016
 *      Author: aschils
 */

#pragma once

#include <boost/program_options.hpp>

class Cli {

public:

	int argc;
	char **argv;

	Cli(int argc, char** argv) {
		this->argc = argc;
		this->argv = argv;
	}

	void interpret_cmd() {
		namespace po = boost::program_options;
		// Declare the supported options.
		po::options_description desc("Allowed options");
		desc.add_options()("help", "produce help message")("compression",
				po::value<int>(), "set compression level");

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			//TODO THROW EXCEPTION
		}

		if (vm.count("compression")) {
			cout << "Compression level was set to "
					<< vm["compression"].as<int>() << ".\n";
		} else {
			cout << "Compression level was not set.\n";
		}

	}

};

