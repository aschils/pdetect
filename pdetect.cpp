/*
 * Author: Arnaud Schils, Simon Lardinois, Catholic University of Louvain, 2015
 */

#include "test/test_potential2D.hpp"
//#include <boost/program_options.hpp>

int main() {
	test_serrated_2D_potential();
	test_serrated_rect_limit_cases();
	test_electric_field();
	test_weighting_potential();
	test_various();

	//namespace po = boost::program_options;
	// Declare the supported options.
	/*po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "produce help message")
	    ("compression", po::value<int>(), "set compression level")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << desc << "\n";
	    return 1;
	}

	if (vm.count("compression")) {
	    cout << "Compression level was set to "
	 << vm["compression"].as<int>() << ".\n";
	} else {
	    cout << "Compression level was not set.\n";
	}*/

	return 0;
}
