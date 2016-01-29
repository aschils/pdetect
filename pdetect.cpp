/*
 * Author: Arnaud Schils, Simon Lardinois, Catholic University of Louvain, 2015
 */

#include "test/test_potential2D.hpp"
//#include "view/cli/Cli.hpp"

int main(int argc, char** argv) {
	test_serrated_2D_potential();
	test_serrated_rect_limit_cases();
	test_electric_field();
	test_weighting_potential();
	test_various();

	//Cli cli(argc, argv);
	//cli.parse_cmd();
	return 0;
}
