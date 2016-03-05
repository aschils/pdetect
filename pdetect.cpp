/*
 * Author: Arnaud Schils, Simon Lardinois, Catholic University of Louvain, 2015
 */

#include "test/test_potential2D.hpp"
#include "test/test_grid_generator.hpp"
#include "test/test_electrode_current.hpp"
#include "test/test_straight_line.hpp"
#include "test/test_cases_projet.hpp"
//#include "view/cli/Cli.hpp"

int main(int argc, char* argv[]) {
    //test_serrated_2D_potential();
	//test_serrated_rect_limit_cases();
//	test_electric_field();
	//test_weighting_potential();
//	test_various();
	//test_mid_circle_rect2D_det();
	//test_mid_rect_rect_2D_det();
//
//	test_rectangle_with_circular_hole();
//	test_rectangle_width_rectangular_holes();
	//test_serrated_rectangle_gaussian();

	test_electrode_current_serrated();
	//test_electrode_current_mid_rect_rect();
	//gen_comparison_data();
	//test_straight_line();

	//test_case_silicium();
	//test_case_helium();
	//test_case_helium_scd();

	//Cli cli(argc, argv);
	//cli.parse_cmd();
	return 0;
}
