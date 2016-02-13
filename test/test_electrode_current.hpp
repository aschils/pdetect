/*
 * test_electrode_current.hpp
 *
 *  Created on: 11 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "../model/ElectrodeCurrent.hpp"

void test_electrode_current() {

	double strip_potential = 1;
	unsigned nbr_of_strips = 1;
	unsigned width = 300;
	unsigned strip_length = 100;
	unsigned strip_width = 50;
	unsigned half_pitch = 50;
	unsigned refine_level = 5;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	//std::string output_dir = "tests_electrode_current";
	//Utils::create_directory_if_not_exists(output_dir);

	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, strip_potential, refine_level, max_iter, stop_accuracy);
	srdd.compute();
	srdd.compute_weight();

	Solution<2> potential = srdd.get_solution_potential();
	Solution<2> weight_potential = srdd.get_solution_weight_potential();
	MyGeometryInfo *geo_info = srdd.get_geometry_info();

	std::vector<
			std::pair<typename DoFHandler<2>::active_cell_iterator,
					std::vector<Tensor<1, 2> > > > electric_field =
			srdd.get_electric_field();

	std::vector<
			std::pair<typename DoFHandler<2>::active_cell_iterator,
					std::vector<Tensor<1, 2> > > > electric_field_weight =
			srdd.get_electric_field_weight();


	Point<2> p1(100,0.0);
	Point<2> p2(100, 100);
	Segment seg(p1,p2);

	Line particle_traj(-1, 300);

	ElectrodeCurrent<2> ec(geo_info, &potential, &weight_potential, &electric_field,
			&electric_field_weight, &particle_traj, 7);
	ec.print_charges();
}

