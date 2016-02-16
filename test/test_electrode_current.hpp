/*
 * test_electrode_current.hpp
 *
 *  Created on: 11 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "../model/ElectrodeCurrent.hpp"

void test_electrode_current() {

	double strip_potential = 90; //typiquement 100V
	unsigned nbr_of_strips = 1;
	unsigned width = 10;
	unsigned strip_length = 100;
	unsigned strip_width = 0;
	unsigned half_pitch = 0;
	unsigned refine_level = 5;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	//std::string output_dir = "tests_electrode_current";
	//Utils::create_directory_if_not_exists(output_dir);

	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, strip_potential, refine_level, max_iter, stop_accuracy);
	srdd.compute();
	srdd.compute_weight();
	srdd.draw_vtk_graph_potential("electrode_pot.vtk");
	srdd.draw_vtk_graph_weight_potential("electrode_pot_weight.vtk");
	srdd.draw_vtk_graph_gradient_of_potential("electrode_pot_grad.vtk");

	Solution<2> solution = srdd.get_solution();
	Solution<2> weight_solution = srdd.get_solution_weight();
	MyGeometryInfo *geo_info = srdd.get_geometry_info();

	Point<2> p1(50,0.0);
	Point<2> p2(50, 100);
	Segment seg(p1,p2);

	Line particle_traj(seg);

	ElectrodeCurrent<2> ec(geo_info, &solution, &weight_solution,
			&particle_traj, 10);
	//ec.print_charges();
	double delta_t = 0.0000000000001; //100ps p.76, V_b = 100V, v_d = 30V
	std::vector<std::pair<double, double> > current_vs_time;
	ec.compute_current(delta_t, current_vs_time);

	std::cout << current_vs_time[0].second << " " << std::endl;

	std::string output_graph = "current";

	Utils::write_gnu_data_file<2>(output_graph, current_vs_time);

}

