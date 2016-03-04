/*
 * test_straight_line.hpp
 *
 *  Created on: 15 févr. 2016
 *      Author: slardinois
 */

#pragma once

#include "../model/StraightLine.hpp"

void test_straight_line() {

	double strip_potential = 1;
	unsigned nbr_of_strips = 1;
	unsigned width = 300;
	unsigned strip_length = 100;
	unsigned strip_width = 50;
	unsigned half_pitch = 150;
	unsigned max_iter = 100000;
	double stop_accuracy = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	SerratedRect2DDetector srdd(nbr_of_strips,
			width, strip_length, strip_width, half_pitch, strip_potential,
			refine_accuracy, max_iter, stop_accuracy);
	srdd.compute_weight();
	srdd.draw_vtk_graph_weight_potential("weighting_pot.vtk");

	Solution<2> solution;
	srdd.get_solution_weight(solution);
	SerratedRectGeoInfo *geo_info =
			(SerratedRectGeoInfo*) srdd.get_geometry_info();

	double x = half_pitch + strip_length / 2;
	Point<2> pass(x, 0);
	double precision = 0.01;
	StraightLine<2> line(M_PI / 2, pass, &solution, geo_info, precision);

	std::vector<std::pair<double, double>> data = line.get_data();
	std::vector<std::pair<double, double>> exact_data = line.get_exact_data();
	//std::vector<std::pair<double, double>> ratio = line.get_ratio();

	Utils::write_gnu_data_file<2>("sol_on_line", data);
	Utils::write_gnu_data_file<2>("exact_sol_on_line", exact_data);
	//Utils::write_gnu_data_file<2>("ratio_on_line", ratio);
}
