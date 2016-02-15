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
	unsigned half_pitch = 50;
	unsigned refine_level = 5;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, strip_potential, refine_level, max_iter, stop_accuracy);
	srdd.compute();

	Solution<2> solution = srdd.get_solution();
	MyGeometryInfo *geo_info = srdd.get_geometry_info();

	Point<2> pass(100, 0);
	double precision = 0.1;
	StraightLine<2> line(PI/2, pass, &solution, geo_info, precision);
}