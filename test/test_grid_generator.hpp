/*
 * test_grid_generator.hpp
 *
 *  Created on: 29 janv. 2016
 *      Author: aschils
 */

#pragma once

#include "../model/MyGridGenerator.hpp"
#include "../model/Utils.hpp"

void test_rectangle_with_circular_holes() {

	std::string output_dir = "tests_output_grid_generator/";
	Utils::create_directory_if_not_exists(output_dir);

	for (unsigned nbr_of_holes = 1; nbr_of_holes <= 3; nbr_of_holes++) {
		for (unsigned width = 100; width <= 300; width += 100) {
			for (unsigned holes_radius = 10; holes_radius <= 30; holes_radius +=
					10) {
				for (unsigned inter_holes_centers_dist = 3 * holes_radius;
						inter_holes_centers_dist <= 3 * holes_radius;
						inter_holes_centers_dist += holes_radius) {

					std::string output_file = "nbr_of_holes_"
							+ std::to_string(nbr_of_holes) + "_width_ "
							+ std::to_string(width) + "_holes_radius_"
							+ std::to_string(holes_radius)
							+ "_inter_holes_centers_dist_"
							+ std::to_string(inter_holes_centers_dist) + ".eps";

					try {
						Triangulation<2> tria;
						MyGridGenerator<2>::rectangle_with_circular_holes(tria,
								width, holes_radius, inter_holes_centers_dist,
								nbr_of_holes);
						GridOut grid_out;

						std::ofstream out(output_dir + output_file);
						grid_out.write_eps(tria, out);
					}

					catch (int e) {
						std::cout << output_file << " throws exception " << e
								<< std::endl;
					}
				}

			}
		}
	}
}

void test_rectangle_width_rectangular_holes() {

	unsigned width = 30;
	unsigned hole_length = 50;
	unsigned hole_width = 10;
	unsigned inter_holes_dist = 20;
	unsigned nbr_of_holes = 3;

	Triangulation<2> tria;
	MyGridGenerator<2>::rectangle_with_rectangular_holes(tria, width,
			hole_length, hole_width, inter_holes_dist, nbr_of_holes);
	GridOut grid_out;

	std::string output_dir = "tests_output_grid_generator/";
	Utils::create_directory_if_not_exists(output_dir);

	std::ofstream out(output_dir + "out.eps");
	grid_out.write_eps(tria, out);
}

