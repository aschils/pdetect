/*
 * potential2D.cpp
 *
 *  Created on: 8 nov. 2015
 *      Author: aschils
 */

#include <sys/stat.h>

#include "../model/ZeroRightHandSide.hpp"
#include "../model/SerratedRect2DDetector.hpp"

void test_serrated_2D_potential() {

	ZeroRightHandSide<2> rhs;
	double strip_potential = 1;
	unsigned strip_length = 100;
	unsigned strip_width = 50;
	unsigned pitch = 100;
	unsigned refine_level = 5;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	std::string output_dir = "tests_output";

	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 10; nbr_of_strips++) {
		std::string output_file = output_dir + "/" + "pitch_"
				+ std::to_string(pitch) + "_" + std::to_string(nbr_of_strips)
				+ ".vtk";
		SerratedRect2DDetector srdd(nbr_of_strips, strip_length, strip_width,
				pitch, strip_potential, refine_level, max_iter, stop_accuracy);
		srdd.compute_potential(output_file);
	}
}

void test_serrated_rect_limit_cases() {

	unsigned max_nbr_of_strips = 5;

	ZeroRightHandSide<2> rhs;
	double strip_potential = 1;
	unsigned strip_length = 100;
	unsigned width = 300;
	unsigned strip_width = 300;
	unsigned pitch = 50;
	unsigned refine_level = 7;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	std::string output_dir = "tests_output/";
	/*
	 //Limit case: hole_width = width

	 for (unsigned nbr_of_strips = 0; nbr_of_strips <= max_nbr_of_strips;
	 nbr_of_strips++) {
	 std::string output_file = output_dir + "width_equals_hole_width_"
	 + std::to_string(nbr_of_strips) + ".vtk";
	 SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length,
	 strip_width, pitch, strip_potential, refine_level, max_iter,
	 stop_accuracy);
	 srdd.compute_potential(output_file);
	 }

	 //Limit case: strip_length = 0

	 strip_width = 80;
	 strip_length = 0;
	 for (unsigned nbr_of_strips = 0; nbr_of_strips <= max_nbr_of_strips;
	 nbr_of_strips++) {
	 std::string output_file = output_dir + "strip_length_0_"
	 + std::to_string(nbr_of_strips) + ".vtk";
	 SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length,
	 strip_width, pitch, strip_potential, refine_level, max_iter,
	 stop_accuracy);
	 srdd.compute_potential(output_file);
	 }

	 //Limit case: strip_width = 0
	 refine_level = 7;
	 strip_length = 100;
	 strip_width = 0;
	 for (unsigned nbr_of_strips = 0; nbr_of_strips <= max_nbr_of_strips;
	 nbr_of_strips++) {
	 std::string output_file = output_dir + "strip_width_0_"
	 + std::to_string(nbr_of_strips) + ".vtk";
	 SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length,
	 strip_width, pitch, strip_potential, refine_level, max_iter,
	 stop_accuracy);
	 srdd.compute_potential(output_file);
	 }

	 //pitch = 0

	 strip_width = 80;
	 pitch = 0;
	 for (unsigned nbr_of_strips = 0; nbr_of_strips <= max_nbr_of_strips;
	 nbr_of_strips++) {
	 std::string output_file = output_dir + "pitch_0_"
	 + std::to_string(nbr_of_strips) + ".vtk";
	 SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length,
	 strip_width, pitch, strip_potential, refine_level, max_iter,
	 stop_accuracy);
	 srdd.compute_potential(output_file);
	 }*/

	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 3; nbr_of_strips++) {
		for (unsigned strip_length = 0; strip_length <= 50; strip_length +=
				50) {
			for (unsigned strip_width = 0; strip_width <= 30; strip_width +=
					30) {
				for (unsigned pitch = 0; pitch <= 50; pitch += 50) {

					SerratedRect2DDetector srdd(nbr_of_strips, width,
							strip_length, strip_width, pitch, strip_potential,
							refine_level, max_iter, stop_accuracy);
					std::string output_file = output_dir
							+ srdd.params_to_string() + ".vtk";
					srdd.compute_potential(output_file);

				}
			}
		}
	}

	//Limit case: hole_width = width && strip_length = 0

	//Limit case: hole_width = width && strip_length = 0
}

void test_electric_field() {
	ZeroRightHandSide<2> rhs;
	double strip_potential = 1;
	unsigned strip_length = 100;
	unsigned strip_width = 50;
	unsigned pitch = 100;
	unsigned refine_level = 5;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;
	unsigned nbr_of_strips = 3;
	unsigned width = 300;

	std::string output_dir = "tests_output/";

	std::string output_file = output_dir + "pitch_" + std::to_string(pitch)
			+ "_" + std::to_string(nbr_of_strips) + ".vtk";
	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			pitch, strip_potential, refine_level, max_iter, stop_accuracy);
	srdd.compute_potential(output_file);
	srdd.compute_electric_field(output_dir + "EE.vtk");
}

void test_weighting_potential() {
	ZeroRightHandSide<2> rhs;
	double strip_potential = 100;
	unsigned strip_length = 100;
	unsigned strip_width = 50;
	unsigned pitch = 100;
	unsigned refine_level = 5;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	std::string output_dir = "tests_output";

	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 10; nbr_of_strips++) {
		std::string output_file = output_dir + "/" + "pitch_"
				+ std::to_string(pitch) + "_" + std::to_string(nbr_of_strips)
				+ ".vtk";
		SerratedRect2DDetector srdd(nbr_of_strips, strip_length, strip_width,
				pitch, strip_potential, refine_level, max_iter, stop_accuracy);
		srdd.compute_weighting_potential(output_file);
	}
}

