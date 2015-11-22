/*
 * potential2D.cpp
 *
 *  Created on: 8 nov. 2015
 *      Author: aschils
 */

#include "../model/Rect2DBoundaryValues.hpp"
#include "../model/ZeroRightHandSide.hpp"
#include "../model/Rect2DDetector.hpp"
#include "../model/SerratedRect2DDetector.hpp"

void test_2D_potential() {

	double strip_potential = 1;
	/*unsigned strip_length = 100;
	 unsigned pitch = 100;*/
	unsigned refine_level = 3;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	std::string output_dir = "tests_output";

	/*
	 for (unsigned nbr_of_strips = 0; nbr_of_strips <= 5; nbr_of_strips++) {
	 Rect2DDetector rdd(nbr_of_strips, strip_length, pitch, strip_potential,
	 refine_level, max_iter, stop_accuracy);
	 rdd.compute_potential();
	 }
	 */
	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 3; nbr_of_strips++) {
		for (unsigned strip_length = 0; strip_length <= 3; strip_length++) {
			for (unsigned pitch = 0; pitch <= 3; pitch++) {
				std::string output_file = output_dir+"/"+"nbr_strips_"+std::to_string(nbr_of_strips) +"_strip_length_"+
						std::to_string(strip_length)+"_pitch_"+std::to_string(pitch)+".vtk";
				Rect2DDetector rdd(nbr_of_strips, strip_length*100, pitch*100,
						strip_potential, refine_level, max_iter, stop_accuracy, output_file);
				rdd.compute_potential();
			}
		}
	}
}

void test_serrated_2D_potential() {

	ZeroRightHandSide<2> rhs;
	double strip_potential = 1;
	unsigned strip_length = 100;
	unsigned strip_width = 100;
	unsigned pitch = 100;

	/*
	 * SerratedRect2DDetector(unsigned nbr_of_strips,
			unsigned strip_length, unsigned strip_width, unsigned pitch,
			double strip_potential, unsigned refine_level, unsigned max_iter,
			double stop_accuracy, std::string ouput_file)
	 */
	unsigned refine_level = 3;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	std::string output_dir = "tests_output";

	for (unsigned nbr_of_strips = 1; nbr_of_strips <= 10; nbr_of_strips++) {
		std::string output_file = output_dir+"/"+std::to_string(nbr_of_strips)+".vtk";
		SerratedRect2DDetector srdd(nbr_of_strips, strip_length, strip_width,
				pitch, strip_potential, refine_level, max_iter, stop_accuracy, output_file);
		srdd.compute_potential();
	}
}

