/*
 * potential2D.cpp
 *
 *  Created on: 8 nov. 2015
 *      Author: aschils
 */

#include <deal.II/lac/vector.h>

#include "../model/detectors/MidCircleRect2DDetector.hpp"
#include "../model/detectors/SerratedRect2DDetector.hpp"
#include "../model/detectors/MidRectRect2DDetector.hpp"
#include "../model/ZeroRightHandSide.hpp"
#include "../model/Utils.hpp"
#include "../model/Constants.hpp"

void test_serrated_2D_potential() {

	double strip_potential = 200;
	unsigned strip_length = 100;
	unsigned strip_width = 20;
	unsigned half_pitch = 50;
	unsigned max_iter = 100000;
	double stop_accuracy = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	std::string output_dir = "tests_output";
	Utils::create_directory_if_not_exists(output_dir);

	for (unsigned nbr_of_strips = 1; nbr_of_strips <= 10; nbr_of_strips++) {

		std::cout << "Computing potential for " << nbr_of_strips << " strips"
				<< std::endl;

		std::string output_file = output_dir + "/" + "half-pitch_"
				+ std::to_string(half_pitch) + "_"
				+ std::to_string(nbr_of_strips);
		SerratedRect2DDetector srdd(nbr_of_strips, strip_length, strip_width,
				half_pitch, strip_potential, refine_accuracy, max_iter,
				stop_accuracy, TYPE_SILICON);
		srdd.comp_potential();
		srdd.comp_weight_potential();
		srdd.draw_vtk_graph_potential(output_file+ ".vtk");
		srdd.draw_vtk_graph_weight_potential(output_file+ "weight.vtk");
		srdd.draw_vtk_graph_gradient_of_potential(output_file+ "grad.vtk");
		srdd.draw_vtk_graph_gradient_of_weight_potential(output_file+ "grad_weight.vtk");
	}
}

void test_serrated_rect_limit_cases() {

	double strip_potential = 1;
	unsigned width = 300;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	std::string output_dir = "tests_output_limit_cases/";
	Utils::create_directory_if_not_exists(output_dir);

	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 3; nbr_of_strips++) {

		std::cout << "Computing potential for " << nbr_of_strips
				<< " strips in limit case" << std::endl;

		for (unsigned strip_length = 0; strip_length <= 50; strip_length +=
				50) {
			for (unsigned strip_width = 0; strip_width <= 30; strip_width +=
					30) {
				for (unsigned half_pitch = 0; half_pitch <= 25; half_pitch +=
						25) {

					SerratedRect2DDetector srdd(nbr_of_strips, width,
							strip_length, strip_width, half_pitch,
							strip_potential, refine_accuracy, max_iter,
							stop_accuracy, TYPE_SILICON);
					srdd.comp_potential();
					std::string output_file = output_dir
							+ srdd.params_to_string() + ".vtk";
					srdd.draw_vtk_graph_potential(output_file);
				}
			}
		}
	}
}

void test_electric_field() {

	double strip_potential = 1;
	unsigned strip_length = 100;
	unsigned strip_width = 30;
	unsigned half_pitch = 50;
	unsigned max_iter = 100000;
	double stop_accuracy = 10e-12;
	unsigned nbr_of_strips = 4;
	unsigned width = 300;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	std::string output_dir = "tests_output_electric_field/";

	std::cout << "Computing electric field for " << nbr_of_strips << " strips"
			<< std::endl;

	Utils::create_directory_if_not_exists(output_dir);

	std::string output_file = output_dir + "half-pitch_"
			+ std::to_string(half_pitch) + "_" + std::to_string(nbr_of_strips)
			+ ".vtk";
	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, strip_potential, refine_accuracy, max_iter, stop_accuracy);
	srdd.comp_potential();
	srdd.draw_vtk_graph_potential(output_file);
	srdd.draw_vtk_graph_gradient_of_potential(output_dir + "gradient.vtk");
}

void test_weighting_potential() {

	double strip_potential = 100;
	unsigned strip_length = 100;
	unsigned strip_width = 50;
	unsigned half_pitch = 50;
	unsigned max_iter = 100000;
	double stop_accuracy = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.015;

	std::string output_dir = "tests_output_weighting_potential";
	Utils::create_directory_if_not_exists(output_dir);

	for (unsigned nbr_of_strips = 0; nbr_of_strips <= 10; nbr_of_strips++) {

		std::cout << "Computing weighting potential for " << nbr_of_strips
				<< " strips" << std::endl;

		std::string output_file = output_dir + "/" + "half-pitch_"
				+ std::to_string(half_pitch) + "_"
				+ std::to_string(nbr_of_strips) + ".vtk";
		SerratedRect2DDetector srdd(nbr_of_strips, strip_length, strip_width,
				half_pitch, strip_potential, refine_accuracy, max_iter,
				stop_accuracy, TYPE_SILICON);
		srdd.comp_weight_potential();

		srdd.draw_vtk_graph_weight_potential(output_file);
		srdd.draw_vtk_graph_gradient_of_weight_potential(
				output_dir + "gradient.vtk");

	}
}

void test_mid_circle_rect2D_det() {

	std::string output_dir = "tests_mid_circle_rect2D_det/";
	Utils::create_directory_if_not_exists(output_dir);

	unsigned half_width = 100;
	unsigned half_inter_potential_srcs_dist = 50;
	double potential = 10;
	unsigned max_iter = 100000;
	double max_error = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.01;

	for (unsigned nbr_of_potential_src = 1; nbr_of_potential_src <= 3;
			nbr_of_potential_src++) {
		for (unsigned potential_src_radius = 1; potential_src_radius < 40;
				potential_src_radius += 10) {

			MidCircleRect2DDetector det(half_width, nbr_of_potential_src,
					potential_src_radius, half_inter_potential_srcs_dist,
					potential, refine_accuracy, max_iter, max_error,
					TYPE_SILICON);
			std::cout << "number of strips: " << nbr_of_potential_src << std::endl;
			det.comp_potential();
			det.comp_weight_potential();
			det.draw_vtk_graph_potential(
					output_dir + det.params_to_string() + ".vtk");
			det.draw_vtk_graph_gradient_of_potential(
					output_dir + det.params_to_string() + "_grad.vtk");
			det.draw_vtk_graph_weight_potential(
					output_dir + det.params_to_string() + "_weight.vtk");

		}
	}
}

void test_mid_rect_rect_2D_det() {

	unsigned half_width = 100;
	unsigned strip_length = 50;
	unsigned half_strip_width = 25;
	unsigned half_inter_potential_srcs_dist = 50;
	unsigned nbr_of_strips = 3;
	unsigned potential = 10;
	unsigned max_iter = 100000;
	double max_error = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	std::string output_dir = "tests_mid_rect_rect_det/";
	Utils::create_directory_if_not_exists(output_dir);

	MidRectRect2DDetector det(half_width, strip_length, half_strip_width,
			half_inter_potential_srcs_dist, nbr_of_strips, potential,
			refine_accuracy, max_iter, max_error, TYPE_SILICON);

	det.comp_potential();
	det.comp_weight_potential();
	det.draw_vtk_graph_potential(output_dir+"out.vtk");
	det.draw_vtk_graph_weight_potential(output_dir+"out_weight.vtk");
	det.draw_vtk_graph_gradient_of_potential(output_dir+"grad.vtk");
	det.draw_vtk_graph_gradient_of_weight_potential(output_dir+"grad_weight.vtk");
}

void test_various() {
	unsigned nbr_of_strips = 1;
	unsigned strip_length = 50;
	unsigned strip_width = 1;
	unsigned half_pitch = 300;
	double strip_potential = 1;
	unsigned max_iter = 100000;
	double max_error = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	SerratedRect2DDetector srdd(nbr_of_strips, strip_length, strip_width,
			half_pitch, strip_potential, refine_accuracy, max_iter, max_error,
			TYPE_SILICON);
	srdd.comp_potential();

	std::string output_file = "various.vtk";
	srdd.draw_vtk_graph_potential(output_file);

}

