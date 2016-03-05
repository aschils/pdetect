/*
 * test_cases_projet.hpp
 *
 *  Created on: 1 mars 2016
 *      Author: aschils
 */

#pragma once

#include "../model/detectors/SerratedRect2DDetector.hpp"
#include "../model/Utils.hpp"

//Semi-conducteur
void test_case_silicium() {

	std::string output_dir = "tests_cases/";
	Utils::create_directory_if_not_exists(output_dir);

	unsigned width = 300; //microm
	unsigned strip_length = 25; //microm
	double potential = 100; //Volt
	unsigned half_pitch = 75 / 2; //microm
	unsigned strip_width = 3;
	unsigned nbr_of_strips = 3;

	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;
	double refine_accuracy = 0.01;
	unsigned ec_refine_level = 7;

	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, potential, refine_accuracy, max_iter, stop_accuracy,
			TYPE_SILICIUM);
	std::cout << "detector params: " << srdd.params_to_string() << std::endl;
	srdd.comp_potential();
	srdd.comp_weight_potential();
	srdd.draw_vtk_graph_potential(output_dir + "silicium_potential.vtk");
	srdd.draw_vtk_graph_weight_potential(
			output_dir + "silicium_weight_potential.vtk");
	srdd.draw_vtk_graph_gradient_of_potential(output_dir + "silicium_grad.vtk");

	std::cout << "Starting to compute current" << std::endl;
	ElectrodeCurrent<2> ec(&srdd, ec_refine_level);

	std::vector<std::pair<double, double> > current_vs_time;
	ec.compute_current(current_vs_time);

	std::string output_graph = output_dir + "silicium_current";
	//gnuplot --persist -e 'set terminal png; plot "current" with dots;' > out.png
	Utils::write_gnu_data_file<2>(output_graph, current_vs_time);
	ResultsOut::write_current_vs_time(
			"../src/plots/silicium_current_vs_time.txt", current_vs_time);
}

void test_case_helium() {

	std::string output_dir = "tests_cases/";
	Utils::create_directory_if_not_exists(output_dir);

	unsigned nbr_of_strips = 3;
	unsigned width = 5000;
	unsigned half_pitch = 5000;
	unsigned strip_length = 9000;
	unsigned strip_width = 0;
	double potential = 10000;

	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;
	double refine_accuracy = 0.01;
	unsigned ec_refine_level = 6;

	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, potential, refine_accuracy, max_iter, stop_accuracy,
			TYPE_HELIUM);
	std::cout << "detector params: " << srdd.params_to_string() << std::endl;
	srdd.comp_potential();
	srdd.comp_weight_potential();
	srdd.draw_vtk_graph_potential(output_dir + "helium_potential.vtk");
	srdd.draw_vtk_graph_weight_potential(
			output_dir + "helium_weight_potential.vtk");
	srdd.draw_vtk_graph_gradient_of_potential(output_dir + "helium_grad.vtk");

	std::cout << "Starting to compute current" << std::endl;
	ElectrodeCurrent<2> ec(&srdd, ec_refine_level);

	std::vector<std::pair<double, double> > current_vs_time;
	ec.compute_current(current_vs_time);

	std::string output_graph = output_dir + "helium_current";
	//gnuplot --persist -e 'set terminal png; plot "current" with dots;' > out.png
	Utils::write_gnu_data_file<2>(output_graph, current_vs_time);
	ResultsOut::write_current_vs_time(
			"../src/plots/helium_current_vs_time.txt", current_vs_time);

}

void test_case_helium_scd() {

	std::string output_dir = "tests_cases/";
	Utils::create_directory_if_not_exists(output_dir);

	unsigned nbr_of_strips = 3;
	unsigned width = 30000;
	unsigned half_pitch = 5000;
	unsigned strip_length = 1000;
	unsigned strip_width = 0;
	double potential = 2000;

	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;
	double refine_accuracy = 0.01;
	unsigned ec_refine_level = 6;

	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, potential, refine_accuracy, max_iter, stop_accuracy,
			TYPE_HELIUM);
	std::cout << "detector params: " << srdd.params_to_string() << std::endl;
	srdd.comp_potential();
	srdd.comp_weight_potential();
	srdd.draw_vtk_graph_potential(output_dir + "helium_potential_scd.vtk");
	srdd.draw_vtk_graph_weight_potential(
			output_dir + "helium_weight_potential_scd.vtk");
	srdd.draw_vtk_graph_gradient_of_potential(output_dir + "helium_grad_scd.vtk");

	std::cout << "Starting to compute current" << std::endl;
	ElectrodeCurrent<2> ec(&srdd, ec_refine_level);

	std::vector<std::pair<double, double> > current_vs_time;
	ec.compute_current(current_vs_time);

	std::string output_graph = output_dir + "helium_current_scd";
	//gnuplot --persist -e 'set terminal png; plot "current" with dots;' > out.png
	Utils::write_gnu_data_file<2>(output_graph, current_vs_time);
	ResultsOut::write_current_vs_time(
			"../src/plots/helium_current_vs_time_scd.txt", current_vs_time);

}


