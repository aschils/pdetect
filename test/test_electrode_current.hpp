/*
 * test_electrode_current.hpp
 *
 *  Created on: 11 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "../model/ElectrodeCurrent.hpp"
#include "../model/ResultsOut.hpp"

void test_electrode_current_serrated() {

	double strip_potential = 90; //typiquement 100V
	unsigned nbr_of_strips = 1;
	unsigned width = 300;
	unsigned strip_length = 100;
	unsigned strip_width = 2;
	unsigned half_pitch = 50;
	unsigned max_iter = 100000;
	double stop_accuracy = 10e-12;
	//Use to define the number of charges along the trajectory => 2^refine_level
	unsigned refine_level = 9;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	std::string output_dir = "tests_electrode_current/";
	Utils::create_directory_if_not_exists(output_dir);

	SerratedRect2DDetector srdd(nbr_of_strips, width, strip_length, strip_width,
			half_pitch, strip_potential, refine_accuracy, max_iter,
			stop_accuracy, TYPE_SILICON);
	srdd.comp_potential();
	std::cout << "after comp pot" << std::endl;
	srdd.comp_weight_potential();
	std::cout << "after comp weight pot" << std::endl;
	srdd.draw_vtk_graph_potential(output_dir + "electrode_pot.vtk");
	srdd.draw_vtk_graph_weight_potential(
			output_dir + "electrode_pot_weight.vtk");
	srdd.draw_vtk_graph_gradient_of_potential(
			output_dir + "electrode_pot_grad.vtk");

	double x = half_pitch + (double) strip_length / 2.0;

	Line particle_traj(0, 149);

	ElectrodeCurrent<2> ec(&srdd, refine_level);
	//ec.print_charges();
	//double delta_t = 0.0000000000001; //100ps p.76, V_b = 100V, v_d = 30V
	std::vector<std::pair<double, double> > current_vs_time;
	ec.compute_current(current_vs_time);

	//std::cout << current_vs_time[0].second << " " << std::endl;

	std::string output_graph = output_dir + "current";
	//gnuplot --persist -e 'set terminal png; plot "current" with dots;' > out.png
	Utils::write_gnu_data_file<2>(output_graph, current_vs_time);
	ResultsOut::write_current_vs_time("../src/plots/pdetect_I_vs_t_5_strip.txt",
			current_vs_time);
}

void test_electrode_current_mid_rect_rect() {

	double strip_potential = 90; //typiquement 100V
	unsigned nbr_of_strips = 3;
	unsigned half_width = 150;
	unsigned strip_length = 100;
	unsigned half_strip_width = 10;
	unsigned half_inter_strip_dist = 50;
	unsigned max_iter = 100000;
	double stop_accuracy = 10e-12;
	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	std::string output_dir = "tests_electrode_current/";
	Utils::create_directory_if_not_exists(output_dir);

	MidRectRect2DDetector *mrr = new MidRectRect2DDetector(half_width,
			strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips, strip_potential, refine_accuracy, max_iter,
			stop_accuracy, TYPE_SILICON);
	mrr->comp_potential();
	mrr->comp_weight_potential();
	mrr->draw_vtk_graph_potential(
			output_dir + "electrode_pot_mid_rect_rect.vtk");
	mrr->draw_vtk_graph_weight_potential(
			output_dir + "electrode_pot_weight_mid_rect_rect.vtk");
	mrr->draw_vtk_graph_gradient_of_potential(
			output_dir + "electrode_pot_grad_mid_rect_rect.vtk");

	//Point<2> p1(0.0,half_width*2);
	//Point<2> p2(nbr_of_strips*(strip_length+2*half_inter_strip_dist), 0);
	/*Point<2> p1(0, 20);
	Point<2> p2(0, 40);
	Line particle_traj(p1, p2);*/

	ElectrodeCurrent<2> ec(mrr, 10);
	//ec.print_charges();
	//double delta_t = 0.0000000000001; //100ps p.76, V_b = 100V, v_d = 30V
	std::vector<std::pair<double, double> > current_vs_time;
	ec.compute_current(current_vs_time);

	//std::cout << current_vs_time[0].second << " " << std::endl;

	std::string output_graph = output_dir + "current_mrr";
	//gnuplot --persist -e 'set terminal png; plot "current" with dots;' > out.png
	Utils::write_gnu_data_file<2>(output_graph, current_vs_time);
	//ResultsOut::write_current_vs_time("../src/plots/pdetect_I_vs_t_mrr.txt",
	//		current_vs_time);

	std::cout << "detector params: " << mrr->params_to_string() << std::endl;

	delete mrr;
}

void gen_comparison_data() {

	unsigned nbr_of_det = 5;
	unsigned max_iter = 100000;
	double stop_accuracy = 10e-12;

	double strip_potential[5] = { 90, 90, 90, 500, 10 }; //typiquement 100V
	unsigned nbr_of_strips = 1;
	unsigned width[5] = { 300, 500, 100, 300, 300 };
	unsigned strip_length = 100;
	unsigned strip_width = 2;
	unsigned half_pitch = 50;

	/*
	 * refine_accuracy suggested value between:
	 * 				0.009 => Quick but not precise
	 *				0.008 ; 0.005 ; 0.0045 ; 0.003
	 *				0.0025 => Slow but very precise
	 */
	double refine_accuracy = 0.009;

	for (unsigned i = 0; i < nbr_of_det; i++) {

		SerratedRect2DDetector srdd(nbr_of_strips, width[i], strip_length,
				strip_width, half_pitch, strip_potential[i], refine_accuracy,
				max_iter, stop_accuracy);
		srdd.comp_potential();
		srdd.comp_weight_potential();
		ElectrodeCurrent<2> ec(&srdd, 10);
		std::vector<std::pair<double, double> > current_vs_time;
		ec.compute_current(current_vs_time);

		std::string output_dir = "../src/plots/" + srdd.params_to_string()
				+ "/";
		Utils::create_directory_if_not_exists(output_dir);

		ResultsOut::write_current_vs_time(output_dir + "pdetect_I_vs_t.txt",
				current_vs_time);
		//std::cout << "detector params: " << srdd->params_to_string() << std::endl;
	}
}

