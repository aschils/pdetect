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
	unsigned refine_level = 6;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	std::string output_dir = "tests_electrode_current/";
	Utils::create_directory_if_not_exists(output_dir);

	SerratedRect2DDetector *srdd = new SerratedRect2DDetector(nbr_of_strips,
			width, strip_length, strip_width, half_pitch, strip_potential,
			refine_level, max_iter, stop_accuracy);
	srdd->compute();
	srdd->compute_weight();
	srdd->draw_vtk_graph_potential(output_dir + "electrode_pot.vtk");
	srdd->draw_vtk_graph_weight_potential(
			output_dir + "electrode_pot_weight.vtk");
	srdd->draw_vtk_graph_gradient_of_potential(
			output_dir + "electrode_pot_grad.vtk");

	Solution<2> *solution = new Solution<2>();

	srdd->get_solution(*solution);
	Solution<2> *weight_solution = new Solution<2>();
	srdd->get_solution_weight(*weight_solution);
	MyGeometryInfo *geo_info = srdd->get_geometry_info();

	Point<2> p1(100, 0.0);
	Point<2> p2(100, 100);
	Segment seg(p1, p2);

	Line particle_traj(0, 149);

	ElectrodeCurrent<2> ec(strip_potential, geo_info, solution, weight_solution,
			10);
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

	std::cout << "detector params: " << srdd->params_to_string() << std::endl;

	delete solution;
	delete weight_solution;
	delete srdd;
}

void test_electrode_current_mid_rect_rect() {

	double strip_potential = 90; //typiquement 100V
	unsigned nbr_of_strips = 3;
	unsigned half_width = 150;
	unsigned strip_length = 100;
	unsigned half_strip_width = 10;
	unsigned half_inter_strip_dist = 50;
	unsigned refine_level = 6;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	std::string output_dir = "tests_electrode_current/";
	Utils::create_directory_if_not_exists(output_dir);

	MidRectRect2DDetector *mrr = new MidRectRect2DDetector(half_width,
			strip_length, half_strip_width, half_inter_strip_dist,
			nbr_of_strips, strip_potential, refine_level,
			max_iter, stop_accuracy);
	mrr->compute();
	mrr->compute_weight();
	mrr->draw_vtk_graph_potential(output_dir + "electrode_pot_mid_rect_rect.vtk");
	mrr->draw_vtk_graph_weight_potential(
			output_dir + "electrode_pot_weight_mid_rect_rect.vtk");
	mrr->draw_vtk_graph_gradient_of_potential(
			output_dir + "electrode_pot_grad_mid_rect_rect.vtk");

	Solution<2> *solution = new Solution<2>();

	mrr->get_solution(*solution);
	Solution<2> *weight_solution = new Solution<2>();
	mrr->get_solution_weight(*weight_solution);
	MyGeometryInfo *geo_info = mrr->get_geometry_info();

	//Point<2> p1(0.0,half_width*2);
	//Point<2> p2(nbr_of_strips*(strip_length+2*half_inter_strip_dist), 0);
	Point<2> p1(0,20);
	Point<2> p2(20, 20);
	Line particle_traj(p1,p2);

	ElectrodeCurrent<2> ec(strip_potential, geo_info, solution, weight_solution,
			particle_traj, 13);
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

	delete solution;
	delete weight_solution;
	delete mrr;
}


void gen_comparison_data() {

	unsigned nbr_of_det = 5;
	unsigned max_iter = 10000;
	double stop_accuracy = 10e-12;

	double strip_potential[5] = {90, 90, 90, 500, 10}; //typiquement 100V
	unsigned nbr_of_strips = 1;
	unsigned width[5] = {300, 500, 100, 300, 300};
	unsigned strip_length = 100;
	unsigned strip_width = 2;
	unsigned half_pitch = 50;
	unsigned refine_level = 6;

	for(unsigned i=0; i<nbr_of_det; i++){

		SerratedRect2DDetector srdd(nbr_of_strips, width[i], strip_length,
				strip_width, half_pitch, strip_potential[i],
				refine_level, max_iter, stop_accuracy);
		srdd.compute();
		srdd.compute_weight();
		Solution<2> solution;
		srdd.get_solution(solution);
		Solution<2> weight_solution;
		srdd.get_solution_weight(weight_solution);
		MyGeometryInfo *geo_info = srdd.get_geometry_info();
		ElectrodeCurrent<2> ec(strip_potential[i], geo_info, &solution, &weight_solution,
				10);
		std::vector<std::pair<double, double> > current_vs_time;
		ec.compute_current(current_vs_time);

		std::string output_dir = "../src/plots/"+srdd.params_to_string()+"/";
		Utils::create_directory_if_not_exists(output_dir);

		ResultsOut::write_current_vs_time(output_dir+"pdetect_I_vs_t.txt",
				current_vs_time);
	//std::cout << "detector params: " << srdd->params_to_string() << std::endl;
	}
}

