/*
 * Detector2D.cpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */

#include "Detector2D.hpp"

void Detector2D::compute() {
	compute_solution(potential_solver, *solution_potential);
}

void Detector2D::compute_weight() {
	compute_solution(potential_solver_weight, *solution_weight_potential);
}

void Detector2D::compute_solution(LaplaceSolver<2> *potential_solver,
		Solution<2> &solution_potential) {
	potential_solver->compute_solution();
	potential_solver->get_solution(solution_potential);
	solution_potential.sort_cells_by_coord();
}

void Detector2D::draw_vtk_graph_potential(std::string output_file) {
	solution_potential->draw_vtk_graph_fun(output_file);
}

void Detector2D::draw_vtk_graph_weight_potential(std::string output_file) {
	solution_weight_potential->draw_vtk_graph_fun(output_file);
}

void Detector2D::draw_vtk_graph_gradient_of_potential(std::string output_file) {
	solution_potential->draw_vtk_graph_derivatives(output_file);
}

void Detector2D::draw_vtk_graph_gradient_of_weight_potential(
		std::string output_file) {
	solution_weight_potential->draw_vtk_graph_derivatives(output_file);
}

MyGeometryInfo* Detector2D::get_geometry_info() {
	return geo_info;
}

void Detector2D::get_solution(Solution<2> &sol) {
	sol = *solution_potential;
}

void Detector2D::get_solution_weight(Solution<2> &sol) {
	sol = *solution_weight_potential;
}

Detector2D::~Detector2D() {
	delete zero_right_hand_side;
	delete boundary_conditions;
	delete solution_potential;
	delete potential_solver;
	delete triangulation;

	delete boundary_conditions_weight;
	delete solution_weight_potential;
	delete potential_solver_weight;
	delete triangulation_weight;

	delete geo_info;
}

