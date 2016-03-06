/*
 * Detector2D.cpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */

#include "Detector2D.hpp"

Detector2D::Detector2D(unsigned max_iter, double strip_potential,
		double stop_accuracy, double refine_accuracy, unsigned material_id) :
		hole(material_id), electron(material_id) {
	this->max_iter = max_iter;
	this->strip_potential = strip_potential;
	this->stop_accuracy = stop_accuracy;
	this->refine_accuracy = refine_accuracy;
	this->material_id = material_id;
}

void Detector2D::comp_potential() {
	compute_solution(potential_solver, *solution_potential);
}

void Detector2D::comp_weight_potential() {
	compute_solution(potential_solver_weight, *solution_weight_potential);
}

void Detector2D::compute_solution(LaplaceSolver<2> *potential_solver,
		Solution<2> &solution_potential) {
	potential_solver->compute_solution();
	potential_solver->get_solution(solution_potential);
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

Hole Detector2D::get_hole(){
	return this->hole;
}

Electron Detector2D::get_electron(){
	return this->electron;
}

double Detector2D::get_strip_potential(){
	return strip_potential;
}

/**
 * Return the number of hole-electron pairs per microm created by
 * the traversal of a particle.
 */
double Detector2D::get_hole_pairs_nbr_per_lgth(){

	switch(material_id){
	{case TYPE_SILICIUM:
		return 75;}
	{case TYPE_HELIUM:
		double energy_per_surface = 1.60217656535e-4; //microJ
		double temperature = 300; //°K
		return energy_per_surface*ATMOSPHERIC_PRESSURE*MOLAR_MASS_HELIUM/
				(GAS_CONSTANT*temperature);}
	default:
		return 75;
	}
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
