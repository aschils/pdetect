/*
 * Detector2D.hpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

#pragma once

#include <deal.II/numerics/data_out.h>

#include "../LaplaceSolver.hpp"
#include "../MyGridGenerator.hpp"
#include "../ZeroRightHandSide.hpp"
#include "../StraightLine.hpp"
#include "../boundary_conditions/BoundaryConditions.hpp"
#include "../geometry_info/MyGeometryInfo.hpp"

class Detector2D {

public:
	void compute();

	void compute_weight();

	void draw_vtk_graph_potential(std::string output_file);

	void draw_vtk_graph_weight_potential(std::string output_file);

	void draw_vtk_graph_gradient_of_potential(std::string output_file);

	void draw_vtk_graph_gradient_of_weight_potential(std::string output_file);

	MyGeometryInfo* get_geometry_info();

	Solution<2> get_solution();

	Solution<2> get_solution_weight();

	virtual std::string params_to_string() = 0;

	virtual ~Detector2D();

protected:
	unsigned refine_level, max_iter = 1;
	double strip_potential = 1.0;
	double stop_accuracy = 1.0;
	double weight_strip_potential = 1.0;

	Triangulation<2> *triangulation = new Triangulation<2>();
	ZeroRightHandSide<2> *zero_right_hand_side = new ZeroRightHandSide<2>();
	BoundaryConditions<2> *boundary_conditions;
	LaplaceSolver<2> *potential_solver;

	Triangulation<2> *triangulation_weight = new Triangulation<2>();
	BoundaryConditions<2> *boundary_conditions_weight;
	LaplaceSolver<2> *potential_solver_weight;

	Solution<2> solution_potential, solution_weight_potential;

	MyGeometryInfo *geo_info;

	void compute_solution(LaplaceSolver<2> *potential_solver,
			Solution<2> &solution_potential);
};
