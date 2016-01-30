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
#include "../SerratedRect2DBoundaryValues.hpp"
#include "../SerratedRect2DBoundaryValuesWeight.hpp"
#include "../StraightLine.hpp"
#include "../BoundaryConditions.hpp"

class Detector2D {

public:
	void compute();

	void compute_weight();

	void draw_vtk_graph_potential(std::string output_file);

	void draw_vtk_graph_weight_potential(std::string output_file);

	void draw_vtk_graph_gradient_of_potential(
			std::string output_file);

	void draw_vtk_graph_gradient_of_weight_potential(
			std::string output_file);

	virtual std::string params_to_string() = 0;

	virtual ~Detector2D();

protected:
	unsigned refine_level, max_iter = 1;
	double strip_potential = 1.0;
	double stop_accuracy = 1.0;

	Triangulation<2> *triangulation;
	ZeroRightHandSide<2> *zero_right_hand_side;
	//Function<2> *boundary_val;
	BoundaryConditions<2> *boundary_conditions;
	LaplaceSolver<2> *rect_potential_solver;

	Triangulation<2> *triangulation_weight;
	//Function<2> *boundary_val_weight;
	BoundaryConditions<2> *boundary_conditions_weight;
	LaplaceSolver<2> *rect_potential_solver_weight;

	Solution<2> solution_potential, solution_weight_potential;
	std::vector<
			std::pair<typename DoFHandler<2>::active_cell_iterator,
					std::vector<Tensor<1, 2> > > > electric_field;
	std::vector<
			std::pair<typename DoFHandler<2>::active_cell_iterator,
					std::vector<Tensor<1, 2> > > > electric_field_weight;

	StraightLine<2> *line;

	void compute_electric_field(Solution<2> &potential,
				std::vector<std::pair<typename DoFHandler<2>::active_cell_iterator,
				std::vector<Tensor<1, 2> > > > &electric_field);
};
