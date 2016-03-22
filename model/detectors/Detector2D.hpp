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
#include "../boundary_conditions/BoundaryConditions.hpp"
#include "../geometry_info/MyGeometryInfo.hpp"
#include "../Charge.hpp"

class Detector2D {

public:
	void comp_potential();

	void comp_weight_potential();

	void draw_vtk_graph_potential(std::string output_file);

	void draw_vtk_graph_weight_potential(std::string output_file);

	void draw_vtk_graph_gradient_of_potential(std::string output_file);

	void draw_vtk_graph_gradient_of_weight_potential(std::string output_file);

	MyGeometryInfo* get_geometry_info();

	void get_solution(Solution<2> &sol);

	void get_solution_weight(Solution<2> &sol);

	Hole get_hole();

	Electron get_electron();

	double get_hole_pairs_nbr_per_lgth();

	double get_strip_potential();

	double get_first_townsend_coefficient(Point<2> &pos,
			PhysicalValues<2> &values_at_pos);
	
	/*unsigned electric_charge_multiplicator(Point<2> &pos,
			Charge *charge,
			PhysicalValues<2> &values_at_pos,
			double &displacement);*/

	virtual std::string params_to_string() = 0;

	virtual ~Detector2D();

protected:

	const double WEIGHT_STRIP_POTENTIAL = 1.0;

	unsigned max_iter = 1;
	double strip_potential = 1.0;
	double stop_accuracy = 1.0;
	double refine_accuracy;
	unsigned material_id = 0; //Silicium, gaz,...

	Hole hole;
	Electron electron;

	Detector2D(unsigned max_iter, double strip_potential, double stop_accuracy,
			double refine_accuracy, unsigned material_id);

	Triangulation<2> *triangulation = new Triangulation<2>();
	ZeroRightHandSide<2> *zero_right_hand_side = new ZeroRightHandSide<2>();
	BoundaryConditions<2> *boundary_conditions;
	LaplaceSolver<2> *potential_solver;

	Triangulation<2> *triangulation_weight = new Triangulation<2>();
	BoundaryConditions<2> *boundary_conditions_weight;
	LaplaceSolver<2> *potential_solver_weight;

	Solution<2> *solution_potential = new Solution<2>(); 
	Solution<2> *solution_weight_potential = new Solution<2>();

	MyGeometryInfo *geo_info;

	void compute_solution(LaplaceSolver<2> *potential_solver,
			Solution<2> &solution_potential);

	/*unsigned diethorn(Point<2> &pos, PhysicalValues<2> &values_at_pos, double &displacement);

	unsigned townsend_electron_mult(Point<2> &pos, Charge *charge,
			PhysicalValues<2> &values_at_pos, double &displacement);*/

};
