/*
 * LaplaceSolver.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#pragma once

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/fe_field_function.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_refinement.h>

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <chrono>

#include "Derivatives.hpp"
#include "Solution.hpp"
#include "Utils.hpp"
#include "boundary_conditions/BoundaryConditions.hpp"
#include "Constants.hpp"

using namespace dealii;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


template<unsigned dim>
class LaplaceSolver {

public:
	LaplaceSolver(Triangulation<dim> *triangulation, double refine_accuracy,
			unsigned max_iter, double stop_accuracy,
			const Function<dim> *right_hand_side,
			BoundaryConditions<dim> *boundary_conditions,
			bool constraints_are_periodic);

	void compute_solution();
	void get_solution(Solution<dim> &sol);

	~LaplaceSolver();

private:

	typedef bg::model::point<double, 2, bg::cs::cartesian> bpoint;
	typedef bg::model::box<bpoint> box;
	typedef std::pair<box, std::pair<
			typename DoFHandler<dim>::active_cell_iterator, float>>
	cell_coord_pair;

	bool constraints_are_periodic;
	unsigned max_iter;
	double refine_accuracy, stop_accuracy, max_uncertainty;

	Triangulation<dim> *triangulation;
	FE_Q<dim> fe;
	FEValues<dim> *fe_values;
	QGauss<dim> *quadrature_formula;
	DoFHandler<dim> dof_handler;
	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;
	Vector<double> solution_vec;
	ConstraintMatrix constraints;
	Derivatives<dim> derivatives;

	Vector<double> system_rhs;
	const Function<dim> *right_hand_side;
	BoundaryConditions<dim> *boundary_conditions;

	Vector<float> uncertainty_per_cell;

	void setup_system();
	void assemble_system();
	void solve();
	void compute_uncertainties();
	void refine_grid();
	void output_results(std::string result_file_path) const;
	void build_solution(bgi::rtree<cell_coord_pair, bgi::quadratic<16> > &values_at_cells);
};

template<unsigned dim>
LaplaceSolver<dim>::LaplaceSolver(Triangulation<dim> *triangulation,
		double refine_accuracy, unsigned max_iter, double stop_accuracy,
		const Function<dim> *right_hand_side,
		BoundaryConditions<dim> *boundary_conditions,
		bool constraints_are_periodic) :
		fe(1), dof_handler(*triangulation) {

	this->constraints_are_periodic = constraints_are_periodic;
	this->triangulation = triangulation;
	this->right_hand_side = right_hand_side;
	this->boundary_conditions = boundary_conditions;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;
	this->refine_accuracy = refine_accuracy;

	quadrature_formula = new QGauss<dim>(2);
	fe_values = new FEValues<dim>(fe, *quadrature_formula,
			update_values | update_gradients | update_quadrature_points
					| update_JxW_values | update_second_derivatives);

	if (constraints_are_periodic) {
		for (typename Triangulation<dim>::active_cell_iterator cell =
				triangulation->begin_active(); cell != triangulation->end();
				++cell) {

			boundary_conditions->set_periodicity_constraints(cell);

		}
	}
	this->triangulation->refine_global(1);
}

template<unsigned dim>
void LaplaceSolver<dim>::setup_system() {

	dof_handler.distribute_dofs(fe);

	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);
	solution_vec.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

template<unsigned dim>
void LaplaceSolver<dim>::assemble_system() {

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula->size();
	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);
	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();

	for (; cell != endc; ++cell) {

		fe_values->reinit(cell);
		cell_matrix = 0;
		cell_rhs = 0;

		for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
			for (unsigned int i = 0; i < dofs_per_cell; ++i) {
				for (unsigned int j = 0; j < dofs_per_cell; ++j)

					cell_matrix(i, j) += (fe_values->shape_grad(i, q_index)
							* fe_values->shape_grad(j, q_index)
							* fe_values->JxW(q_index));
				cell_rhs(i) += (fe_values->shape_value(i, q_index)
						* (*right_hand_side).value(
								fe_values->quadrature_point(q_index))
						* fe_values->JxW(q_index));
			}

		cell->get_dof_indices(local_dof_indices);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j) {
				system_matrix.add(local_dof_indices[i], local_dof_indices[j],
						cell_matrix(i, j));
			}
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}

	std::map<types::global_dof_index, double> boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler, 0,
			*(boundary_conditions->get_values()), boundary_values);
	MatrixTools::apply_boundary_values(boundary_values, system_matrix,
			solution_vec, system_rhs);
}

template<unsigned dim>
void LaplaceSolver<dim>::solve() {
	SolverControl solver_control(max_iter, stop_accuracy);
	SolverCG<> solver(solver_control);
	solver.solve(system_matrix, solution_vec, system_rhs,
			PreconditionIdentity());
}

template<unsigned dim>
void LaplaceSolver<dim>::compute_uncertainties() {
	Vector<float> estimated_error_per_cell(triangulation->n_active_cells());
	KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(3),
			typename FunctionMap<dim>::type(), solution_vec,
			estimated_error_per_cell);
	uncertainty_per_cell = estimated_error_per_cell;
	max_uncertainty = uncertainty_per_cell.linfty_norm();
}

template<unsigned dim>
void LaplaceSolver<dim>::refine_grid() {
	GridRefinement::refine_and_coarsen_fixed_number(*triangulation,
			uncertainty_per_cell, 0.3, 0);
	triangulation->execute_coarsening_and_refinement();
}

template<unsigned dim>
void LaplaceSolver<dim>::compute_solution() {
	bool refine = false;
	while (max_uncertainty > refine_accuracy || !refine) {
		if (refine)
			refine_grid();
		setup_system();
		assemble_system();
		solve();
		compute_uncertainties();
		refine = true;
	}
}

template<unsigned dim>
void LaplaceSolver<dim>::build_solution(
		bgi::rtree<cell_coord_pair, bgi::quadratic<16> > &values_at_cells) {

	const unsigned int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;

	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();

	unsigned pos = 0;
	for (; cell != endc; cell++) {

		//Point<dim> bottom_left_dealii = cell->vertex(0);
		//Point<dim> top_right_dealii = cell->vertex(vertices_per_cell - 1);

		//From deal.ii cell to boost geometry polygon
		std::vector<bpoint> boost_points(vertices_per_cell);
		Utils::cell_to_bpoints<dim>(vertices_per_cell, cell, boost_points);


		//bpoint bottom_left(bottom_left_dealii[0], bottom_left_dealii[1]);
		//bpoint top_right(top_right_dealii[0], top_right_dealii[1]);

		//box rect(bottom_left, top_right);

		bg::model::polygon<bpoint> polygon;
		bg::assign_points(polygon, boost_points);

		//std::cout << bg::wkt<bg::model::polygon<bpoint>>(polygon) << std::endl;

		box b = bg::return_envelope<box>(polygon);

		/*bpoint max_corner = b.max_corner();
		bpoint min_corner = b.min_corner();

		std::cout <<"Rectangle: " << std::endl;
		std::cout << "(" << (max_corner.get<0>()) << ", " << (max_corner.get<1>()) << ")" << std::endl;
		std::cout << "(" << (min_corner.get<0>()) << ", " << (min_corner.get<1>()) << ")" << std::endl;*/


		//std::cout << bg::wkt<bg::model::box<bpoint>>(b) << std::endl;


		//boost::geometry::model::polygon<bpoint> polygon;
		//boost::geometry::read_wkt("POLYGON((0 0,1.123 9.987,8.876 2.234,0 0),(3.345 4.456,7.654 8.765,9.123 5.432,3.345 4.456))", polygon);
		//boost::geometry::for_each_point(poly, round_coordinates<bpoint>(0.1));
		//std::cout << "Rounded: " << boost::geometry::wkt(polygon) << std::endl;

		values_at_cells.insert(
				std::make_pair(b,
				//std::make_pair(rect,
				std::make_pair(cell, uncertainty_per_cell[pos])));
		pos++;
	}
}

template<unsigned dim>
void LaplaceSolver<dim>::get_solution(Solution<dim> &sol) {

	//Used only to output vtk file
	DataOut<2> fun_drawer;
	fun_drawer.attach_dof_handler(dof_handler);
	fun_drawer.add_data_vector(solution_vec, "potential");
	fun_drawer.build_patches();

	//PostProcessor, used only to output vtk file
	DataOut<dim> derivatives_drawer;
	derivatives_drawer.attach_dof_handler(dof_handler);
	derivatives_drawer.add_data_vector(solution_vec, derivatives);
	derivatives_drawer.build_patches();

	sol.solution_vec = solution_vec;
	sol.uncertainty_per_cell = uncertainty_per_cell;

	sol.set_fun_drawer(fun_drawer);
	sol.set_derivatives_drawer(derivatives_drawer);

	build_solution(sol.values_at_cells);
}

template<unsigned dim>
LaplaceSolver<dim>::~LaplaceSolver() {
	delete fe_values;
	delete quadrature_formula;
}

