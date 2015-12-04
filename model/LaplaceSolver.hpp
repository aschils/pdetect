/*
 * LaplaceSolver.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#ifndef __LAPLACE_SOLVER_HPP__
#define __LAPLACE_SOLVER_HPP__

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

#include <fstream>
#include <iostream>
#include <unordered_map>

#include "Gradient.hpp"
#include "Solution.hpp"
#include "Utils.hpp"

using namespace dealii;

template<int dim>
class LaplaceSolver {

public:
	LaplaceSolver(Triangulation<dim> *triangulation,
					double rect_length_fe,
					double rect_width_fe,
					unsigned refine_level,
					unsigned max_iter,
					double stop_accuracy,
					const Function<dim> *right_hand_side,
					Function<dim> *boundary_values,
					bool constraints_are_periodic);

	void compute_solution();
	SolutionScalar<dim> get_solution();
	SolutionVector<dim> compute_gradient_of_solution();

private:

	bool constraints_are_periodic;
	unsigned max_iter;
	double stop_accuracy;

	Triangulation<dim> *triangulation;
	FE_Q<dim> fe;
	DoFHandler<dim> dof_handler;
	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;
	Vector<double> solution_vec;
	ConstraintMatrix constraints;

	Vector<double> system_rhs;
	const Function<dim> *right_hand_side;
	Function<dim> *boundary_values_fun;

	void setup_system();
	void assemble_system();
	void solve();
	void output_results(std::string result_file_path) const;
	void make_periodicity_constraints();
};

#include "LaplaceSolver.cpp"

#endif
