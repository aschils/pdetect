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
#include <chrono>

#include "Derivatives.hpp"
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
	void get_solution(Solution<dim> &sol);

	SolutionData<dim> extrapolate_data_at_point(
					std::vector<std::pair
					<std::vector<double>, SolutionData<dim> > >
					&coord_and_data, double pos);
	SolutionData<dim> get_solution_at_point(Point<dim> &point,
					std::vector<std::pair
					<std::vector<double>, SolutionData<dim> > >
					&coord_and_data);

	~LaplaceSolver();

private:

	bool constraints_are_periodic;
	unsigned max_iter;
	double stop_accuracy;

	double rect_length_fe, rect_width_fe;

	Triangulation<dim> *triangulation;
	FE_Q<dim> fe;
	FEValues <dim> *fe_values;
	QGauss < dim > *quadrature_formula;
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

	void set_solution_at_that_point_as_already_known(
			std::unordered_map<Point<2>, bool> &already_known, Point<dim> &point);
	void save_solution_at_that_point(Point<dim> &point,
	double &fun_at_point,
	Tensor<1, dim> &gradient_at_point,
	Tensor<2, dim> &hessian_at_point,
	std::vector<std::pair<std::vector<double>, SolutionData<dim> > >
	&coord_and_data);
	void build_solution(
			std::vector<std::pair<std::vector<double>, SolutionData<dim> > >
			&coord_and_data);
};

#include "LaplaceSolver.cpp"

#endif
