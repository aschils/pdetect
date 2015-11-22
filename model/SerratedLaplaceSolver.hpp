/*
 * SerratedLaplaceSolver.hpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#ifndef __SERRATED_LAPLACE_SOLVER_HPP__
#define __SERRATED_LAPLACE_SOLVER_HPP__

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <deal.II/base/logstream.h>

#include "MyGridGenerator.hpp"
#include "Rect2DBoundaryValues.hpp"

using namespace dealii;

template<int dim>
class SerratedLaplaceSolver {

public:
	SerratedLaplaceSolver(double rect_length_fe,
			double rect_width_fe, unsigned nbr_of_strips, unsigned strip_length,
			unsigned strip_width, unsigned pitch,
			const Function<dim> *right_hand_side,
			Function<dim> *boundary_values_fun, std::string result_file_path);
	void run();

private:

	double rect_length_fe, rect_width_fe = 1.0;
	unsigned nbr_of_strips, strip_length, strip_width, pitch = 1;

	Triangulation<dim> triangulation;
	FE_Q<dim> fe;
	DoFHandler<dim> dof_handler;
	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;
	Vector<double> solution;
	Vector<double> system_rhs;
	const Function<dim> *right_hand_side;
	Function<dim> *boundary_values_fun;
	std::string result_file_path;

	void make_grid();
	void setup_system();
	void assemble_system();
	void solve();
	void output_results() const;
};

#include "SerratedLaplaceSolver.cpp"

#endif
