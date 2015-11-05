/*
 * LaplaceSolver.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#ifndef __LAPLACE_SOLVER__HPP__
#define __LAPLACE_SOLVER__HPP__

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

#include "RightHandSide.hpp"
#include "BoundaryValues.hpp"

using namespace dealii;

template<int dim>
class LaplaceSolver {

public:
	LaplaceSolver();
	void run();

private:

	Triangulation<dim> triangulation;
	FE_Q<dim> fe;
	DoFHandler<dim> dof_handler;
	SparsityPattern sparsity_pattern;
	SparseMatrix<double> system_matrix;
	Vector<double> solution;
	Vector<double> system_rhs;

	void make_grid();
	void setup_system();
	void assemble_system();
	void solve();
	void output_results() const;
};

#include "LaplaceSolver.cpp"

#endif
