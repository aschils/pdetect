#ifndef __PERIODIC_CONSTRAINTS_HPP__
#define __PERIODIC_CONSTRAINTS_HPP__

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>

using namespace dealii;

template<int dim>
class PeriodicConstraints {

public:
	PeriodicConstraints();
	void make_periodicity_constraints(ConstraintMatrix *constraints,
            						DoFHandler<dim> *dof_handler);
};

#include "PeriodicConstraints.cpp"

#endif
