/*
 * DrawGraph.hpp
 *
 *  Created on: 29 nov. 2015
 *      Author: aschils
 *
 * Code copied form deal.II documentation, tutorial step 13.
 */

#ifndef __GRAPH_DRAWER_HPP__
#define __GRAPH_DRAWER_HPP__

#include <deal.II/numerics/data_out.h>
#include <deal.II/base/table_handler.h>

using namespace dealii;

template<int dim>
class EvaluationBase {
public:
	virtual ~EvaluationBase();
	void set_refinement_cycle(const unsigned int refinement_cycle);
	virtual void operator ()(const DoFHandler<dim> &dof_handler,
			const Vector<double> &solution) const = 0;
protected:
	unsigned int refinement_cycle;
};

template<int dim>
class PointValueEvaluation: public EvaluationBase<dim> {
public:
	PointValueEvaluation(const Point<dim> &evaluation_point,
			TableHandler &results_table);
	virtual void operator ()(const DoFHandler<dim> &dof_handler,
			const Vector<double> &solution) const;

	DeclException1 (ExcEvaluationPointNotFound,
			Point<dim>,
			<< "The evaluation point " << arg1
			<< " was not found among the vertices of the present grid.");
private:
	const Point<dim> evaluation_point;
	TableHandler &results_table;
};

template<int dim>
class SolutionOutput: public EvaluationBase<dim> {

public:
	SolutionOutput(const std::string &output_name_base,
			const DataOutBase::OutputFormat output_format);
	virtual void operator ()(const DoFHandler<dim> &dof_handler,
			const Vector<double> &solution) const;

private:
	const std::string output_name_base;
	const DataOutBase::OutputFormat output_format;
};

#include "EvaluationBase.cpp"

#endif
