/*
 * DrawGraph.cpp
 *
 *  Created on: 29 nov. 2015
 *      Author: aschils
 */

template<int dim>
EvaluationBase<dim>::~EvaluationBase() {
}

template<int dim>
void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step) {
	refinement_cycle = step;
}

template<int dim>
PointValueEvaluation<dim>::PointValueEvaluation(
		const Point<dim> &evaluation_point, TableHandler &results_table) :
		evaluation_point(evaluation_point), results_table(results_table) {
}

template<int dim>
void PointValueEvaluation<dim>::operator ()(const DoFHandler<dim> &dof_handler,
		const Vector<double> &solution) const {

	double point_value = 1e20;
	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();

	bool evaluation_point_found = false;

	for (; (cell != endc) && !evaluation_point_found; ++cell)
		for (unsigned int vertex = 0;
				vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
			if (cell->vertex(vertex) == evaluation_point) {
				point_value = solution(cell->vertex_dof_index(vertex, 0));
				evaluation_point_found = true;
				break;
			};
	AssertThrow(evaluation_point_found,
			ExcEvaluationPointNotFound(evaluation_point));
	results_table.add_value("DoFs", dof_handler.n_dofs());
	results_table.add_value("u(x_0)", point_value);
}

template<int dim>
SolutionOutput<dim>::SolutionOutput(const std::string &output_name_base,
		const DataOutBase::OutputFormat output_format) :
		output_name_base(output_name_base), output_format(output_format) {
}

template<int dim>
void SolutionOutput<dim>::operator ()(const DoFHandler<dim> &dof_handler,
		const Vector<double> &solution) const {

	DataOut < dim > data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");
	data_out.build_patches();
	std::ostringstream filename;
	filename << output_name_base << "-" << this->refinement_cycle
			<< data_out.default_suffix(output_format) << std::ends;
	std::ofstream out(filename.str().c_str());
	data_out.write(out, output_format);
}

