/*
 * LaplaceSolver.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template<int dim>
LaplaceSolver<dim>::LaplaceSolver(double rect_length_fe, double rect_width_fe,
									const Function<dim> *right_hand_side,
									Function<dim> *boundary_values_fun, 
									PeriodicConstraints<dim> *periodic_constraints,
									std::string result_file_path) :
									fe(1), dof_handler(triangulation) {
	this->right_hand_side = right_hand_side;
	this->periodic_constraints = periodic_constraints;
	this->boundary_values_fun = boundary_values_fun;
	this->result_file_path = result_file_path;
	this->rect_length_fe = rect_length_fe;
	this->rect_width_fe = rect_width_fe;
}

template<int dim>
void LaplaceSolver<dim>::make_grid() {

	Point<dim> point_bot(-this->rect_length_fe/2, -rect_width_fe/2);
	Point<dim> point_top(this->rect_length_fe/2, rect_width_fe/2);

	GridGenerator::hyper_rectangle(triangulation, point_bot, point_top);
	triangulation.begin_active()->face(0)->set_boundary_id(1);
    triangulation.begin_active()->face(1)->set_boundary_id(1);
	triangulation.refine_global(7);
}

template<int dim>
void LaplaceSolver<dim>::setup_system() {
	dof_handler.distribute_dofs(fe);
 	
	constraints.clear();
	periodic_constraints->set_utilities(&constraints, &dof_handler);
    periodic_constraints->make_periodicity_constraints();

    VectorTools::interpolate_boundary_values(dof_handler, 1,
                                             ZeroFunction<2>(),
                                             constraints);
    constraints.close ();


	DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
	dsp.compress();

	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

template<int dim>
void LaplaceSolver<dim>::assemble_system() {

	QGauss<dim> quadrature_formula(2);

	FEValues<dim> fe_values(fe, quadrature_formula,
					update_values | update_gradients | update_quadrature_points
							| update_JxW_values);
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs(dofs_per_cell);
	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();

	for (; cell != endc; ++cell) {

		fe_values.reinit(cell);
		cell_matrix = 0;
		cell_rhs = 0;

		for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
			for (unsigned int i = 0; i < dofs_per_cell; ++i) {
				for (unsigned int j = 0; j < dofs_per_cell; ++j)

					cell_matrix(i, j) += (fe_values.shape_grad(i, q_index)
							* fe_values.shape_grad(j, q_index)
							* fe_values.JxW(q_index));
				cell_rhs(i) += (fe_values.shape_value(i, q_index)
						* (*right_hand_side).value(
								fe_values.quadrature_point(q_index))
						* fe_values.JxW(q_index));
			}

		cell->get_dof_indices(local_dof_indices);

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			for (unsigned int j = 0; j < dofs_per_cell; ++j)
			{
				system_matrix.add(local_dof_indices[i], local_dof_indices[j],
									cell_matrix(i, j));
			}
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}

	std::map<types::global_dof_index, double> boundary_values;
	VectorTools::interpolate_boundary_values(dof_handler, 0,
			*boundary_values_fun, boundary_values);
	MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution,
			system_rhs);
}

template<int dim>
void LaplaceSolver<dim>::solve() {
	SolverControl solver_control(10000, 1e-12);
	SolverCG<> solver(solver_control);
	solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

template<int dim>
void LaplaceSolver<dim>::output_results() const {
	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");
	data_out.build_patches();
	std::ofstream output(result_file_path);
	data_out.write_vtk(output);
}

template<int dim>
void LaplaceSolver<dim>::run() {
	make_grid();
	setup_system();
	assemble_system();
	solve();
	output_results();
}

