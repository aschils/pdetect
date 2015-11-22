/*
 * LaplaceSolver.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template<int dim>
LaplaceSolver<dim>::LaplaceSolver(Triangulation<dim> *triangulation,
		unsigned refine_level, unsigned max_iter, double stop_accuracy,
		const Function<dim> *right_hand_side, Function<dim> *boundary_values,
		std::string result_file_path, bool constraints_are_periodic) :
		fe(1), dof_handler(*triangulation) {

	this->constraints_are_periodic = constraints_are_periodic;
	this->triangulation = triangulation;
	this->right_hand_side = right_hand_side;
	this->boundary_values_fun = boundary_values;
	this->result_file_path = result_file_path;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;

	if(constraints_are_periodic){
		this->triangulation->begin_active()->face(0)->set_boundary_id(1);
		this->triangulation->begin_active()->face(1)->set_boundary_id(1);
	}

	this->triangulation->refine_global(refine_level);
}

template<int dim>
void LaplaceSolver<dim>::make_periodicity_constraints(
		ConstraintMatrix *constraints, DoFHandler<dim> *dof_handler) {
    std::map<unsigned int, double> dof_locations;

    for(DoFHandler<2>::active_cell_iterator cell = dof_handler->begin_active();
    			cell != dof_handler->end(); ++cell) {
		if(cell->at_boundary() && cell->face(1)->at_boundary()) {

			dof_locations[cell->face(1)->vertex_dof_index(0, 0)]
            				= cell->face(1)->vertex(0)[1];
          	dof_locations[cell->face(1)->vertex_dof_index(1, 0)]
            				= cell->face(1)->vertex(1)[1];
        }
    }

    for (DoFHandler<2>::active_cell_iterator cell = dof_handler->begin_active ();
         			cell != dof_handler->end (); ++cell) {
		if (cell->at_boundary() && cell->face(0)->at_boundary()) {
          	for (unsigned int face_vertex = 0; face_vertex<2; ++face_vertex) {

				constraints->add_line(cell->face(0)->vertex_dof_index(face_vertex, 0));
				std::map<unsigned int, double>::const_iterator p = dof_locations.begin();
              	for(; p != dof_locations.end(); ++p) {
                	if(std::fabs(p->second-cell->face(0)->vertex(face_vertex)[1])
                					< 1e-8) {
						constraints->add_entry(cell->face(0)->vertex_dof_index(face_vertex, 0),
                                    p->first, 1.0);
                    	break;
                  	}
                }
              	Assert (p != dof_locations.end(),
                    ExcMessage ("No corresponding degree of freedom was found!"));
            }
        }
    }
}

template<int dim>
void LaplaceSolver<dim>::setup_system() {

	dof_handler.distribute_dofs(fe);

	if(constraints_are_periodic){
		ConstraintMatrix constraints;
		std::cout << "before clear" << std::endl;
		constraints.clear();
		std::cout << "before period constr" << std::endl;
		make_periodicity_constraints(&constraints, &dof_handler);
		std::cout << "before interpolate bv" << std::endl;
		VectorTools::interpolate_boundary_values(dof_handler, 1, ZeroFunction<2>(),
				constraints);
		std::cout << "before const close" << std::endl;
		constraints.close();
		std::cout << "before dsp" << std::endl;
		DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
		std::cout << "before make_sparsit" << std::endl;
		DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
		dsp.compress();
		sparsity_pattern.copy_from(dsp);
	}

	else{
		DynamicSparsityPattern dsp(dof_handler.n_dofs());
		DoFTools::make_sparsity_pattern(dof_handler, dsp);
		sparsity_pattern.copy_from(dsp);
	}

	system_matrix.reinit(sparsity_pattern);
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

template<int dim>
void LaplaceSolver<dim>::assemble_system() {

	QGauss < dim > quadrature_formula(2);

	FEValues < dim
			> fe_values(fe, quadrature_formula,
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
			for (unsigned int j = 0; j < dofs_per_cell; ++j) {
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
	SolverControl solver_control(max_iter, stop_accuracy); //10000 1e-12
	SolverCG<> solver(solver_control);
	solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

template<int dim>
void LaplaceSolver<dim>::output_results() const {
	DataOut < dim > data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");
	data_out.build_patches();
	std::ofstream output(result_file_path);
	data_out.write_vtk(output);
}

template<int dim>
void LaplaceSolver<dim>::run() {
	std::cout << "before setup" << std::endl;
	setup_system();
	std::cout << "after setup" << std::endl;
	assemble_system();
	std::cout << "assemble system" << std::endl;
	solve();
	std::cout << "after sovle" << std::endl;
	output_results();
}

