/*
 * LaplaceSolver.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template<int dim>
LaplaceSolver<dim>::LaplaceSolver(Triangulation<dim> *triangulation,
		double rect_length_fe, double rect_width_fe, unsigned refine_level,
		unsigned max_iter, double stop_accuracy,
		const Function<dim> *right_hand_side, Function<dim> *boundary_values,
		bool constraints_are_periodic) :
		fe(1), dof_handler(*triangulation) {

	this->constraints_are_periodic = constraints_are_periodic;
	this->triangulation = triangulation;
	this->right_hand_side = right_hand_side;
	this->boundary_values_fun = boundary_values;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;
	this->rect_length_fe = rect_length_fe;
	this->rect_width_fe = rect_width_fe;

	quadrature_formula = new QGauss<dim>(2);;
	fe_values = new FEValues<dim>(fe, *quadrature_formula,
				update_values | update_gradients | update_quadrature_points
						| update_JxW_values | update_second_derivatives);

	if (constraints_are_periodic) {
		for (typename Triangulation<dim>::active_cell_iterator cell =
				triangulation->begin_active(); cell != triangulation->end();
				++cell) {
			for (unsigned int f = 0; f < GeometryInfo < dim > ::faces_per_cell;
					++f) {
				if (cell->face(f)->at_boundary()) {
					if (Utils::equals_double(cell->face(f)->center()[0], 0.0,
							0.000001)
							|| Utils::equals_double(cell->face(f)->center()[0],
									rect_length_fe, 0.000001)
							|| Utils::equals_double(cell->face(f)->center()[1],
									rect_width_fe, 0.000001))

						cell->face(f)->set_boundary_id(1);
				}
			}
		}
	}
	this->triangulation->refine_global(refine_level);
}

template<int dim>
void LaplaceSolver<dim>::make_periodicity_constraints() {
	std::map<unsigned int, double> dof_locations;

	for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
			cell != dof_handler.end(); ++cell) {
		if (cell->at_boundary() && cell->face(1)->at_boundary()) {

			dof_locations[cell->face(1)->vertex_dof_index(0, 0)] =
					cell->face(1)->vertex(0)[1];
			dof_locations[cell->face(1)->vertex_dof_index(1, 0)] =
					cell->face(1)->vertex(1)[1];
		}
	}

	for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
			cell != dof_handler.end(); ++cell) {
		if (cell->at_boundary() && cell->face(0)->at_boundary()) {
			for (unsigned int face_vertex = 0; face_vertex < 2; ++face_vertex) {

				constraints.add_line(
						cell->face(0)->vertex_dof_index(face_vertex, 0));
				std::map<unsigned int, double>::const_iterator p =
						dof_locations.begin();
				for (; p != dof_locations.end(); ++p) {
					if (std::fabs(
							p->second - cell->face(0)->vertex(face_vertex)[1])
							< 1e-8) {
						constraints.add_entry(
								cell->face(0)->vertex_dof_index(face_vertex, 0),
								p->first, 1.0);
						break;
					}
				}
				Assert(p != dof_locations.end(),
						ExcMessage(
								"No corresponding degree of freedom was found!"));
			}
		}
	}
}

template<int dim>
void LaplaceSolver<dim>::setup_system() {

	dof_handler.distribute_dofs(fe);

	if (constraints_are_periodic) {
		ConstraintMatrix constraints;
		constraints.clear();
		make_periodicity_constraints();
		VectorTools::interpolate_boundary_values(dof_handler, 1,
				ZeroFunction<2>(), constraints);
		constraints.close();
		DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
		DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
		dsp.compress();
		sparsity_pattern.copy_from(dsp);
	}

	else {
		DynamicSparsityPattern dsp(dof_handler.n_dofs());
		DoFTools::make_sparsity_pattern(dof_handler, dsp);
		sparsity_pattern.copy_from(dsp);
	}

	system_matrix.reinit(sparsity_pattern);
	solution_vec.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}

template<int dim>
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
			*boundary_values_fun, boundary_values);
	MatrixTools::apply_boundary_values(boundary_values, system_matrix,
			solution_vec, system_rhs);
}

template<int dim>
void LaplaceSolver<dim>::solve() {
	SolverControl solver_control(max_iter, stop_accuracy);
	SolverCG<> solver(solver_control);
	solver.solve(system_matrix, solution_vec, system_rhs,
			PreconditionIdentity());
}

template<int dim>
void LaplaceSolver<dim>::compute_solution() {
	setup_system();
	assemble_system();
	solve();
}

namespace std {

template<>
struct hash<dealii::Point<2> > {
public:
	size_t operator()(Point<2> x) const throw () {
		return hash<double>()(x[0]) ^ hash<double>()(x[1]);
	}
};

}

//template<int dim>
//void LaplaceSolver<dim>::set_solution_at_that_point_as_already_known(
//		std::unordered_map<Point<2>, bool> &already_known, Point<dim> &point) {
//	std::pair<Point<dim>, bool> pair;
//	pair.first = point;
//	pair.second = true;
//	already_known.insert(pair);
//}
//
//template<int dim>
//void LaplaceSolver<dim>::save_solution_at_that_point(Point<dim> &point,
//		double &fun_at_point,
//		Tensor<1, dim> &gradient_at_point,
//		Tensor<2, dim> &hessian_at_point,
//		std::vector<std::pair<std::vector<double>, SolutionData<dim> > >
//		&coord_and_data) {
//
//	std::pair<std::vector<double>, SolutionData<dim> > coord_and_sol_at_one_point;
//	std::vector<double> coord(dim);
//
//	for (unsigned i = 0; i < dim; i++)
//		coord[i] = point[i];
//
//	coord_and_sol_at_one_point.first = coord;
//	SolutionData<dim> sol(fun_at_point, gradient_at_point, hessian_at_point);
//	coord_and_sol_at_one_point.second = sol;
//	coord_and_data.push_back(coord_and_sol_at_one_point);
//}

template<int dim>
void LaplaceSolver<dim>::build_solution(
		std::vector<std::pair<typename DoFHandler<dim>::active_cell_iterator,
		ValuesAtCell<dim> > >
		&values_at_cells) {

	//std::unordered_map<Point<2>, bool> already_known;
	const unsigned int vertices_per_cell = GeometryInfo < dim
			> ::vertices_per_cell;
	std::vector<double> fun_at_pts_of_one_cell(vertices_per_cell);
	std::vector<Tensor<1, dim> > gradient_at_pts_of_one_cell(vertices_per_cell);
	std::vector<Tensor<2, dim> > hessian_at_pts_of_one_cell(vertices_per_cell);

	typename DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();

	for (; cell != endc; cell++) {

		fe_values->reinit(cell);
		fe_values->get_function_values(solution_vec, fun_at_pts_of_one_cell);
		fe_values->get_function_gradients(solution_vec,
				gradient_at_pts_of_one_cell);
		fe_values->get_function_hessians(solution_vec,
				hessian_at_pts_of_one_cell);
		ValuesAtCell<dim> values(fun_at_pts_of_one_cell,
				gradient_at_pts_of_one_cell, hessian_at_pts_of_one_cell);
		std::pair<typename DoFHandler<dim>::active_cell_iterator,
			ValuesAtCell<dim> > values_at_cell;
		values_at_cell.first = cell;
		values_at_cell.second = values;
		values_at_cells.push_back(values_at_cell);
/*
		for (unsigned int v = 0; v < GeometryInfo < 2 > ::vertices_per_cell;
				v++) {

			Point < dim > point = cell->vertex(v);

			bool sol_unknown_at_point = already_known.find(point)
					== already_known.end();
			if (sol_unknown_at_point) {
				save_solution_at_that_point(point, fun_at_pts_of_one_cell[v],
						gradient_at_pts_of_one_cell[v],
						hessian_at_pts_of_one_cell[v],
						coord_and_data);
				set_solution_at_that_point_as_already_known(already_known,
						point);
			}
		}*/
	}
}

template<int dim>
void LaplaceSolver<dim>::get_solution(Solution<dim> &sol) {

	//Used only to output vtk file
//	DataOut<2> fun_drawer;
//	fun_drawer.attach_dof_handler(dof_handler);
//	fun_drawer.add_data_vector(solution_vec, "solution");
//	fun_drawer.build_patches();
//
//	//PostProcessor, used only to output vtk file
//	Derivatives<dim> derivatives;
//	DataOut<dim> derivatives_drawer;
//	derivatives_drawer.attach_dof_handler(dof_handler);
//	derivatives_drawer.add_data_vector(solution_vec, derivatives);
//	derivatives_drawer.build_patches();
//
//	sol.set_fun_drawer(fun_drawer);
//	sol.set_derivatives_drawer(derivatives_drawer);
	build_solution(sol.values_at_cells);

	/*
	 std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	 std::vector<std::pair<std::vector<double>, std::vector<double> > > gradient_at_all_points;
	 std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	 auto duration = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
	 std::cout << "temps: " << duration << std::endl;*/
	//VectorUtils::print_vec_of_pair_of_vec(gradient_at_all_points);
	//std::cout << solution_vec.size() << std::endl;
	//std::cout << gradient_at_all_points.size() << std::endl;
}

/*template<int dim>
ValuesAtCell<dim> LaplaceSolver<dim>::get_solution_at_point(Point<dim> &point,
		std::vector<std::pair<std::vector<double>, ValuesAtCell<dim> > >
		&coord_and_data) {

	//Here we are only checking if the point is in our finite elements domain.
	for(int i = 0; i < dim; i++) {
		if(point[i] > rect_width_fe || point[i] > rect_length_fe) {
			std::cout << "Point is not in the detector" << std::endl;
			return NULL;
		}
	}

	int pos = 0;
	bool notFound = true;
	while(notFound) {
		int coord = 0;
		for(int j = 0; j < dim; j++) {
			if(greater_than_or_equals_double(coord_and_data[pos+1].first[j],
					point[j], 0.000001))
				coord++;
		}
		if(coord == dim)
			notFound = false;
		pos++;
	}
	ValuesAtCell<dim> data_at_point = extrapolate_data_at_point(coord_and_data, pos-1);

	return data_at_point;
}*/


template<int dim>
LaplaceSolver<dim>::~LaplaceSolver() {
	delete fe_values;
	delete quadrature_formula;
}
