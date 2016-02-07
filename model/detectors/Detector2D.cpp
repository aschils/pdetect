/*
 * Detector2D.cpp
 *
 *  Created on: 30 janv. 2016
 *      Author: aschils
 */


#include "Detector2D.hpp"

void Detector2D::compute() {
	potential_solver->compute_solution();

	if(0)
		potential_solver->get_solution(solution_potential);

	solution_potential.sort_cells_by_coord();
	compute_electric_field(solution_potential, electric_field);

	/*
	Point<2> p;
	p[0] = 0;
	p[1] = rect_width_fe-0.9;
	line = new StraightLine<2>(90, p, &solution_potential, 0.1);*/
	//solution_potential.get_values(p);
	//solution_potential.print();
}

void Detector2D::compute_weight(){
	potential_solver_weight->compute_solution();

	if(0)
		potential_solver_weight->get_solution(solution_weight_potential);
	
	solution_weight_potential.sort_cells_by_coord();
	compute_electric_field(solution_weight_potential, electric_field_weight);
}

void Detector2D::compute_electric_field(Solution<2> &potential,
		std::vector<std::pair<typename DoFHandler<2>::active_cell_iterator,
		std::vector<Tensor<1, 2> > > > &electric_field) {

	std::vector<std::pair<typename DoFHandler<2>::active_cell_iterator,
		ValuesAtCell<2> > > values_at_cells =
			potential.values_at_cells;
	electric_field.resize(values_at_cells.size());

	for (unsigned i = 0; i < values_at_cells.size(); i++) {
		std::pair<
			typename DoFHandler<2>::active_cell_iterator,
			std::vector<Tensor<1, 2> >
		> EF_at_one_cell;
		EF_at_one_cell.first = values_at_cells[i].first;
		//Electric field is -grad V
		EF_at_one_cell.second = TensorUtils::opposite_vector_of_tensors<1,2>(
				values_at_cells[i].second.gradient);
		electric_field[i] = EF_at_one_cell;
	}
	//std::cout << electric_field.size() << std::endl;
	//VectorUtils::print_vec_of_pair_of_vec(electric_field);
}

void Detector2D::draw_vtk_graph_potential(
		std::string output_file) {
	solution_potential.draw_vtk_graph_fun(output_file);
}

void Detector2D::draw_vtk_graph_weight_potential(
		std::string output_file){
	solution_weight_potential.draw_vtk_graph_fun(output_file);
}

void Detector2D::draw_vtk_graph_gradient_of_potential(
		std::string output_file) {
	solution_potential.draw_vtk_graph_derivatives(output_file);
}

void Detector2D::draw_vtk_graph_gradient_of_weight_potential(
		std::string output_file) {
	solution_weight_potential.draw_vtk_graph_derivatives(output_file);
}

Detector2D::~Detector2D(){
	delete zero_right_hand_side;
	delete boundary_conditions;
	delete potential_solver;
	delete triangulation;

	//delete line;

	delete boundary_conditions_weight;
	delete potential_solver_weight;
	delete triangulation_weight;
}

