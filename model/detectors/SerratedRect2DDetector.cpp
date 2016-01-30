/*
 * SerratedRect2DDetector.cpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#include "SerratedRect2DDetector.hpp"

/**
 * Compute the length of the rectangle (detector) in the domain language unit
 * (microm) depending on the number of strips, the strip length and the pitch.
 */
double SerratedRect2DDetector::compute_total_length() {
	if (nbr_of_strips == 0)
		return pitch;
	else
		return nbr_of_strips * (strip_length + pitch);
}

/**
 * Rectangle width should always be 300microm. Thus finite element rectangle
 * width should be adapted depending on the length of the rectangle in microm.
 */
double SerratedRect2DDetector::compute_rect_width_fe() {
	double one_micron_fe =
			(total_length == 0) ? 1 : rect_length_fe / total_length;
	return rect_width * one_micron_fe;
}

/**
 * Translate lengths expressed in the domain langage (i.e. microm,...) in
 * length in the finite elements domain.
 */
void SerratedRect2DDetector::compute_and_set_fe_values() {

	rect_length_fe = (total_length == 0) ? 0 : DEFAULT_RECT_LENGTH_FE;
	unsigned total_strips_length = nbr_of_strips * strip_length;
	double total_strips_length_fe =
			(total_length == 0) ?
					0 :
					total_strips_length / (double) total_length
							* rect_length_fe;
	double total_pitches_length_fe = rect_length_fe - total_strips_length_fe;

	strip_length_fe =
			(nbr_of_strips == 0) ? 0 : total_strips_length_fe / nbr_of_strips;
	strip_width_fe =
			(strip_length == 0) ?
					0 : strip_width / (double) strip_length * strip_length_fe;
	pitch_length_fe =
			(nbr_of_strips == 0) ?
					rect_length_fe : total_pitches_length_fe / nbr_of_strips;
	rect_width_fe = compute_rect_width_fe();
}

SerratedRect2DDetector::SerratedRect2DDetector(unsigned nbr_of_strips,
		unsigned strip_length, unsigned strip_width, unsigned pitch,
		double strip_potential, unsigned refine_level, unsigned max_iter,
		double stop_accuracy) :
		SerratedRect2DDetector(nbr_of_strips, DEFAULT_RECT_WIDTH, strip_length,
				strip_width, pitch, strip_potential, refine_level, max_iter,
				stop_accuracy) {
}

SerratedRect2DDetector::SerratedRect2DDetector(unsigned nbr_of_strips,
		unsigned width, unsigned strip_length, unsigned strip_width,
		unsigned pitch, double strip_potential, unsigned refine_level,
		unsigned max_iter, double stop_accuracy) {

	this->rect_width = width;
	this->nbr_of_strips = nbr_of_strips;
	this->strip_length = strip_length;
	this->strip_width = strip_width;
	this->pitch = pitch;
	this->strip_potential = strip_potential;
	this->refine_level = refine_level;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;
	total_length = compute_total_length();
	compute_and_set_fe_values();

	triangulation = new Triangulation<2>();

	MyGridGenerator<2>::serrated_hyper_rectangle(*triangulation, rect_width_fe,
			nbr_of_strips, strip_length_fe, strip_width_fe, pitch_length_fe);
	zero_right_hand_side = new ZeroRightHandSide<2>();
	/*boundary_val = new SerratedRect2DBoundaryValues<2>(nbr_of_strips,
			rect_length_fe, rect_width_fe, strip_potential, pitch_length_fe,
			strip_length_fe, strip_width_fe);*/
	boundary_conditions = new SerratedRect2DBoundaryCond<2>(nbr_of_strips,
			rect_length_fe, rect_width_fe, strip_potential, pitch_length_fe,
			strip_length_fe, strip_width_fe);


	rect_potential_solver = new LaplaceSolver<2>(triangulation, refine_level,
			max_iter, stop_accuracy, zero_right_hand_side, boundary_conditions,
			true);

	boundary_conditions_weight = new SerratedRect2DBoundaryCondWeight<2>(
			nbr_of_strips, rect_length_fe, rect_width_fe, strip_potential,
			pitch_length_fe, strip_length_fe, strip_width_fe);
	triangulation_weight = new Triangulation<2>();
	MyGridGenerator<2>::serrated_hyper_rectangle(*triangulation_weight,
			rect_width_fe, nbr_of_strips, strip_length_fe, strip_width_fe,
			pitch_length_fe);
	rect_potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
			refine_level, max_iter, stop_accuracy, zero_right_hand_side,
			boundary_conditions_weight, true);
	//nbr_of_points_along_axes();
}

/**
 * The number of points along an axis is simply the number
 * of cells along this axes plus one. 
 * (The bottom left corner of each cells plus de bottom right 
 * corner of the last cell)
 */
/*void SerratedRect2DDetector::nbr_of_points_along_axes() {

	Triangulation<2>::active_cell_iterator cell = triangulation->begin_active(),
			endc = triangulation->end();
	for(; cell != endc; ++cell) {
		for(unsigned int v = 0; v < GeometryInfo<2>::faces_per_cell; ++v) {
			if(cell->face(v)->at_boundary()){
				Point<2> p = cell->face(v)->center();

				if(Utils::equals_double(p[1], 0, 0.000001))
					nbr_of_pts_along_x++;
				else if(Utils::equals_double(p[0], 0, 0.000001))
					nbr_of_pts_along_y++;
			}
		}
	}

	nbr_of_pts_along_x++;
	nbr_of_pts_along_y++;
}*/

/*SerratedRect2DDetector::~SerratedRect2DDetector() {
	delete zero_right_hand_side;
	delete boundary_val;
	delete rect_potential_solver;
	delete triangulation;

	//delete line;

	delete boundary_val_weight;
	delete rect_potential_solver_weight;
	delete triangulation_weight;
}*/

//void SerratedRect2DDetector::compute() {
//	rect_potential_solver->compute_solution();
//	rect_potential_solver->get_solution(solution_potential);
//	solution_potential.sort_cells_by_coord();
//	compute_electric_field(solution_potential, electric_field);
//
//	/*
//	Point<2> p;
//	p[0] = 0;
//	p[1] = rect_width_fe-0.9;
//	line = new StraightLine<2>(90, p, &solution_potential, 0.1);*/
//	//solution_potential.get_values(p);
//	//solution_potential.print();
//}

//void SerratedRect2DDetector::compute_weight(){
//	rect_potential_solver_weight->compute_solution();
//	rect_potential_solver_weight->get_solution(solution_weight_potential);
//	solution_weight_potential.sort_cells_by_coord();
//	compute_electric_field(solution_weight_potential, electric_field_weight);
//}

//void SerratedRect2DDetector::compute_electric_field(Solution<2> &potential,
//		std::vector<std::pair<typename DoFHandler<2>::active_cell_iterator,
//		std::vector<Tensor<1, 2> > > > &electric_field) {
//
//	std::vector<std::pair<typename DoFHandler<2>::active_cell_iterator,
//		ValuesAtCell<2> > > values_at_cells =
//			potential.values_at_cells;
//	electric_field.resize(values_at_cells.size());
//
//	for (unsigned i = 0; i < values_at_cells.size(); i++) {
//		std::pair<
//			typename DoFHandler<2>::active_cell_iterator,
//			std::vector<Tensor<1, 2> >
//		> EF_at_one_cell;
//		EF_at_one_cell.first = values_at_cells[i].first;
//		//Electric field is -grad V
//		EF_at_one_cell.second = TensorUtils::opposite_vector_of_tensors<1,2>(
//				values_at_cells[i].second.gradient);
//		electric_field[i] = EF_at_one_cell;
//	}
//	//std::cout << electric_field.size() << std::endl;
//	//VectorUtils::print_vec_of_pair_of_vec(electric_field);
//}

/*std::vector<double> SerratedRect2DDetector::get_electric_field(Point<2> p) {

}*/

//void SerratedRect2DDetector::draw_vtk_graph_potential(
//		std::string output_file) {
//	solution_potential.draw_vtk_graph_fun(output_file);
//}
//
//void SerratedRect2DDetector::draw_vtk_graph_weight_potential(
//		std::string output_file){
//	solution_weight_potential.draw_vtk_graph_fun(output_file);
//}
//
//void SerratedRect2DDetector::draw_vtk_graph_gradient_of_potential(
//		std::string output_file) {
//	solution_potential.draw_vtk_graph_derivatives(output_file);
//}
//
//void SerratedRect2DDetector::draw_vtk_graph_gradient_of_weight_potential(
//		std::string output_file) {
//	solution_weight_potential.draw_vtk_graph_derivatives(output_file);
//}

std::string SerratedRect2DDetector::params_to_string() {

	std::string str = "width" + std::to_string(rect_width) + "_nbr_of_strips_"
			+ std::to_string(nbr_of_strips) + +"_strip_length_"
			+ std::to_string(strip_length) + "_strip_width_"
			+ std::to_string(strip_width) + "_pitch_" + std::to_string(pitch)
			+ "_strip_potential_" + std::to_string(strip_potential)
			+ "_refine_level_" + std::to_string(refine_level) + "_max_iter_"
			+ std::to_string(max_iter) + "_stop_accuracy_"
			+ std::to_string(stop_accuracy);
	return str;
}
