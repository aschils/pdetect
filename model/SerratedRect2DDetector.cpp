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
	boundary_val = new SerratedRect2DBoundaryValues<2>(nbr_of_strips,
			rect_length_fe, rect_width_fe, strip_potential, pitch_length_fe,
			strip_length_fe, strip_width_fe);
	rect_potential_solver = new LaplaceSolver<2>(triangulation, rect_length_fe, 
			rect_width_fe, refine_level, max_iter, stop_accuracy, 
			zero_right_hand_side, boundary_val, true);

	boundary_val_weight = new SerratedRect2DBoundaryValuesWeight<2>(
			nbr_of_strips, rect_length_fe, rect_width_fe, strip_potential,
			pitch_length_fe, strip_length_fe, strip_width_fe);
	triangulation_weight = new Triangulation<2>();
	MyGridGenerator<2>::serrated_hyper_rectangle(*triangulation_weight,
			rect_width_fe, nbr_of_strips, strip_length_fe, strip_width_fe,
			pitch_length_fe);
	rect_potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
			rect_length_fe, rect_width_fe, refine_level, max_iter, stop_accuracy, 
			zero_right_hand_side, boundary_val_weight, true);
}

SerratedRect2DDetector::~SerratedRect2DDetector() {
	delete zero_right_hand_side;
	delete boundary_val;
	delete rect_potential_solver;
	delete triangulation;

	delete boundary_val_weight;
	delete rect_potential_solver_weight;
	delete triangulation_weight;
}

SolutionScalar<2> SerratedRect2DDetector::compute_potential() {
	rect_potential_solver->compute_solution();
	return rect_potential_solver->get_solution();
}

SolutionVector<2> SerratedRect2DDetector::compute_electric_field() {
	return rect_potential_solver->compute_gradient_of_solution();
}

SolutionScalar<2> SerratedRect2DDetector::compute_weighting_potential() {
	rect_potential_solver_weight->compute_solution();
	return rect_potential_solver_weight->get_solution();
}

std::string SerratedRect2DDetector::params_to_string() {

	std::string str = "width"+ std::to_string(rect_width) +"_nbr_of_strips_" +
			std::to_string(nbr_of_strips)
			+ +"_strip_length_" + std::to_string(strip_length) + "_strip_width_"
			+ std::to_string(strip_width) + "_pitch_" + std::to_string(pitch)
			+ "_strip_potential_" + std::to_string(strip_potential)
			+ "_refine_level_" + std::to_string(refine_level) + "_max_iter_"
			+ std::to_string(max_iter) + "_stop_accuracy_"
			+ std::to_string(stop_accuracy);
	return str;
}

