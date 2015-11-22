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
double SerratedRect2DDetector::compute_total_length(){
	if(nbr_of_strips == 0 || strip_length == 0)
		return RECT_WIDTH;
	return nbr_of_strips * (strip_length + pitch);
}

/**
 * Rectangle width should always be 300microm. Thus finite element rectangle
 * width should be adapted depending on the length of the rectangle in microm.
 */
double SerratedRect2DDetector::compute_rect_width_fe() {
	double one_micron_fe = RECT_LENGTH_FE/total_length;
	return RECT_WIDTH*one_micron_fe;
}

void SerratedRect2DDetector::compute_fe_values(){
	//if(rect_length_fe <= 0)
	//	throw ZERO_OR_NEGATIVE_DOMAIN_LENGTH_ERROR;

	/*
	this->rect_length_fe = rect_length_fe;
	this->strip_potential = strip_potential;

	if(nbr_of_strip == 0 || strip_length == 0){
		strip_length_fe = 0;
		pitch_length_fe = 0;
		strip_pitch_pair_length_fe = 0;
		return;
	}*/

	double total_strips_length = nbr_of_strips*strip_length;
	double total_pitches_length = nbr_of_strips*pitch;
	double total_length = total_strips_length + total_pitches_length;

	double total_strips_length_fe =
			total_strips_length/total_length*RECT_LENGTH_FE;
	double total_pitches_length_fe = RECT_LENGTH_FE - total_strips_length_fe;

	strip_length_fe = total_strips_length_fe/nbr_of_strips;

	this->strip_width_fe = strip_width/strip_length*strip_length_fe;
	pitch_length_fe = total_pitches_length_fe/nbr_of_strips;
}

SerratedRect2DDetector::SerratedRect2DDetector(unsigned nbr_of_strips,
		unsigned strip_length, unsigned strip_width, unsigned pitch,
		double strip_potential, unsigned refine_level, unsigned max_iter,
		double stop_accuracy, std::string ouput_file) {

	this->nbr_of_strips = nbr_of_strips;
	this->strip_length = strip_length;
	this->strip_width = strip_width;
	this->pitch = pitch;
	this->strip_potential = strip_potential;
	total_length = compute_total_length();
	rect_width_fe = compute_rect_width_fe();

	compute_fe_values();

	triangulation = new Triangulation<2>();

	MyGridGenerator<2>::serrated_hyper_rectangle(*triangulation, RECT_LENGTH_FE,
			rect_width_fe, nbr_of_strips, strip_length_fe, strip_width_fe, pitch_length_fe);

	zero_right_hand_side = new ZeroRightHandSide<2>();
	boundary_val = new SerratedRect2DBoundaryValues<2>(nbr_of_strips, strip_length,
			strip_width,
			pitch, RECT_LENGTH_FE, rect_width_fe, strip_potential);

	rect_potential_solver = new LaplaceSolver<2>(triangulation, refine_level,
			max_iter, stop_accuracy,
			zero_right_hand_side, boundary_val,
			ouput_file, true);
}

SerratedRect2DDetector::~SerratedRect2DDetector() {
	delete zero_right_hand_side;
	delete boundary_val;
	delete rect_potential_solver;
	delete triangulation;
}

void SerratedRect2DDetector::compute_potential() {
	rect_potential_solver->run();
}

