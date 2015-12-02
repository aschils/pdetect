/*
 * SerratedRect2DDetector.hpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#ifndef __SERRATED_RECT_2D_DETECTOR_HPP__
#define __SERRATED_RECT_2D_DETECTOR_HPP__

#include "Detector2D.hpp"
#include "LaplaceSolver.hpp"
#include "MyGridGenerator.hpp"
#include "ZeroRightHandSide.hpp"
#include "SerratedRect2DBoundaryValues.hpp"
#include "SerratedRect2DBoundaryValuesWeight.hpp"

#define DEFAULT_RECT_WIDTH 300.0 //i.e. in domain language (microm,..)


class SerratedRect2DDetector : public Detector2D {

public:
	SerratedRect2DDetector(unsigned nbr_of_strips,
			unsigned strip_length, unsigned strip_width, unsigned pitch,
			double strip_potential, unsigned refine_level, unsigned max_iter,
			double stop_accuracy);

	SerratedRect2DDetector(unsigned nbr_of_strips, unsigned width,
				unsigned strip_length, unsigned strip_width, unsigned pitch,
				double strip_potential, unsigned refine_level, unsigned max_iter,
				double stop_accuracy);

	~SerratedRect2DDetector();

	Solution<2> compute_potential();

	DataOut<2> compute_electric_field();

	Solution<2> compute_weighting_potential();

	std::string params_to_string();

	std::pair<double, double> get_electric_field(Point<2> p);

private:

	unsigned nbr_of_strips, strip_length, strip_width, pitch, total_length = 1;
	unsigned refine_level, max_iter = 1;
	double strip_potential = 1.0;
	double rect_width = 1.0;
	const double DEFAULT_RECT_LENGTH_FE = 10000.0;
	double rect_length_fe = 1.0;
	double rect_width_fe = 1.0;
	double stop_accuracy = 1.0;

	double strip_length_fe, strip_width_fe, pitch_length_fe = 1.0;

	Triangulation<2> *triangulation;
	ZeroRightHandSide<2> *zero_right_hand_side;
	SerratedRect2DBoundaryValues<2> *boundary_val;
	LaplaceSolver<2> *rect_potential_solver;

	Triangulation<2> *triangulation_weight;
	SerratedRect2DBoundaryValuesWeight<2> *boundary_val_weight;
	LaplaceSolver<2> *rect_potential_solver_weight;

	DataOut<2> electric_field_data_container;

	std::unordered_map<std::pair<double, double>, std::pair<double, double> >
			 electric_field;

	double compute_total_length();
	double compute_rect_width_fe();
	void compute_and_set_fe_values();

};

#endif
