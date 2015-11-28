/*
 * SerratedRect2DDetector.hpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#ifndef __SERRATED_RECT_2D_DETECTOR_HPP__
#define __SERRATED_RECT_2D_DETECTOR_HPP__

#include "LaplaceSolver.hpp"
#include "MyGridGenerator.hpp"
#include "ZeroRightHandSide.hpp"
#include "SerratedRect2DBoundaryValues.hpp"
#include "Detector.hpp"

#define DEFAULT_RECT_WIDTH 300.0 //i.e. in domain language (microm,..)


class SerratedRect2DDetector : public Detector {

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

	void compute_potential(std::string result_file_path);

	std::string params_to_string();

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

	double compute_total_length();
	double compute_rect_width_fe();
	void compute_and_set_fe_values();

};

#endif
