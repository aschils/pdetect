/*
 * SerratedRect2DDetector.hpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#ifndef __SERRATED_RECT_2D_DETECTOR_HPP__
#define __SERRATED_RECT_2D_DETECTOR_HPP__

#include <unordered_map>

#include "Detector2D.hpp"
#include "LaplaceSolver.hpp"
#include "MyGridGenerator.hpp"
#include "ZeroRightHandSide.hpp"
#include "SerratedRect2DBoundaryValues.hpp"
#include "SerratedRect2DBoundaryValuesWeight.hpp"
#include "VectorUtils.hpp"

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

	SolutionScalar<2> compute_potential();

	SolutionVector<2> compute_gradient_of_potential();

	void compute_electric_field();

	std::vector<double> get_electric_field(Point<2> p);

	SolutionScalar<2> compute_weighting_potential();

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

	unsigned nbr_of_pts_along_x, nbr_of_pts_along_y = 0;

	Triangulation<2> *triangulation;
	ZeroRightHandSide<2> *zero_right_hand_side;
	SerratedRect2DBoundaryValues<2> *boundary_val;
	LaplaceSolver<2> *rect_potential_solver;

	Triangulation<2> *triangulation_weight;
	SerratedRect2DBoundaryValuesWeight<2> *boundary_val_weight;
	LaplaceSolver<2> *rect_potential_solver_weight;

	SolutionScalar<2> potential, weighting_potential;
	SolutionVector<2> gradient_of_potential;
	std::vector<std::pair<std::vector<double>, std::vector<double> > >
	electric_field;

	double compute_total_length();
	double compute_rect_width_fe();
	void compute_and_set_fe_values();
	void nbr_of_points_along_axes();
};

#endif
