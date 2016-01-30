/*
 * SerratedRect2DDetector.hpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#pragma once

#include <unordered_map>

#include "Detector2D.hpp"
#include <functional>
#include "../TensorUtils.hpp"
#include "../SerratedRect2DBoundaryCond.hpp"
#include "../SerratedRect2DBoundaryCondWeight.hpp"

#define DEFAULT_RECT_WIDTH 300.0 //i.e. in domain language (microm,..)

using namespace std;

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

	//~SerratedRect2DDetector();

	//void compute();

	//void compute_weight();

//	void draw_vtk_graph_potential(std::string output_file);
//
//	void draw_vtk_graph_weight_potential(std::string output_file);
//
//	void draw_vtk_graph_gradient_of_potential(std::string output_file);
//
//	void draw_vtk_graph_gradient_of_weight_potential(std::string output_file);

	std::vector<double> get_electric_field(Point<2> p);

	std::string params_to_string();

private:

	unsigned nbr_of_strips, strip_length, strip_width, pitch, total_length = 1;
	double rect_width = 1.0;
	const double DEFAULT_RECT_LENGTH_FE = 10000.0;
	double rect_length_fe = 1.0;
	double rect_width_fe = 1.0;

	double strip_length_fe, strip_width_fe, pitch_length_fe = 1.0;

	unsigned nbr_of_pts_along_x = 0, nbr_of_pts_along_y = 0;

	double compute_total_length();
	double compute_rect_width_fe();
	void compute_and_set_fe_values();
	//void nbr_of_points_along_axes();

	//void compute_electric_field(Solution<2> &potential,
	//				std::vector<std::pair<typename DoFHandler<2>::active_cell_iterator,
	//				std::vector<Tensor<1, 2> > > > &electric_field);
};

