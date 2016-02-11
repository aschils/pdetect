/*
 * SerratedRect2DDetector.cpp
 *
 *  Created on: 12 nov. 2015
 *      Author: aschils
 */

#include "SerratedRect2DDetector.hpp"

SerratedRect2DDetector::SerratedRect2DDetector(unsigned nbr_of_strips,
		unsigned strip_length, unsigned strip_width, unsigned half_pitch,
		double strip_potential, unsigned refine_level, unsigned max_iter,
		double stop_accuracy) :
		SerratedRect2DDetector(nbr_of_strips, DEFAULT_RECT_WIDTH, strip_length,
				strip_width, half_pitch, strip_potential, refine_level,
				max_iter, stop_accuracy) {
}

SerratedRect2DDetector::SerratedRect2DDetector(unsigned nbr_of_strips,
		unsigned width, unsigned strip_length, unsigned strip_width,
		unsigned half_pitch, double strip_potential, unsigned refine_level,
		unsigned max_iter, double stop_accuracy) {

	this->strip_potential = strip_potential;
	this->refine_level = refine_level;
	this->max_iter = max_iter;
	this->stop_accuracy = stop_accuracy;

	serr_geo_info = new SerratedRectGeoInfo(2, nbr_of_strips, width,
			strip_length, strip_width, half_pitch);
	geo_info = serr_geo_info;

	MyGridGenerator<2>::serrated_rectangle(*triangulation, width, nbr_of_strips,
			strip_length, strip_width, half_pitch);
	boundary_conditions = new SerratedRect2DBoundaryCond<2>(serr_geo_info,
			strip_potential);
	potential_solver = new LaplaceSolver<2>(triangulation, refine_level,
			max_iter, stop_accuracy, zero_right_hand_side, boundary_conditions,
			true);

	boundary_conditions_weight = new SerratedRect2DBoundaryCondWeight<2>(
			serr_geo_info, strip_potential);
	MyGridGenerator<2>::serrated_rectangle(*triangulation_weight, width,
			nbr_of_strips, strip_length, strip_width, half_pitch);
	potential_solver_weight = new LaplaceSolver<2>(triangulation_weight,
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

std::string SerratedRect2DDetector::params_to_string() {

	std::string str = "width" + std::to_string(geo_info->get_width())
			+ "_nbr_of_strips_" + std::to_string(geo_info->get_nbr_of_strips())
			+ "_strip_length_"
			+ std::to_string(serr_geo_info->get_strip_length())
			+ "_strip_width_" + std::to_string(serr_geo_info->get_strip_width())
			+ "_half-pitch_" + std::to_string(serr_geo_info->get_half_pitch())
			+ "_strip_potential_" + std::to_string(strip_potential)
			+ "_refine_level_" + std::to_string(refine_level) + "_max_iter_"
			+ std::to_string(max_iter) + "_stop_accuracy_"
			+ std::to_string(stop_accuracy);
	return str;
}
