/*
 * ElectrodeCurrent.hpp
 *
 *  Created on: 11 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <forward_list>
#include <deal.II/dofs/dof_handler.h>

#include "geometry_info/MyGeometryInfo.hpp"
#include "Utils.hpp"

using namespace dealii;

template<unsigned dim>
class ElectrodeCurrent {

public:

	ElectrodeCurrent(MyGeometryInfo *geo_info, Solution<dim> *potential,
			Solution<dim> *weight_potential,
			std::vector<
					std::pair<typename DoFHandler<dim>::active_cell_iterator,
							std::vector<Tensor<1, dim> > > > *electric_field,
			std::vector<
					std::pair<typename DoFHandler<dim>::active_cell_iterator,
							std::vector<Tensor<1, dim> > > > *electric_field_weight,
			Line *particle_trajectory, unsigned refine_level) {
		this->geo_info = geo_info;
		this->potential = potential;
		this->weight_potential = weight_potential;
		this->electric_field = electric_field;
		this->electric_field_weight = electric_field_weight;
		this->particle_trajectory = particle_trajectory;
		this->refine_level = refine_level;

		place_initial_charges();
	}

private:

	double hole_pairs_nbr_per_lgth = 80; //per microm
	unsigned refine_level;

	MyGeometryInfo *geo_info;
	Solution<dim> *potential, *weight_potential;
	std::vector<
			std::pair<typename DoFHandler<dim>::active_cell_iterator,
					std::vector<Tensor<1, dim> > > > *electric_field,
			*electric_field_weight;
	Line *particle_trajectory;

	//bool type: true if hole, false if electron
	std::forward_list<std::pair<Point<2>, bool>> charges;

	/**
	 * Compute total distance covered by the particle INSIDE the detector.
	 *
	 * @pre: intersect.size() is even.
	 *
	 */
	double dist_covered_inside_det(std::vector<Point<2>> intersect) {

		double total_dist = 0.0;
		for (unsigned i = 0; i < intersect.size(); i += 2) {
			Segment seg(intersect[i], intersect[i + 1]);
			total_dist += seg.compute_length();
		}
		return total_dist;
	}

	void place_initial_charges() {

		//Get all the points of intersection between the line
		//(particle trajectory) and the boundaries of the detector
		std::vector<Point<2>> intersect = geo_info->boundaries_intersections(
				*particle_trajectory);

		//
		for(unsigned i=0; i<intersect.size(); i++){
			std::cout << "(" << intersect[i][0] << "," << intersect[i][1] << ")";
		}

		//

		//TODO: remove latter
		if (intersect.size() % 2 != 0) {
			std::cout << "Warning odd intersections number (initial_charges): "
				<< intersect.size()	<< std::endl;

		}

		double covered_dist = dist_covered_inside_det(intersect);
		unsigned nbr_of_punctual_charges = std::pow(2, refine_level);
		double dist_between_punctual_charges = covered_dist
				/ (nbr_of_punctual_charges + 1);

		for (unsigned i = 0; i < intersect.size(); i += 2) {
			Point<2> particle_entry = intersect[i];
			Point<2> particle_exit = intersect[1];
			Segment seg(particle_entry, particle_exit);
			place_initial_charges_on_segment(seg, dist_between_punctual_charges,
					covered_dist, nbr_of_punctual_charges);
		}
	}

	void place_initial_charges_on_segment(Segment seg,
			double dist_between_punctual_charges, double covered_dist,
			unsigned nbr_of_punctual_charges) {

		unsigned nbr_of_charges_on_seg = ceil(
				seg.compute_length() / covered_dist * nbr_of_punctual_charges);

		double delta_x;
		double delta_y;

		if (seg.is_vertical()) {
			delta_x = 0.0;
			delta_y = dist_between_punctual_charges;
		} else {
			double slope = seg.get_slope();
			delta_x = dist_between_punctual_charges
					/ std::sqrt(1 + slope * slope);
			delta_y = slope * delta_x;
		}

		double x_coord = seg.p1[0] + delta_x;
		double y_coord = seg.p1[1] + delta_y;

		for (unsigned i = 0; i < nbr_of_charges_on_seg; i++) {
			Point<2> p(x_coord, y_coord);
			std::pair<Point<2>, bool> holes(p, true);
			std::pair<Point<2>, bool> electrons(p, false);
			charges.push_front(holes);
			charges.push_front(electrons);
			x_coord += delta_x;
			y_coord += delta_y;
		}
	}

};
