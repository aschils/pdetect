/*
 * ElectrodeCurrent.hpp
 *
 *  Created on: 11 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <queue>
#include <deal.II/dofs/dof_handler.h>

#include "geometry_info/MyGeometryInfo.hpp"
#include "Utils.hpp"
#include "Charge.hpp"

using namespace dealii;

template<unsigned dim>
class ElectrodeCurrent {

public:

	ElectrodeCurrent(MyGeometryInfo *geo_info, Solution<dim> *solution,
			Solution<dim> *solution_weight, Line *particle_trajectory,
			unsigned refine_level): e(TYPE_SILICIUM), h(TYPE_SILICIUM) {
		this->geo_info = geo_info;
		this->laplace_sol = solution;
		this->laplace_sol_weight = solution_weight;
		this->particle_trajectory = particle_trajectory;
		this->refine_level = refine_level;

		place_initial_charges();
	}

	void print_charges() {
		for (unsigned i = 0; i < charges.size(); i++) {
			std::pair<Point<2>, Charge*> el = charges.front();
			charges.pop();
			Point<2> p = el.first;
			std::cout << "(" << p[0] << "," << p[1] << ") ";
			charges.push(el);
		}
		std::cout << std::endl;
	}

	/**
	 * @pre delta_t > 0
	 */
	void compute_current(double delta_t,
			std::vector<std::pair<double, double> > &current_vs_time) {

		if (delta_t < 0.0)
			throw PRECONDITIONS_VIOLATED;

		double time = 0.0;

		while (!charges.empty()) {
			double current = move_charges(delta_t);
			std::pair<double, double> point(time, current);
			current_vs_time.push_back(point);
			time += delta_t;
		}
	}

private:

	double hole_pairs_nbr_per_lgth = 80; //per microm
	unsigned refine_level;
	MyGeometryInfo *geo_info;
	Solution<dim> *laplace_sol, *laplace_sol_weight;
	Line *particle_trajectory;
	std::queue<std::pair<Point<2>, Charge*>> charges;
	double ponctual_charge;// C/e
	Electron e;
	Hole h;
	double electron_charge = 1.60217662e-13; // microCoulomb


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
		std::cout << "intersections:" << std::endl;
		for (unsigned i = 0; i < intersect.size(); i++) {
			std::cout << "(" << intersect[i][0] << "," << intersect[i][1]
					<< ")";
		}
		std::cout << std::endl;
		//

		//Should not happen
		if (intersect.size() % 2 != 0) {
			std::cout << "Warning odd intersections number (initial_charges): "
					<< intersect.size() << std::endl;
		}

		double covered_dist = dist_covered_inside_det(intersect);
		unsigned nbr_of_punctual_charges = std::pow(2, refine_level);
		double dist_between_punctual_charges = covered_dist
				/ (nbr_of_punctual_charges + 1);
		double total_charge = covered_dist * hole_pairs_nbr_per_lgth;
		ponctual_charge = total_charge / nbr_of_punctual_charges;

		for (unsigned i = 0; i < intersect.size(); i += 2) {
			Point<2> particle_entry = intersect[i];
			Point<2> particle_exit = intersect[i + 1];
			Segment seg(particle_entry, particle_exit);
			place_initial_charges_on_segment(seg, dist_between_punctual_charges,
					covered_dist, nbr_of_punctual_charges);
		}
	}

	void place_initial_charges_on_segment(Segment &seg,
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

		//std::cout << "delta_x " << delta_x << " delta_y " << delta_y << std::endl;

		double x_coord = seg.p1[0] + delta_x;
		double y_coord = seg.p1[1] + delta_y;

		for (unsigned i = 0; i < nbr_of_charges_on_seg; i++) {
			Point<2> p(x_coord, y_coord);
			std::pair<Point<2>, Charge*> holes(p, &h);
			std::pair<Point<2>, Charge*> electrons(p, &e);
			charges.push(holes);
			charges.push(electrons);
			x_coord += delta_x;
			y_coord += delta_y;
		}
	}

	Tensor<1,2> poncutal_charge_speed(Point<2> pos, Charge *charge){
		//Formula in latex: \vec{v} = \mu \vec{E}

		ValuesAtPoint<2> val = laplace_sol->get_values(pos);
		Tensor<1,2> electric_field = val.gradient;
		return charge->get_charge_sign()*electric_field*charge->get_mobility();
	}

	double current(Point<2> pos, Tensor<1,2> speed, Charge *charge){
		//Formula in latex: \vec{i} = -q \vec{v} \cdot \vec{E_w}
		//But current is negative with this, formula therefore
		//we implement \vec{i} = +q \vec{v} \cdot \vec{E_w}
		ValuesAtPoint<2> val = laplace_sol_weight->get_values(pos);
		Tensor<1,2> weighting_electric_field = val.gradient;
		double current = -charge->get_charge_sign()
				*ponctual_charge*(speed*weighting_electric_field);
		return current;
	}

	Point<2> compute_new_pos(Point<2> prev_pos, Tensor<1,2> speed,
			double delta_t){
		double new_x = prev_pos[0]+speed[0]*delta_t;
		double new_y = prev_pos[1]+speed[1]*delta_t;
		Point<2> new_pos(new_x, new_y);
		return new_pos;
	}

	double move_charges(double delta_t) {

		double current_tot = 0.0;

		for(unsigned i=0; i<charges.size(); i++){
			std::pair<Point<2>, Charge*> ponct_charge = charges.front();
			charges.pop();
			Point<2> pos = ponct_charge.first;
			Charge *charge = ponct_charge.second;
			Tensor<1,2> speed = poncutal_charge_speed(pos, charge);
			current_tot += current(pos, speed, charge);
			Point<2> new_pos = compute_new_pos(pos, speed, delta_t);

			if(geo_info->is_point_inside_geometry(
					Utils::point_to_std_vec<2>(new_pos))){
				std::pair<Point<2>, Charge*> new_charge(new_pos, charge);
				charges.push(new_charge);
			}
		}

		return current_tot*electron_charge;
	}

};
