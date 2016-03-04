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
#include "Constants.hpp"

using namespace dealii;

template<unsigned dim>
class ElectrodeCurrent {

public:

	ElectrodeCurrent(unsigned strip_potential, MyGeometryInfo *geo_info,
			Solution<dim> *solution, Solution<dim> *solution_weight,
			Line particle_trajectory, unsigned refine_level) :
			e(TYPE_SILICIUM), h(TYPE_SILICIUM) {
		this->particle_trajectory = particle_trajectory;
		common_constructor(strip_potential, geo_info, solution, solution_weight,
				refine_level);
	}

	ElectrodeCurrent(unsigned strip_potential, MyGeometryInfo *geo_info,
			Solution<dim> *solution, Solution<dim> *solution_weight,
			unsigned refine_level) :
			e(TYPE_SILICIUM), h(TYPE_SILICIUM) {

		this->particle_trajectory = geo_info->get_mid_length_vertical_line();
		common_constructor(strip_potential, geo_info, solution, solution_weight,
				refine_level);
	}

	void print_charges() {
		unsigned nbr_of_charges = punctual_charges.size();
		for (unsigned i = 0; i < nbr_of_charges; i++) {
			std::pair<Point<2>, Charge*> el = punctual_charges.front();
			punctual_charges.pop();
			Point<2> p = el.first;
			std::cout << "(" << p[0] << "," << p[1] << ") ";
			punctual_charges.push(el);
		}
		std::cout << std::endl;
	}

	/**
	 * @pre: delta_t > 0
	 *
	 * @return:
	 * - true if all charges have left the detector at the end of simultion.
	 * - false if the simulation has stopped because the speed of every
	 * charges with zero. It may happen, for example, if the detector potential
	 * difference is 0V.
	 *
	 */
	bool compute_current(double delta_t,
			std::vector<std::pair<double, double> > &current_vs_time) {

		if (delta_t <= 0.0)
			throw PRECONDITIONS_VIOLATED;

		std::cout << "Delta time: " << delta_t << "s" << std::endl;

		double time = 0.0;
		bool no_moves = false;

		while (!punctual_charges.empty() && !no_moves) {
			double current = move_charges(delta_t, no_moves);
			std::pair<double, double> point(time, current);
			current_vs_time.push_back(point);
			time += delta_t;
		}
		return no_moves;
	}

	bool compute_current(
			std::vector<std::pair<double, double> > &current_vs_time) {
		double delta_t = compute_delta_t();
		return compute_current(delta_t, current_vs_time);
	}

private:

	unsigned refine_level;
	MyGeometryInfo *geo_info;
	Solution<dim> *laplace_sol, *laplace_sol_weight;
	Line particle_trajectory;
	double hole_pairs_nbr_per_lgth = 75; //per microm
	double strip_potential;
	double covered_dist;

	void common_constructor(unsigned strip_potential, MyGeometryInfo *geo_info,
			Solution<dim> *solution, Solution<dim> *solution_weight,
			unsigned refine_level) {
		this->geo_info = geo_info;
		this->laplace_sol = solution;
		this->laplace_sol_weight = solution_weight;
		this->refine_level = refine_level;
		this->strip_potential = strip_potential;

		std::vector<Point<2>> intersect = get_trajectory_intersect();
		covered_dist = dist_covered_inside_det(intersect);
		place_initial_charges(intersect, covered_dist);
	}

	/**
	 * A punctual charge is defined as a point where resides an electric
	 * charge whose value is total charge inside the detector resulting
	 * from the particle crossing divided by the number of punctual charges.
	 *
	 * Punctual charges are spread uniformly on particle trajectory.
	 *
	 */
	std::queue<std::pair<Point<2>, Charge*>> punctual_charges;
	double punctual_electric_charge; // electric charge at each punctual charge (C/e)
	Electron e;
	Hole h;

	/**
	 * Compute total distance covered by the particle INSIDE the detector.
	 *
	 * @pre: intersect.size() is even.
	 *
	 */
	double dist_covered_inside_det(std::vector<Point<2>> intersect) {

		if (Utils::is_odd(intersect.size()))
			throw PRECONDITIONS_VIOLATED;

		double total_dist = 0.0;
		for (unsigned i = 0; i < intersect.size(); i += 2) {
			Segment seg(intersect[i], intersect[i + 1]);
			total_dist += seg.compute_length();
		}
		return total_dist;
	}

	std::vector<Point<2>> get_trajectory_intersect() {
		//Get all the points of intersection between the line
		//(particle trajectory) and the boundaries of the detector
		std::vector<Point<2>> intersect = geo_info->boundaries_intersections(
				particle_trajectory);

		//for (unsigned i = 0; i < intersect.size(); i++) {
		//	Utils::print_point<2>(intersect[i]);
		//}
		//Should not happen
		if (Utils::is_odd(intersect.size())) {
			std::cout << "Warning odd intersections number (initial_charges): "
					<< intersect.size() << std::endl;
			throw PRECONDITIONS_VIOLATED;
		}
		return intersect;
	}

	void place_initial_charges(std::vector<Point<2>> intersect,
			double dist_covered_by_particle) {

		std::cout << "dist_covered_by_particle " << dist_covered_by_particle << std::endl;

		unsigned nbr_of_punctual_charges = std::pow(2, refine_level);
		double dist_between_punctual_charges = dist_covered_by_particle
				/ (nbr_of_punctual_charges + 1);

		std::cout << "dist_between_punctual_charges: " << dist_between_punctual_charges << std::endl;

		double total_hole_pairs_nbr = dist_covered_by_particle
				* hole_pairs_nbr_per_lgth;
		std::cout << "Hole pairs number: " << total_hole_pairs_nbr << std::endl;
		punctual_electric_charge = total_hole_pairs_nbr
				/ nbr_of_punctual_charges;

		std::cout << "punctual_electric_charge: " << punctual_electric_charge << std::endl;

		for (unsigned i = 0; i < intersect.size(); i += 2) {
			Point<2> particle_entry = intersect[i];
			Point<2> particle_exit = intersect[i + 1];
			Segment seg(particle_entry, particle_exit);
			place_initial_charges_on_segment(seg, dist_between_punctual_charges,
					dist_covered_by_particle, nbr_of_punctual_charges);
		}
	}

	void place_initial_charges_on_segment(Segment &seg,
			double dist_between_punctual_charges,
			double dist_covered_by_particle, unsigned nbr_of_punctual_charges) {

		unsigned nbr_of_punctual_charges_on_seg = ceil(
				seg.compute_length() / dist_covered_by_particle
						* nbr_of_punctual_charges);

		std::cout << "nbr_of_charges_on_seg " << nbr_of_punctual_charges_on_seg << std::endl;

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

		for (unsigned i = 0; i < nbr_of_punctual_charges_on_seg; i++) {
			Point<2> p(x_coord, y_coord);
			std::pair<Point<2>, Charge*> holes(p, &h);
			std::pair<Point<2>, Charge*> electrons(p, &e);
			punctual_charges.push(holes);
			punctual_charges.push(electrons);
			x_coord += delta_x;
			y_coord += delta_y;
		}
	}

	Tensor<1, 2> puncutal_charge_speed(Point<2> pos, Charge *charge) {
		//Formula in latex: \vec{v} = \mu \vec{E}
		PhysicalValues<2> val = laplace_sol->get_values(pos);
		Tensor<1, 2> electric_field = val.electric_field;
		//std::cout << "electric_field (" << electric_field[0] << "," << electric_field[1] << ") ";
		int electric_charge_sign = charge->get_charge_sign();
		double mobility = charge->get_mobility_saturation(electric_field);
		return electric_charge_sign * electric_field * mobility;		//*0.86;
	}

	double current(Point<2> pos, Tensor<1, 2> speed, Charge *charge) {
		//Formula in latex: \vec{i} = -q \vec{v} \cdot \vec{E_w}
		PhysicalValues<2> val = laplace_sol_weight->get_values(pos);
		Tensor<1, 2> weighting_electric_field = val.electric_field;
		//std::cout << "WF (" << weighting_electric_field[0] << "," << weighting_electric_field[1] << ") " << std::endl;
		int electric_charge_sign = charge->get_charge_sign();
		double current = -electric_charge_sign * punctual_electric_charge
				* speed * weighting_electric_field;
		return current;
	}

	Point<2> compute_new_pos(Point<2> prev_pos, Tensor<1, 2> speed,
			double delta_t) {
		double new_x = prev_pos[0] + speed[0] * delta_t;
		double new_y = prev_pos[1] + speed[1] * delta_t;
		Point<2> new_pos(new_x, new_y);
		return new_pos;
	}

	double move_charges(double delta_t, bool &no_moves) {

		double current_tot = 0.0;
		unsigned init_nbr_punctual_charges = punctual_charges.size();
		no_moves = true;
		for (unsigned i = 0; i < init_nbr_punctual_charges; i++) {

			std::pair<Point<2>, Charge*> punct_charge =
					punctual_charges.front();
			punctual_charges.pop();
			Point<2> pos = punct_charge.first;
			Charge *charge = punct_charge.second;
			Tensor<1, 2> speed = puncutal_charge_speed(pos, charge);
			if (!TensorUtils::is_zero_tensor<2>(speed))
				no_moves = false;

			current_tot += current(pos, speed, charge);
			Point<2> new_pos = compute_new_pos(pos, speed, delta_t);

			if (geo_info->is_point_inside_geometry(new_pos)) {
				std::pair<Point<2>, Charge*> new_charge(new_pos, charge);
				punctual_charges.push(new_charge);
			}
		}

		return current_tot * ELECTRON_CHARGE;
	}

	double estimate_collection_time() {
		//t = d^2/(mu V)
		double max_mobility = h.get_mobility_300K();
		double dist_to_strip = geo_info->get_width()-geo_info->get_strip_width();
		return std::pow(dist_to_strip,2) / (max_mobility * strip_potential);
	}

	double compute_delta_t() {
		double col_time = estimate_collection_time();
		std::cout << "Predicted collection time: " << col_time << std::endl;
		return col_time / std::pow(2, refine_level);
	}
};
