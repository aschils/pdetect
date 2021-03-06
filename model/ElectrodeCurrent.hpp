/*
 * ElectrodeCurrent.hpp
 *
 *  Created on: 11 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <queue>
#include <deal.II/dofs/dof_handler.h>
#include <boost/geometry.hpp>

#include "geometry_info/MyGeometryInfo.hpp"
#include "Utils.hpp"
#include "Charge.hpp"
#include "Constants.hpp"

using namespace dealii;

template<unsigned dim>
class ElectrodeCurrent {

public:

	typedef bg::model::point<double, dim, bg::cs::cartesian> bpoint;

	ElectrodeCurrent(Detector2D *det, Line particle_trajectory,
			unsigned refine_level) {
		this->particle_trajectory = particle_trajectory;
		common_constructor(det, refine_level);
	}

	ElectrodeCurrent(Detector2D *det, unsigned refine_level) {
		this->particle_trajectory =
				det->get_geometry_info()->get_mid_length_vertical_line();
		common_constructor(det, refine_level);
	}

	void print_charges() {
		unsigned nbr_of_charges = punctual_charges.size();
		for (unsigned i = 0; i < nbr_of_charges; i++) {
			std::tuple<bpoint, Charge*, double> el = punctual_charges.front();
			punctual_charges.pop();
			bpoint p = std::get<0>(el);
			std::cout << "(" << bg::get<0>(p) << "," << bg::get<1>(p) << ") ";
			punctual_charges.push(el);
		}
		std::cout << std::endl;
	}

	/**
	 * @pre: delta_t > 0
	 *
	 * @return:
	 * - true if all charges have left the detector at the end of simulation.
	 * - false if the simulation has stopped because the speed of every
	 * charges is zero. It may happen, for example, if the detector potential
	 * difference is 0V.
	 *
	 */
	bool compute_current(double delta_t,
			std::vector<std::pair<double, double> > &current_vs_time) {
		if (delta_t <= 0.0)
			throw PRECONDITIONS_VIOLATED;
		return compute_current(true, delta_t, current_vs_time);
	}

	bool compute_current(
			std::vector<std::pair<double, double> > &current_vs_time) {
		return compute_current(false, 0.0, current_vs_time);
	}

private:

	unsigned refine_level;
	Detector2D *det;
	MyGeometryInfo *geo_info;
	Solution<dim> laplace_sol, laplace_sol_weight;
	Line particle_trajectory;
	double hole_pairs_nbr_per_lgth; //per microm
	double strip_potential;
	double covered_dist;
	double coll_t_to_delta_t_div;
	double tot_electrons;

	Electron electron;
	Hole hole;

	void common_constructor(Detector2D *det, unsigned refine_level) {
		this->det = det;
		coll_t_to_delta_t_div = std::pow(2, refine_level);
		this->geo_info = det->get_geometry_info();
		det->get_solution(laplace_sol);
		det->get_solution_weight(laplace_sol_weight);
		this->refine_level = refine_level;
		this->strip_potential = det->get_strip_potential();
		electron = det->get_electron();
		hole = det->get_hole();
		this->hole_pairs_nbr_per_lgth = det->get_hole_pairs_nbr_per_lgth();

		std::vector<bpoint> intersect = get_trajectory_intersect();
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
	//std::queue<std::pair<Point<2>, Charge*>> punctual_charges;
	std::queue<std::tuple<bpoint, Charge*, double>> punctual_charges;
	//double punctual_electric_charge; // electric charge at each punctual charge (C/e)

	bool compute_current(bool user_delta_t_on, double user_delta_t,
			std::vector<std::pair<double, double> > &current_vs_time) {

		bool no_moves = false;
		double time = 0.0;
		double delta_t =
				(user_delta_t_on) ?
						user_delta_t : adaptive_delta_t(LIGHT_SPEED);

		while (!punctual_charges.empty() && !no_moves) {

			double max_speed_y = 0.0;
			double current = move_charges(delta_t, no_moves, max_speed_y);
			std::pair<double, double> point(time, current);
			current_vs_time.push_back(point);

			time += delta_t;
			if (!user_delta_t_on)
				delta_t = adaptive_delta_t(max_speed_y);
		}

		return no_moves;
	}

	/**
	 * Compute total distance covered by the particle INSIDE the detector.
	 *
	 * @pre: intersect.size() is even.
	 *
	 */
	double dist_covered_inside_det(std::vector<bpoint> intersect) {

		if (Utils::is_odd(intersect.size()))
			throw PRECONDITIONS_VIOLATED;

		double total_dist = 0.0;
		for (unsigned i = 0; i < intersect.size(); i += 2) {
			bg::model::segment<bpoint> seg(intersect[i], intersect[i + 1]);
			total_dist += bg::length(seg);
		}
		return total_dist;
	}

	std::vector<bpoint> get_trajectory_intersect() {
		//Get all the points of intersection between the line
		//(particle trajectory) and the boundaries of the detector
		std::vector<bpoint> intersect = geo_info->boundaries_intersections(
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

	void place_initial_charges(std::vector<bpoint> intersect,
			double dist_covered_by_particle) {

		std::cout << "dist_covered_by_particle " << dist_covered_by_particle
				<< std::endl;

		unsigned nbr_of_punctual_charges = std::pow(2, refine_level);
		double dist_between_punctual_charges = dist_covered_by_particle
				/ (nbr_of_punctual_charges + 1);

		std::cout << "dist_between_punctual_charges: "
				<< dist_between_punctual_charges << std::endl;

		double total_hole_pairs_nbr = dist_covered_by_particle
				* hole_pairs_nbr_per_lgth;
		std::cout << "Hole pairs number: " << total_hole_pairs_nbr << std::endl;
		tot_electrons = total_hole_pairs_nbr;
		double init_punctual_electric_charge = total_hole_pairs_nbr
				/ nbr_of_punctual_charges;

		std::cout << "init_punctual_electric_charge: "
				<< init_punctual_electric_charge << std::endl;

		for (unsigned i = 0; i < intersect.size(); i += 2) {
			bpoint particle_entry = intersect[i];
			bpoint particle_exit = intersect[i + 1];
			bg::model::segment<bpoint> seg(particle_entry, particle_exit);
			place_initial_charges_on_segment(seg, dist_between_punctual_charges,
					dist_covered_by_particle, nbr_of_punctual_charges,
					init_punctual_electric_charge);
		}
	}

	void place_initial_charges_on_segment(bg::model::segment<bpoint> &seg,
			double dist_between_punctual_charges,
			double dist_covered_by_particle, unsigned nbr_of_punctual_charges,
			double init_punctual_electric_charge) {

		std::cout << "seg length: " << bg::length(seg) << std::endl;

		unsigned nbr_of_punctual_charges_on_seg = ceil(
				bg::length(seg) / dist_covered_by_particle
						* nbr_of_punctual_charges);

		std::cout << "nbr_of_charges_on_seg " << nbr_of_punctual_charges_on_seg
				<< std::endl;

		for (unsigned i = 0; i < nbr_of_punctual_charges_on_seg; i++) {
			double dist_from_seg_first_pt = dist_between_punctual_charges
					* (i + 1);
			bpoint p = point_on_seg(seg, dist_from_seg_first_pt);
			//std::cout << "(" << bg::get<0>(p) << "," << bg::get<1>(p) << ")" << std::endl;
			std::tuple<bpoint, Charge*, double> holes(p, &hole,
					init_punctual_electric_charge);
			std::tuple<bpoint, Charge*, double> electrons(p, &electron,
					init_punctual_electric_charge);
			punctual_charges.push(holes);
			punctual_charges.push(electrons);
		}

		/*
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
		 std::tuple<Point<2>, Charge*, double> holes(p, &hole,
		 init_punctual_electric_charge);
		 std::tuple<Point<2>, Charge*, double> electrons(p, &electron,
		 init_punctual_electric_charge);
		 punctual_charges.push(holes);
		 punctual_charges.push(electrons);
		 x_coord += delta_x;
		 y_coord += delta_y;
		 }*/

	}

	/**
	 * Let seg be the segment between the points p1 and p2:
	 *
	 * p1 .--------. p2
	 *
	 * Find coordinates of point p located on the segment seg and at distance
	 * dist_from_p1 from p1
	 *
	 */
	bpoint point_on_seg(bg::model::segment<bpoint> &seg, double dist_from_p1) {

		double seg_length = bg::length(seg);

		assert(dist_from_p1 <= seg_length);
		assert(seg_length > 0.0);
		assert(dist_from_p1 >= 0.0);

		bpoint p1 = seg.first;
		bpoint p2 = seg.second;
		std::cout << "p2: (" << bg::get<0>(p2) << "," << bg::get<1>(p2) << ")" << std::endl;
		// direction vector v = p1-p2
		bpoint direction_vec = p1;
		std::cout << "p1: (" << bg::get<0>(direction_vec) << "," << bg::get<1>(direction_vec) << ")" << std::endl;
		bg::subtract_point(direction_vec, p2);
		std::cout << "dir_vec: (" << bg::get<0>(direction_vec) << "," << bg::get<1>(direction_vec) << ")" << std::endl;
		double dir_vec_norm =
				sqrt(Utils::DimRecursion<dim, dim>::norm_squared(direction_vec));
		std::cout << "norm squared dir vec: " << dir_vec_norm << std::endl;

		/*
		 * parametric equations for lines:
		 * P_x = lambda * v_x + P0_x
		 * P_y = lambda * v_y + P0_y
		 *
		 * Find lambda pour P0 = p2 and P = p1
		 *
		 * lambda = (x + y - P0_x - P0_y)/(v_x+v_y)
		 *
		 * Uses both equations to isolate lambda avoids dividing by 0 when
		 * either v_x = 0 or v_y = 0.
		 *
		 * Note: lambda_p1 is either +1 or -1
		 */
		double lambda_p1 = bg::get<0>(p1)-bg::get<0>(p2)+bg::get<1>(p1)-bg::get<1>(p2);
		lambda_p1 = lambda_p1/(bg::get<0>(direction_vec)+bg::get<1>(direction_vec));

		//TODO may not work in all cases: must chose automatically correct sign
		//lambda_p1 +/- dist_from_p1/dir_vec_norm
		double lambda = lambda_p1 - dist_from_p1/dir_vec_norm;
		double x = lambda * bg::get<0>(direction_vec) + bg::get<0>(p2);
		double y = lambda * bg::get<1>(direction_vec) + bg::get<1>(p2);
		bpoint point_on_seg(x, y);
		return point_on_seg;
	}

	double adaptive_delta_t(const double &max_speed_y) {
		return geo_info->get_width() / max_speed_y / coll_t_to_delta_t_div;
	}

	Tensor<1, 2> puncutal_charge_speed(Tensor<1, 2> &electric_field,
			Charge *charge) {
		//Formula in latex: \vec{v} = \mu \vec{E}
		//std::cout << "electric_field (" << electric_field[0] << "," << electric_field[1] << ") ";
		int electric_charge_sign = charge->get_charge_sign();
		double mobility = charge->get_mobility_saturation(electric_field);
		return electric_charge_sign * electric_field * mobility;		//*0.86;
	}

	double current(bpoint &pos, Tensor<1, 2> &speed, Charge *charge,
			double &punctual_electric_charge) {
		//Formula in latex: \vec{i} = -q \vec{v} \cdot \vec{E_w}
		Point<dim> dpos;
		Utils::DimRecursion<dim, dim>::bpoint_to_deailii_point(pos, dpos);
		PhysicalValues<2> val = laplace_sol_weight.get_values(dpos);
		Tensor<1, 2> weighting_electric_field = val.electric_field;
		//std::cout << "WF (" << weighting_electric_field[0] << "," << weighting_electric_field[1] << ") " << std::endl;
		int electric_charge_sign = charge->get_charge_sign();
		double current = -electric_charge_sign * punctual_electric_charge
				* speed * weighting_electric_field;
		return current;
	}

	bpoint compute_new_pos(bpoint &prev_pos, Tensor<1, 2> &speed,
			double &delta_t) {
		double new_x = bg::get<0>(prev_pos) + speed[0] * delta_t;
		double new_y = bg::get<1>(prev_pos) + speed[1] * delta_t;
		bpoint new_pos(new_x, new_y);
		return new_pos;
	}

	void generate_new_charges(Charge *charge, bpoint &pos, bpoint &new_pos,
			PhysicalValues<2> &values_at_pos, double &electric_charge) {

		double displacement = (new_pos - pos).norm();
		double first_townsend_coefficient = det->get_first_townsend_coefficient(
				pos, values_at_pos);
		double new_electric_charge = electric_charge
				* exp(first_townsend_coefficient * displacement);

		double new_electrons = (new_electric_charge - electric_charge);
		tot_electrons += new_electrons;

		std::tuple<bpoint, Charge*, double> new_charge(new_pos, charge,
				new_electric_charge);
		punctual_charges.push(new_charge);
	}

	double move_charges(double &delta_t, bool &no_moves, double &max_speed_y) {

		double current_tot = 0.0;
		unsigned init_nbr_punctual_charges = punctual_charges.size();
		no_moves = true;
		max_speed_y = 0.0;

		for (unsigned i = 0; i < init_nbr_punctual_charges; i++) {

			std::cout << "loop in move charges " << i << std::endl;

			std::tuple<bpoint, Charge*, double> punct_charge =
					punctual_charges.front();
			punctual_charges.pop();

			bpoint pos = std::get<0>(punct_charge);
			Charge *charge = std::get<1>(punct_charge);
			double electric_charge = std::get<2>(punct_charge);
			//bool townsend_avalanche_done = std::get<3>(punct_charge);

			Point<dim> dp;
			Utils::DimRecursion<dim, dim>::bpoint_to_deailii_point(pos, dp);

			std::cout << "bpoint: (" << bg::get<0>(pos) << ","
					<< bg::get<1>(pos) << ")" << std::endl;
			std::cout << "dealiipoint: (" << dp[0] << "," << dp[1] << ")"
					<< std::endl;

			std::cout << "after point type conv " << i << std::endl;
			PhysicalValues<2> values_at_pos = laplace_sol.get_values(dp);
			std::cout << "after lapalce sol. get val " << i << std::endl;
			Tensor<1, 2> electric_field = values_at_pos.electric_field;

			std::cout << "after tensor EF" << i << std::endl;

			Tensor<1, 2> speed = puncutal_charge_speed(electric_field, charge);
			if (no_moves && !TensorUtils::is_zero_tensor<2>(speed))
				no_moves = false;

			double abs_speed_y = fabs(speed[1]);
			if (abs_speed_y > max_speed_y)
				max_speed_y = abs_speed_y;

			current_tot += current(pos, speed, charge, electric_charge);

			bpoint new_pos = compute_new_pos(pos, speed, delta_t);

			if (geo_info->is_point_inside_geometry(new_pos)) {

				//int charge_sign = charge->get_charge_sign();
				//if(charge_sign < 0 && tot_electrons < 1e12) { //if electron, townsend avalanche
				//	generate_new_charges(charge, pos, new_pos, values_at_pos,
				//			electric_charge);
				//}
				//else { //if hole
				std::tuple<bpoint, Charge*, double> new_charge(new_pos, charge,
						electric_charge);
				punctual_charges.push(new_charge);
				//}
			}
		}
		//if(current_tot*ELECTRON_CHARGE < -1e-05)
		//		current_tot = -1e-05/ELECTRON_CHARGE;
		return current_tot * ELECTRON_CHARGE;
	}

	/*
	 double estimate_collection_time() {
	 //t = d^2/(mu V)
	 double max_mobility = hole.get_mobility_300K();
	 return std::pow(dist_to_strip, 2) / (max_mobility * strip_potential);
	 }

	 double initial_delta_t() {
	 double col_time = estimate_collection_time();
	 std::cout << "Predicted collection time: " << col_time << std::endl;
	 return col_time / std::pow(2, refine_level);
	 }*/
};
