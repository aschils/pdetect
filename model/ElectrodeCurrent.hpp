/*
 * ElectrodeCurrent.hpp
 *
 *  Created on: 11 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <deal.II/dofs/dof_handler.h>

#include "geometry_info/MyGeometryInfo.hpp"
#include "Utils.hpp"

using namespace dealii;

template<unsigned dim>
class ElectrodeCurrent {

public:

	std::vector<Point<dim>, unsigned> charges;

	ElectrodeCurrent(MyGeometryInfo *geo_info, Solution<dim> *potential,
			Solution<dim> *weight_potential,
			std::vector<
					std::pair<typename DoFHandler<dim>::active_cell_iterator,
							std::vector<Tensor<1, dim> > > > *electric_field,
			std::vector<
					std::pair<typename DoFHandler<dim>::active_cell_iterator,
							std::vector<Tensor<1, dim> > > > *electric_field_weight,
			Line particle_trajectory) {
		this->geo_info = geo_info;
		this->electric_field = electric_field;
		this->electric_field_weight = electric_field_weight;
		this->particle_trajectory = particle_trajectory;
	}

	void initial_charges(){

	}

private:

	double hole_pairs_nbr = 80; //per microm

	MyGeometryInfo *geo_info;
	Solution<dim> *potential, *weight_potential;
	std::vector<
			std::pair<typename DoFHandler<dim>::active_cell_iterator,
					std::vector<Tensor<1, dim> > > > *electric_field,
			*electric_field_weight;
	Line particle_trajectory;

};
