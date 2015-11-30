/*
 * CustomGridGen.hpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

#ifndef __GRID_GENERATOR_HPP__
#define __GRID_GENERATOR_HPP__

#include <deal.II/grid/tria.h>

using namespace dealii;

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/intergrid_map.h>

#include <deal.II/distributed/tria.h>

#include <iostream>
#include <cmath>
#include <limits>

#include "errors.hpp"
#include "Utils.hpp"

template <int dim>
class MyGridGenerator {

public:

	static void hyper_rectangle(dealii::Triangulation<dim> &triangulation,
			double length, double width);

	/**
	 * This method builds a finite element domain of the type:
	 *
	 *
	 *    inter_hole_space
	 *        <---->       hole_length
	 *  __     ____     ____<--->__        _    _
	 * |  |___|    |___|    |___|  |       |    _ -> hole_width
	 * .                           .       |
	 * .                           .     width
	 * .                           .       |
	 * |___________________________|       _
	 *
	 * <--------- length ---------->
	 *
	 * In this example, the holes_nbr is 3.
	 * Note that an half inter_hole_space is set near the left and right sides of
	 * the rectangle.
	 *
	 * @pre:
	 * - tria is an instantiated object.
	 * - 0 <= width
	 * - 0 <= hole_length
	 * - 0 <= hole_width <= width
	 * - 0 <= inter_hole_space
	 *
	 * if pre not respected (except for the tria one), the exception
	 * PRECONDITIONS_VIOLATED is thrown.
	 *
	 */
	//TODO generalize to 3D
	static void serrated_hyper_rectangle(dealii::Triangulation<dim> &tria,
				double width, unsigned holes_nbr, double hole_length,
				double hole_width, double inter_hole_space);

private:

	static void gen_points_for_serrated_hrect(double length, double width,
			double hole_length,	double hole_width, double inter_hole_space,
			unsigned nbr_of_cells_fst_dim, unsigned nbr_of_cells_scd_dim,
			std::vector<Point<2> > &points);

	static void gen_cells_for_serrated_hrect(unsigned nbr_of_cells_fst_dim,
			unsigned nbr_of_cells_scd_dim, std::vector<CellData<2> > &cells);

	static bool are_precond_fullfilled_serr_hrect(double width,
			double hole_length, double hole_width, double inter_hole_space);
};

#include "MyGridGenerator.cpp"

#endif
