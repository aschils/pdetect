/*
 * CustomGridGen.hpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

#ifndef __CUSTOM_GRID_GEN_HPP__
#define __CUSTOM_GRID_GEN_HPP__

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

using namespace dealii;

void serrated_rectangle(Triangulation<2> &tria,
		unsigned nbr_of_strips);

#endif
