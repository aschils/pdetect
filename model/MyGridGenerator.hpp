/*
 * CustomGridGen.hpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

#pragma once

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

template<unsigned dim>
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
	static void serrated_hyper_rectangle(Triangulation<dim> &tria,
			unsigned width, unsigned holes_nbr, unsigned hole_length, unsigned hole_width,
			unsigned inter_hole_space);

	static void rectangle_with_circular_holes(dealii::Triangulation<dim> &tria,
			unsigned width, unsigned holes_radius,
			unsigned inter_holes_centers_dist, int nbr_of_holes);

private:

	static bool are_precond_fullfilled_serr_hrect(unsigned width,
			unsigned hole_length, unsigned hole_width, unsigned inter_hole_space);

	static void rectangle_with_circular_hole(dealii::Triangulation<dim> &tria,
			unsigned width, unsigned length, unsigned holes_radius);
};

template<unsigned dim>
void MyGridGenerator<dim>::hyper_rectangle(
		dealii::Triangulation<dim> &triangulation, double length,
		double width) {
	Point<dim> point_bottom(0, 0);
	Point<dim> point_top(length, width);
	dealii::GridGenerator::hyper_rectangle(triangulation, point_bottom,
			point_top);
}

/**
 * Verify the following conditions on parameters:
 * - hole_width <= width
 */
template<unsigned dim>
bool MyGridGenerator<dim>::are_precond_fullfilled_serr_hrect(unsigned width,
		unsigned hole_length, unsigned hole_width, unsigned inter_hole_space) {

	//return width >= 0.0 && hole_length >= 0 && hole_width >= 0
	//		&& hole_width <= width && inter_hole_space >= 0.0;
	return hole_width <= width;
}

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
 * - hole_width <= width
 * - 0 <= inter_hole_space
 *
 * if pre not respected (except for the tria one), the exception
 * PRECONDITIONS_VIOLATED is throwned.
 *
 */
template<unsigned dim>
void MyGridGenerator<dim>::serrated_hyper_rectangle(Triangulation<dim> &tria,
		unsigned width, unsigned holes_nbr, unsigned hole_length, unsigned hole_width,
		unsigned inter_hole_space) {

	if (!are_precond_fullfilled_serr_hrect(width, hole_length, hole_width,
			inter_hole_space))
		throw PRECONDITIONS_VIOLATED;

	unsigned length =
			(holes_nbr == 0) ?
					inter_hole_space :
					holes_nbr * (hole_length + inter_hole_space);

	if (length == 0 || width == 0 || hole_width == 0 || hole_length == 0
			|| holes_nbr == 0)
		return hyper_rectangle(tria, length, width);

	if (inter_hole_space == 0)
		return hyper_rectangle(tria, length, width - hole_width);

	unsigned half_inter_hole_space = ceil(inter_hole_space / 2.0);

	Triangulation<dim> half_inter_hole, inter_hole, half_rect, hole_width_rect,
	inter_hole_space_width_rect;

	//Build half_rect
	Point<dim> bottom_left_half_rect_pt(0.0, 0.0);
	Point<dim> top_right_half_rect_pt(half_inter_hole_space, width - hole_width);
	GridGenerator::hyper_rectangle(half_rect, bottom_left_half_rect_pt,
			top_right_half_rect_pt);

	//Build half inter hole for rectangle extremities
	Point<dim> bottom_left_half_hole_pt(0, width - hole_width);
	Point<dim> top_right_half_hole_pt(half_inter_hole_space, width);
	GridGenerator::hyper_rectangle(half_inter_hole, bottom_left_half_hole_pt,
			top_right_half_hole_pt);

	//Build rects
	Point<dim> bottom_left_rect_pt(0.0, 0.0);
	Point<dim> top_right_rect_pt(hole_length, width - hole_width);
	GridGenerator::hyper_rectangle(hole_width_rect, bottom_left_rect_pt,
			top_right_rect_pt);

	Point<dim> bottom_left_rect_pt2(0.0, 0.0);
	Point<dim> top_right_rect_pt2(inter_hole_space, width - hole_width);
	GridGenerator::hyper_rectangle(inter_hole_space_width_rect,
			bottom_left_rect_pt2, top_right_rect_pt2);

	//Build inter hole
	Point<dim> bottom_left_hole_pt(0, width - hole_width);
	Point<dim> top_right_hole_pt(inter_hole_space, width);
	GridGenerator::hyper_rectangle(inter_hole, bottom_left_hole_pt,
			top_right_hole_pt);

	//Init tria to half_rect + half non hole left
	GridGenerator::hyper_rectangle(tria, bottom_left_half_rect_pt,
			top_right_half_rect_pt);
	GridGenerator::merge_triangulations(tria, half_inter_hole, tria);

	dealii::Tensor<1, dim> half_inter_hole_space_translation;
	half_inter_hole_space_translation[0] = half_inter_hole_space;
	GridTools::shift(half_inter_hole_space_translation, hole_width_rect);
	GridTools::shift(half_inter_hole_space_translation, inter_hole);
	GridTools::shift(half_inter_hole_space_translation, inter_hole_space_width_rect);

	dealii::Tensor<1, dim> hole_translation;
	hole_translation[0] = hole_length;
	hole_translation[1] = 0;
	dealii::Tensor<1, dim> inter_hole_translation;
	inter_hole_translation[0] = inter_hole_space;
	inter_hole_translation[1] = 0;

	//Add "middle" rectangles and holes to tria
	for(unsigned i=0; i< 2*holes_nbr-1; i++){

		if(i%2 == 0){
			GridGenerator::merge_triangulations(tria, hole_width_rect, tria);
			GridTools::shift(hole_translation, hole_width_rect);
			GridTools::shift(hole_translation, inter_hole);
			GridTools::shift(hole_translation, inter_hole_space_width_rect);
		}
		else{
			GridGenerator::merge_triangulations(tria, inter_hole_space_width_rect, tria);
			GridGenerator::merge_triangulations(tria, inter_hole, tria);
			GridTools::shift(inter_hole_translation, hole_width_rect);
			GridTools::shift(inter_hole_translation, inter_hole);
			GridTools::shift(inter_hole_translation, inter_hole_space_width_rect);
		}
	}

	//Add Half-hole right extremity half rect and half hole
	hole_translation[0] = half_inter_hole_space + hole_length *holes_nbr+
			(holes_nbr-1)*inter_hole_space;
	GridTools::shift(hole_translation, half_inter_hole);
	GridTools::shift(hole_translation, half_rect);
	GridGenerator::merge_triangulations(tria, half_inter_hole, tria);
	GridGenerator::merge_triangulations(tria, half_rect, tria);
}

template<unsigned dim>
void MyGridGenerator<dim>::rectangle_with_circular_hole(
		dealii::Triangulation<dim> &tria, unsigned width, unsigned length,
		unsigned holes_radius) {

	unsigned shorter_side = std::min<unsigned>(width, length);
	double h_shell_outer_radius = ceil(shorter_side / 2.0);

	Assert(inner_radius < h_shell_outer_radius,
			ExcMessage("outer_radius has to be bigger than inner_radius."));

	Point<dim> center;
	unsigned h_shell_nbr_of_cells = 8;

	// We create an hyper_shell in two dimensions, and then we modify it.
	dealii::GridGenerator::hyper_shell(tria, center, holes_radius,
			h_shell_outer_radius, h_shell_nbr_of_cells);

	typename dealii::Triangulation<dim>::active_cell_iterator cell =
			tria.begin_active();
	typename dealii::Triangulation<dim>::active_cell_iterator endc = tria.end();

	int half_length = ceil(length / 2.0);
	int half_width = ceil(width / 2.0);

	/*
	 * What we obtain after transforming hypershell with the following code:
	 *
	 *  3         2        1
	 *  ___________________
	 *  |                  |
	 * 4|   circle here    | 0
	 *  |__________________|
	 *  5         6        7
	 *
	 * The numbers are the vertices indices (helps to understand the code
	 * in the switch below). Before transformation we had these vertices
	 * organized to form an outer circle, thus we have transformed this
	 * outer circle in a rectangle.
	 *
	 */

	std::vector<bool> treated_vertices(tria.n_vertices(), false);
	for (; cell != endc; ++cell) {
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
			if (cell->face(f)->at_boundary()) {
				for (unsigned int v = 0;
						v < GeometryInfo<dim>::vertices_per_face; ++v) {
					unsigned int vv = cell->face(f)->vertex_index(v);
					if (treated_vertices[vv] == false) {
						treated_vertices[vv] = true;

						//TODO cleaner solution

						switch (vv) {

						case 0:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(half_length, 0);
							break;
						case 1:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(half_length, half_width);
							break;
						case 2:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(0, half_width);
							break;
						case 3:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(-half_length, half_width);
							break;
						case 4:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(-half_length, 0);
							break;
						case 5:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(-half_length, -half_width);
							break;
						case 6:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(0, -half_width);
							break;
						case 7:
							cell->face(f)->vertex(v) = center
									+ Point<dim>(half_length, -half_width);
							break;
						default:
							break;
						}
					}
				}
			}
	}
}

template<unsigned dim>
void MyGridGenerator<dim>::rectangle_with_circular_holes(
		dealii::Triangulation<dim> &tria, unsigned width, unsigned holes_radius,
		unsigned inter_holes_centers_dist, int nbr_of_holes) {

	if (nbr_of_holes < 1 || 2 * holes_radius >= width || width == 0
			|| holes_radius == 0 || inter_holes_centers_dist == 0
			|| inter_holes_centers_dist <= 2 * holes_radius)
		throw INVALID_INPUT;

	if (nbr_of_holes == 1) {
		MyGridGenerator<dim>::rectangle_with_circular_hole(tria, width,
				inter_holes_centers_dist, holes_radius);
		return;
	}

	dealii::Triangulation<dim> periodic_struct;
	/*
	 * TODO optimize
	 * Note: I tried to call the function copy_triangulation of
	 * dealii::Triangulation class instead of creating two times the mesh,
	 * but it segfaults...
	 */
	MyGridGenerator<dim>::rectangle_with_circular_hole(periodic_struct, width,
			inter_holes_centers_dist, holes_radius);
	MyGridGenerator<dim>::rectangle_with_circular_hole(tria, width,
			inter_holes_centers_dist, holes_radius);

	dealii::Tensor<1, dim> translation;
	translation[0] = inter_holes_centers_dist;

	for (unsigned i = 1; i < nbr_of_holes; i++) {
		GridTools::shift(translation, periodic_struct);
		GridGenerator::merge_triangulations(tria, periodic_struct, tria);
	}
}
