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
	static void serrated_hyper_rectangle(dealii::Triangulation<dim> &tria,
			double width, unsigned holes_nbr, double hole_length,
			double hole_width, double inter_hole_space);

	static void rectangle_with_circular_holes(dealii::Triangulation<dim> &tria,
			unsigned width, unsigned holes_radius,
			unsigned inter_holes_centers_dist, int nbr_of_holes);

private:

	static void gen_points_for_serrated_hrect(double length, double width,
			double hole_length, double hole_width, double inter_hole_space,
			unsigned nbr_of_cells_fst_dim, unsigned nbr_of_cells_scd_dim,
			std::vector<Point<2> > &points);

	static void gen_cells_for_serrated_hrect(unsigned nbr_of_cells_fst_dim,
			unsigned nbr_of_cells_scd_dim, std::vector<CellData<2> > &cells);

	static bool are_precond_fullfilled_serr_hrect(double width,
			double hole_length, double hole_width, double inter_hole_space);

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

template<unsigned dim>
void MyGridGenerator<dim>::gen_points_for_serrated_hrect(double length,
		double width, double hole_length, double hole_width,
		double inter_hole_space, unsigned nbr_of_cells_fst_dim,
		unsigned nbr_of_cells_scd_dim, std::vector<Point<2> > &points) {

	//End of half dent on left detector side
	double left_half_dent_end_x = inter_hole_space / 2.0;

	double y = 0.0;
	for (unsigned j = 0; j <= nbr_of_cells_scd_dim; j++) {

		Point<2> left_corner_point(0.0, y);
		points.push_back(left_corner_point);

		bool is_hole_left_point = true;

		double x = left_half_dent_end_x;
		for (unsigned i = 0; i < nbr_of_cells_fst_dim - 1; i++) {

			Point<2> new_point(x, y);
			points.push_back(new_point);

			if (is_hole_left_point)
				x += hole_length;
			else
				x += inter_hole_space;

			is_hole_left_point = !is_hole_left_point;
		}

		Point<2> end_point(length, y);
		points.push_back(end_point);

		if (j == 0)
			y = width - (nbr_of_cells_scd_dim - 1) * hole_width;
		else
			y += hole_width;
	}
}

template<unsigned dim>
void MyGridGenerator<dim>::gen_cells_for_serrated_hrect(
		unsigned nbr_of_cells_fst_dim, unsigned nbr_of_cells_scd_dim,
		std::vector<CellData<2> > &cells) {

	unsigned cell_idx = 0;
	for (unsigned y = 0; y < nbr_of_cells_scd_dim; ++y) {
		for (unsigned x = 0; x < nbr_of_cells_fst_dim; ++x) {

			bool is_top_width_cell = (y == nbr_of_cells_scd_dim - 1);
			bool is_hole_cell = (x % 2 == 1);

			if (is_top_width_cell && is_hole_cell)
				continue;

			Assert(cell_idx < cells.size(), ExcInternalError());

			unsigned top_vertices_row = y * (nbr_of_cells_fst_dim + 1);
			unsigned top_left_vertice_idx = top_vertices_row + x;
			unsigned top_right_vertice_idx = top_left_vertice_idx + 1;

			unsigned bottom_vertices_row = top_vertices_row
					+ nbr_of_cells_fst_dim + 1;
			unsigned bottom_left_vertice_idx = bottom_vertices_row + x;
			unsigned bottom_right_vertice_idx = bottom_left_vertice_idx + 1;

			cells[cell_idx].vertices[0] = top_left_vertice_idx;
			cells[cell_idx].vertices[1] = top_right_vertice_idx;
			cells[cell_idx].vertices[2] = bottom_left_vertice_idx;
			cells[cell_idx].vertices[3] = bottom_right_vertice_idx;
			cells[cell_idx].material_id = 0;

			cell_idx++;
		}
	}
}

/**
 * Verify the following conditions on parameters:
 * - 0 <= width
 * - 0 <= hole_length
 * - 0 <= hole_width <= width
 * - 0 <= inter_hole_space
 */
template<unsigned dim>
bool MyGridGenerator<dim>::are_precond_fullfilled_serr_hrect(double width,
		double hole_length, double hole_width, double inter_hole_space) {

	return width >= 0.0 && hole_length >= 0 && hole_width >= 0
			&& hole_width <= width && inter_hole_space >= 0.0;
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
 * - 0 <= hole_width <= width
 * - 0 <= inter_hole_space
 *
 * if pre not respected (except for the tria one), the exception
 * PRECONDITIONS_VIOLATED is throwned.
 *
 */
template<unsigned dim>
void MyGridGenerator<dim>::serrated_hyper_rectangle(Triangulation<dim> &tria,
		double width, unsigned holes_nbr, double hole_length, double hole_width,
		double inter_hole_space) {

	if (!are_precond_fullfilled_serr_hrect(width, hole_length, hole_width,
			inter_hole_space))
		throw PRECONDITIONS_VIOLATED;

	double length =
			(holes_nbr == 0) ?
					inter_hole_space :
					holes_nbr * (hole_length + inter_hole_space);

	if (length == 0.0 || width == 0.0 || hole_width == 0.0 || hole_length == 0.0
			|| holes_nbr == 0)
		return hyper_rectangle(tria, length, width);

	if (inter_hole_space == 0.0)
		return hyper_rectangle(tria, length, width - hole_width);

	unsigned nbr_of_cells_fst_dim = holes_nbr * 2 + 1;

	unsigned nbr_of_cells_scd_dim = ceil(width / hole_width);

	std::vector<Point<2> > points;
	gen_points_for_serrated_hrect(length, width, hole_length, hole_width,
			inter_hole_space, nbr_of_cells_fst_dim, nbr_of_cells_scd_dim,
			points);

	std::vector<CellData<2> > cells;
	cells.resize(nbr_of_cells_fst_dim * nbr_of_cells_scd_dim - holes_nbr);
	gen_cells_for_serrated_hrect(nbr_of_cells_fst_dim, nbr_of_cells_scd_dim,
			cells);

	tria.create_triangulation(points, cells, SubCellData());
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

	if(nbr_of_holes < 1 || 2*holes_radius >= width || width == 0 ||
			holes_radius == 0 || inter_holes_centers_dist == 0
			|| inter_holes_centers_dist <= 2*holes_radius)
		throw INVALID_INPUT_EXCEPTION;

	if(nbr_of_holes == 1){
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
