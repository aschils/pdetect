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

	static void rectangle(dealii::Triangulation<dim> &triangulation,
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
	 * - hole_width <= width
	 * - 0 <= half_inter_hole_space
	 *
	 * if pre not respected (except for the tria one), the exception
	 * PRECONDITIONS_VIOLATED is throwned.
	 *
	 */
	static void serrated_rectangle(Triangulation<dim> &tria,
			unsigned width, unsigned holes_nbr, unsigned hole_length,
			unsigned hole_width, unsigned half_inter_hole_space);

	static void serrated_rectangle_gaussian(Triangulation<dim> &tria,
			unsigned width, unsigned holes_nbr,
			unsigned inter_gauss_borders_lgth, unsigned gaussian_peak);

	static void rectangle_with_circular_holes(Triangulation<dim> &tria,
			unsigned half_width, unsigned holes_radius,
			unsigned half_inter_holes_centers_dist, int nbr_of_holes);

	static void rectangle_with_rectangular_holes(
			Triangulation<dim> &tria, unsigned half_width,
			unsigned hole_length, unsigned half_hole_width,
			unsigned half_inter_holes_dist, unsigned nbr_of_holes);

	static void rectangle_with_circular_hole(Triangulation<dim> &tria,
				unsigned half_width, unsigned half_length, unsigned holes_radius);

private:

	static void gaussian_hole(Triangulation<dim> &tria,
			unsigned width,	unsigned inter_gauss_borders_lgth,
			unsigned gaussian_peak, double standard_deviation,
			unsigned gaussian_sampling_lvl);
};

template<unsigned dim>
void MyGridGenerator<dim>::rectangle(
		dealii::Triangulation<dim> &triangulation, double length,
		double width) {
	Point<dim> point_bottom(0, 0);
	Point<dim> point_top(length, width);
	dealii::GridGenerator::hyper_rectangle(triangulation, point_bottom,
			point_top);
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
 * - 0 <= half_inter_hole_space
 *
 * if pre not respected (except for the tria one), the exception
 * PRECONDITIONS_VIOLATED is throwned.
 *
 */
template<unsigned dim>
void MyGridGenerator<dim>::serrated_rectangle(Triangulation<dim> &tria,
		unsigned width, unsigned holes_nbr, unsigned hole_length,
		unsigned hole_width, unsigned half_inter_hole_space) {

	unsigned inter_hole_space = 2 * half_inter_hole_space;

	if (hole_width > width)
		throw PRECONDITIONS_VIOLATED;

	unsigned length =
			(holes_nbr == 0) ?
					inter_hole_space :
					holes_nbr * (hole_length + inter_hole_space);

	if (length == 0 || width == 0 || hole_width == 0 || hole_length == 0
			|| holes_nbr == 0)
		return rectangle(tria, length, width);

	if (inter_hole_space == 0)
		return rectangle(tria, length, width - hole_width);

	Triangulation<dim> half_inter_hole, inter_hole, half_rect, hole_width_rect,
			inter_hole_space_width_rect;

	//Build half_rect
	Point<dim> bottom_left_half_rect_pt(0.0, 0.0);
	Point<dim> top_right_half_rect_pt(half_inter_hole_space,
			width - hole_width);
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
	GridTools::shift(half_inter_hole_space_translation,
			inter_hole_space_width_rect);

	dealii::Tensor<1, dim> hole_translation;
	hole_translation[0] = hole_length;
	hole_translation[1] = 0;
	dealii::Tensor<1, dim> inter_hole_translation;
	inter_hole_translation[0] = inter_hole_space;
	inter_hole_translation[1] = 0;

	//Add "middle" rectangles and holes to tria
	for (unsigned i = 0; i < 2 * holes_nbr - 1; i++) {

		if (i % 2 == 0) {
			GridGenerator::merge_triangulations(tria, hole_width_rect, tria);
			GridTools::shift(hole_translation, hole_width_rect);
			GridTools::shift(hole_translation, inter_hole);
			GridTools::shift(hole_translation, inter_hole_space_width_rect);
		} else {
			GridGenerator::merge_triangulations(tria,
					inter_hole_space_width_rect, tria);
			GridGenerator::merge_triangulations(tria, inter_hole, tria);
			GridTools::shift(inter_hole_translation, hole_width_rect);
			GridTools::shift(inter_hole_translation, inter_hole);
			GridTools::shift(inter_hole_translation,
					inter_hole_space_width_rect);
		}
	}

	//Add Half-hole right extremity half rect and half hole
	hole_translation[0] = half_inter_hole_space + hole_length * holes_nbr
			+ (holes_nbr - 1) * inter_hole_space;
	GridTools::shift(hole_translation, half_inter_hole);
	GridTools::shift(hole_translation, half_rect);
	GridGenerator::merge_triangulations(tria, half_inter_hole, tria);
	GridGenerator::merge_triangulations(tria, half_rect, tria);
}

template <unsigned dim>
void MyGridGenerator<dim>::gaussian_hole(Triangulation<dim> &tria,
		unsigned width,	unsigned inter_gauss_borders_lgth,
		unsigned gaussian_peak, double standard_deviation,
		unsigned gaussian_sampling_lvl){

	double max_error = gaussian_peak/((double)gaussian_sampling_lvl);
	double start_sampling = 0.0;
	double end_sampling = Utils::gaussian_near_zero_x(gaussian_peak,
			standard_deviation, 0, max_error);
	double interval = end_sampling - start_sampling;

	unsigned nbr_of_gauss_pts = gaussian_sampling_lvl+2;
	std::vector<double> gaussian_x(nbr_of_gauss_pts);
	std::vector<double> gaussian_y(nbr_of_gauss_pts);

	double x = 0.0;
	double step = interval/(nbr_of_gauss_pts - 1);

	for(unsigned i=0; i<nbr_of_gauss_pts; i++){
		gaussian_x[i] = x;
		gaussian_y[i] = Utils::gaussian(gaussian_peak, standard_deviation,
				0.0, x);
		x += step;
	}

	//Now we have the right part of a gaussian. We want to inverse it regarding
	// the vertical axis (thus inverse x coord of each point, x -> -x) to get
	// the left part. Then inverse it regarding horizontal axis (thus invers
	// y coord) to get left border of the gaussian hole.
	TensorUtils::opposite_vector(gaussian_x);
	TensorUtils::opposite_vector(gaussian_y);

	std::vector< Point< dim > > vertices;
	std::vector<CellData<dim> > cells;

		Point<dim> p1(0,0);
		Point<dim> p2(0.2,1);
		Point<dim> p3(2,0);
		Point<dim> p4(2.2,1);

		vertices.push_back(p1);
		vertices.push_back(p2);
		vertices.push_back(p3);
		vertices.push_back(p4);

		CellData<dim> cell;
		for(unsigned i=0; i<4;i++)
			cell.vertices[i] = i;
		cell.material_id = 0;

		cells.push_back(cell);

		tria.create_triangulation(vertices, cells, SubCellData());

	Utils::gaussian(gaussian_peak, standard_deviation, 0, 0);
}

template <unsigned dim>
void MyGridGenerator<dim>::serrated_rectangle_gaussian(
		Triangulation<dim> &tria,
		unsigned width, unsigned holes_nbr,
		unsigned inter_gauss_borders_lgth, unsigned gaussian_peak){

	std::vector< Point< dim > > vertices;
	std::vector<CellData<dim> > cells;

	Point<dim> p1(0,0);
	Point<dim> p2(0.2,1);
	Point<dim> p3(2,0);
	Point<dim> p4(2.2,1);

	vertices.push_back(p1);
	vertices.push_back(p2);
	vertices.push_back(p3);
	vertices.push_back(p4);

	CellData<dim> cell;
	for(unsigned i=0; i<4;i++)
		cell.vertices[i] = i;
	cell.material_id = 0;

	cells.push_back(cell);

	tria.create_triangulation(vertices, cells, SubCellData());
}

template<unsigned dim>
void MyGridGenerator<dim>::rectangle_with_circular_hole(
		dealii::Triangulation<dim> &tria, unsigned half_width_p,
		unsigned half_length_p, unsigned holes_radius) {

	int half_width = half_width_p;
	int half_length = half_length_p;
	double h_shell_outer_radius = std::min<int>(half_width, half_length);

	Assert(holes_radius < h_shell_outer_radius,
			ExcMessage("outer_radius has to be bigger than inner_radius."));

	Point<dim> center(0.0,0.0);
	unsigned h_shell_nbr_of_cells = 8;

	// We create an hyper_shell in two dimensions, and then we modify it.
	dealii::GridGenerator::hyper_shell(tria, center, holes_radius,
			h_shell_outer_radius, h_shell_nbr_of_cells);

	typename dealii::Triangulation<dim>::active_cell_iterator cell =
			tria.begin_active();
	typename dealii::Triangulation<dim>::active_cell_iterator endc = tria.end();

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
					if (!treated_vertices[vv]) {
						treated_vertices[vv] = true;

						//TODO cleaner solution?

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
		dealii::Triangulation<dim> &tria, unsigned half_width,
		unsigned holes_radius, unsigned half_inter_holes_centers_dist,
		int nbr_of_holes) {

	unsigned width = 2*half_width;
	unsigned inter_holes_centers_dist = 2*half_inter_holes_centers_dist;

	if (nbr_of_holes < 1 || 2 * holes_radius >= width || width == 0
			|| holes_radius == 0 || inter_holes_centers_dist == 0
			|| inter_holes_centers_dist <= 2 * holes_radius)
		throw INVALID_INPUT;

	dealii::Tensor<1, dim> half_length_translation;
	half_length_translation[0] = half_inter_holes_centers_dist;

	if (nbr_of_holes == 1) {
		MyGridGenerator<dim>::rectangle_with_circular_hole(tria, half_width,
				half_inter_holes_centers_dist, holes_radius);
		GridTools::shift(half_length_translation, tria);
		return;
	}

	dealii::Triangulation<dim> periodic_struct;
	/*
	 * TODO optimize
	 * Note: I tried to call the function copy_triangulation of
	 * dealii::Triangulation class instead of creating two times the mesh,
	 * but it segfaulted...
	 */
	MyGridGenerator<dim>::rectangle_with_circular_hole(periodic_struct, half_width,
			half_inter_holes_centers_dist, holes_radius);
	MyGridGenerator<dim>::rectangle_with_circular_hole(tria, half_width,
			half_inter_holes_centers_dist, holes_radius);

	dealii::Tensor<1, dim> translation;
	translation[0] = inter_holes_centers_dist;

	for (unsigned i = 1; i < nbr_of_holes; i++) {
		GridTools::shift(translation, periodic_struct);
		GridGenerator::merge_triangulations(tria, periodic_struct, tria);
	}

	GridTools::shift(half_length_translation, tria);
}

/**
 * Build a grid of the form:
 *       _                     _  ________________________________________________________
 *       | half_wo_hole_width  |  |  top_left  | top_around_hole  |  top_between_holes    | ....
 *       |                     -  |--------------------------------------------------------
 * width | hole_width          |  | middle_left|   HOLE           |  middle_between_holes | HOLE
 *       |                     -  |-------------------------------------------------------
 *       | half_wo_hole_width  |  | bottom_left|bottom_around_hole| bottom_between_holes  | ...
 *       _                     _  |________________________________________________________
 *
 *                                |------------|------------------|-----------------------|
 *
 *                       half_inter_holes_dist     hole_length         inter_holes_dist
 *
 * note that {top,middle,bottom}_{left,right} rectangles have
 * length = ceil(inter_holes_dist/2.0)
 *
 */
template<unsigned dim>
void MyGridGenerator<dim>::rectangle_with_rectangular_holes(
		dealii::Triangulation<dim> &tria, unsigned half_width,
		unsigned hole_length, unsigned half_hole_width,
		unsigned half_inter_holes_dist, unsigned nbr_of_holes) {

	unsigned inter_holes_dist = 2*half_inter_holes_dist;
	unsigned width = 2*half_width;
	unsigned hole_width = 2*half_hole_width;

	if (!(hole_width > 0 && width > hole_width) || nbr_of_holes == 0
			|| hole_length == 0)
		throw INVALID_INPUT;

	bool no_inter_holes_dist = inter_holes_dist == 0;

	//Build the various rectangles needed before assembling them
	Triangulation<dim> bottom_left, middle_left, top_left, bottom_around_hole,
			top_around_hole, bottom_between_holes, middle_between_holes,
			top_between_holes;

	unsigned half_wo_hole_width = half_width - half_hole_width;

	//Build bottom left
	Point<dim> p1(0.0, 0.0);
	Point<dim> p2(half_inter_holes_dist, half_wo_hole_width);
	GridGenerator::hyper_rectangle(bottom_left, p1, p2);

	//Define useful vertical translation vectors
	dealii::Tensor<1, dim> half_wo_hole_width_translation;
	half_wo_hole_width_translation[0] = 0;
	half_wo_hole_width_translation[1] = half_wo_hole_width;
	dealii::Tensor<1, dim> hole_width_translation;
	hole_width_translation[0] = 0;
	hole_width_translation[1] = hole_width;
	dealii::Tensor<1, dim> bottom_to_top_translation =
			half_wo_hole_width_translation + hole_width_translation;

	//Build top left
	GridGenerator::hyper_rectangle(top_left, p1, p2);
	GridTools::shift(bottom_to_top_translation, top_left);

	//Build middle_left
	Point<dim> p3(0.0, half_wo_hole_width);
	Point<dim> p4(half_inter_holes_dist, hole_width + half_wo_hole_width);
	GridGenerator::hyper_rectangle(middle_left, p3, p4);

	//Build bottom_around_hole
	Point<dim> p5(half_inter_holes_dist, 0.0);
	Point<dim> p6(half_inter_holes_dist + hole_length, half_wo_hole_width);
	GridGenerator::hyper_rectangle(bottom_around_hole, p5, p6);

	//Build top_around_hole
	GridGenerator::hyper_rectangle(top_around_hole, p5, p6);
	GridTools::shift(bottom_to_top_translation, top_around_hole);

	if (!no_inter_holes_dist) {
		//Build bottom_between_holes
		Point<dim> p7(half_inter_holes_dist + hole_length, 0.0);
		Point<dim> p8(half_inter_holes_dist + hole_length + inter_holes_dist,
				half_wo_hole_width);
		GridGenerator::hyper_rectangle(bottom_between_holes, p7, p8);

		//Build top_between_holes
		GridGenerator::hyper_rectangle(top_between_holes, p7, p8);
		GridTools::shift(bottom_to_top_translation, top_between_holes);

		//Build middle_between_holes
		Point<dim> p9(half_inter_holes_dist + hole_length, half_wo_hole_width);
		Point<dim> p10(half_inter_holes_dist + hole_length + inter_holes_dist,
				half_wo_hole_width + hole_width);
		GridGenerator::hyper_rectangle(middle_between_holes, p9, p10);
	}

	//Put left rectangles in tria
	GridGenerator::hyper_rectangle(tria, p1, p2); //init tria to bottom_left
	GridGenerator::merge_triangulations(tria, middle_left, tria);
	GridGenerator::merge_triangulations(tria, top_left, tria);

	//Define useful horizontal translation vectors
	dealii::Tensor<1, dim> hole_length_translation;
	hole_length_translation[0] = hole_length;
	dealii::Tensor<1, dim> inter_holes_dist_translation;
	inter_holes_dist_translation[0] = inter_holes_dist;
	dealii::Tensor<1, dim> periodic_struct_translation = hole_length_translation
			+ inter_holes_dist_translation;

	//Now put the "middle rectangles" in tria
	unsigned nbr_of_rect_to_put_middle = nbr_of_holes * 2 - 1;
	for (unsigned i = 0; i < nbr_of_rect_to_put_middle; i++) {

		if (i % 2 == 0) {
			GridGenerator::merge_triangulations(tria, bottom_around_hole, tria);
			GridGenerator::merge_triangulations(tria, top_around_hole, tria);
			GridTools::shift(periodic_struct_translation, bottom_around_hole);
			GridTools::shift(periodic_struct_translation, top_around_hole);
		}

		else if (!no_inter_holes_dist) {
			GridGenerator::merge_triangulations(tria, bottom_between_holes,
					tria);
			GridGenerator::merge_triangulations(tria, middle_between_holes,
					tria);
			GridGenerator::merge_triangulations(tria, top_between_holes, tria);
			GridTools::shift(periodic_struct_translation, bottom_between_holes);
			GridTools::shift(periodic_struct_translation, middle_between_holes);
			GridTools::shift(periodic_struct_translation, top_between_holes);
		}
	}

	//Build right rectangles by shifting left rectangles
	Triangulation<dim> bottom_right, middle_right, top_right;
	unsigned to_shift = half_inter_holes_dist + nbr_of_holes * hole_length
			+ (nbr_of_holes - 1) * inter_holes_dist;
	dealii::Tensor<1, dim> left_rect_translation;
	left_rect_translation[0] = to_shift;
	GridTools::shift(left_rect_translation, top_left);
	GridTools::shift(left_rect_translation, middle_left);
	GridTools::shift(left_rect_translation, bottom_left);

	//Put right rectangles in tria
	GridGenerator::merge_triangulations(tria, bottom_left, tria);
	GridGenerator::merge_triangulations(tria, middle_left, tria);
	GridGenerator::merge_triangulations(tria, top_left, tria);
}
