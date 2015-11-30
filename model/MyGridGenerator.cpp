/*
 * GridGenerator.cpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

template<int dim>
void MyGridGenerator<dim>::hyper_rectangle(
		dealii::Triangulation<dim> &triangulation, double length,
		double width) {
	Point < dim > point_bottom(0, 0);
	Point < dim > point_top(length, width);
	dealii::GridGenerator::hyper_rectangle(triangulation, point_bottom,
			point_top);
}

template<int dim>
void MyGridGenerator<dim>::gen_points_for_serrated_hrect(double length,
		double width, double hole_length, double hole_width,
		double inter_hole_space, unsigned nbr_of_cells_fst_dim,
		unsigned nbr_of_cells_scd_dim, std::vector<Point<2> > &points) {

	//End of half dent on left detector side
	double left_half_dent_end_x = inter_hole_space / 2.0;

	double y = 0.0;
	for (unsigned j = 0; j <= nbr_of_cells_scd_dim; j++) {

		Point < 2 > left_corner_point(0.0, y);
		points.push_back(left_corner_point);

		bool is_hole_left_point = true;

		double x = left_half_dent_end_x;
		for (unsigned i = 0; i < nbr_of_cells_fst_dim - 1; i++) {

			Point < 2 > new_point(x, y);
			points.push_back(new_point);

			if (is_hole_left_point)
				x += hole_length;
			else
				x += inter_hole_space;

			is_hole_left_point = !is_hole_left_point;
		}

		Point < 2 > end_point(length, y);
		points.push_back(end_point);

		if(j == 0)
			y = width-(nbr_of_cells_scd_dim-1)*hole_width;
		else
			y += hole_width;
	}
}

template<int dim>
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
template<int dim>
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
template<int dim>
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

	if (length == 0.0 || width == 0.0 || hole_width == 0.0
			|| hole_length == 0.0 || holes_nbr == 0)
		return hyper_rectangle(tria, length, width);

	if (inter_hole_space == 0.0)
		return hyper_rectangle(tria, length, width - hole_width);

	unsigned nbr_of_cells_fst_dim = holes_nbr * 2 + 1;

	double epsilon = 0.000001;

	unsigned nbr_of_cells_scd_dim =	ceil(width / hole_width);

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
