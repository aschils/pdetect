/*
 * GridGenerator.cpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

template<int dim>
void MyGridGenerator<dim>::hyper_rectangle(dealii::Triangulation<dim> &triangulation,
		double length, double width) {
	Point < dim > point_bottom(-length / 2, -width / 2);
	Point < dim > point_top(length / 2, width / 2);
	dealii::GridGenerator::hyper_rectangle(triangulation, point_bottom, point_top);
}

template<int dim>
void MyGridGenerator<dim>::gen_points_for_serrated_hrect(double length,
		double hole_length, double hole_width, double inter_hole_space,
		unsigned nbr_of_cells_fst_dim, unsigned nbr_of_cells_scd_dim,
		std::vector<Point<2> > &points) {

	//End of half dent on left detector side
	double left_half_dent_end_x = inter_hole_space / 2.0;
	double right_half_dent_start_x = length - left_half_dent_end_x;

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
		y += hole_width;
	}
}

template<int dim>
void MyGridGenerator<dim>::gen_cells_for_serrated_hrect(
		unsigned nbr_of_cells_fst_dim, unsigned nbr_of_cells_scd_dim,
		std::vector<CellData<2> > &cells) {

	std::cout << "enters gen_cells " << std::endl;

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
 * @pre:
 * - length > 0
 * - width > 0
 * - 0 <= hole_length <= length
 * - 0 <= hole_width <= width
 * - 0 <= inter_hole_space <= length
 * - holes_nbr*(hole_length+inter_hole_space) = length
 */
template <int dim>
bool MyGridGenerator<dim>::are_precond_fullfilled_serr_hrect(double length,
		double width, unsigned holes_nbr, double hole_length, double hole_width,
		double inter_hole_space) {

	return length > 0.0 && width > 0.0 && hole_length >= 0
			&& hole_length <= length && hole_width >= 0 && hole_width <= width
			&& inter_hole_space >= 0.0 && inter_hole_space <= length
			&& holes_nbr * (hole_length + inter_hole_space) == length;
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
 * - length > 0
 * - width > 0
 * - 0 <= hole_length <= length
 * - 0 <= hole_width <= width
 * - 0 <= inter_hole_space <= length
 * - holes_nbr*(hole_length+inter_hole_space) = length
 *
 * if pre not respected (except for the tria one), a standard rectangle of length
 * "length" and width "width" is built.
 *
 */
template<int dim>
void MyGridGenerator<dim>::serrated_hyper_rectangle(Triangulation<dim> &tria,
		double length, double width, unsigned holes_nbr, double hole_length,
		double hole_width, double inter_hole_space) {

	//void serrated_hyper_rectangle(Triangulation<dim> &tria, double length_fe,
	//		double width_fe, unsigned nbr_of_strips, unsigned strip_length,
	//		unsigned strip_width, unsigned pitch);

	if(!are_precond_fullfilled_serr_hrect(length, width, holes_nbr,
			hole_length, hole_width, inter_hole_space))
		return hyper_rectangle(tria, abs(length), abs(width));

	unsigned nbr_of_cells_fst_dim = holes_nbr * 2 + 1;
	unsigned nbr_of_cells_scd_dim = ceil(width / hole_width);

	std::vector<Point<2> > points;
	gen_points_for_serrated_hrect(length, hole_length, hole_width,
			inter_hole_space, nbr_of_cells_fst_dim, nbr_of_cells_scd_dim,
			points);

	std::vector<CellData<2> > cells;
	std::cout << nbr_of_cells_fst_dim * nbr_of_cells_scd_dim - holes_nbr << std::endl;
	cells.resize(nbr_of_cells_fst_dim * nbr_of_cells_scd_dim - holes_nbr);
	std::cout << "before cells" << std::endl;
	gen_cells_for_serrated_hrect(nbr_of_cells_fst_dim, nbr_of_cells_scd_dim,
			cells);

	std::cout << "after cells" << std::endl;


	tria.create_triangulation(points, cells, SubCellData());
}

/*
 * GridGenerator.cpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

//template <int dim>
//void GridGenerator::serrated_hyper_rectangle(Triangulation<dim> &tria, double length,
//		double width, unsigned hole_nbr, unsigned hole_length,
//		unsigned hole_width, unsigned inter_hole_space) {
//
//	//void serrated_hyper_rectangle(Triangulation<dim> &tria, double length_fe,
//		//		double width_fe, unsigned nbr_of_strips, unsigned strip_length,
//		//		unsigned strip_width, unsigned pitch);
//
//	/*
//	 * Pre:
//	 * total_length > 0
//	 * nbr_of_strips > 0
//	 * total_strip_length > 0
//	 * strip_length > 0
//	 *
//	 */
//
//	//Length "domain language" (microm...)
//	//double total_strip_length = nbr_of_strips * strip_length;
//	//double total_pitch_length = nbr_of_strips * pitch;
//	//double total_length = total_strip_length + total_pitch_length;
//
//	//length_fe are length in the finite element domain
//	//double total_strip_length_fe = total_strip_length / total_length
//	//		* length_fe;
//	//double strip_length_fe = total_strip_length_fe / nbr_of_strips;
//	//double strip_width_fe = strip_width/(double)strip_length*strip_length_fe;
//
//	//double total_pitch_length_fe = length_fe - total_strip_length_fe;
//	//double pitch_fe = total_pitch_length_fe / nbr_of_strips;
//
//	//One cell for each strip, one for each of the two half pitch at left and
//	//right and one cell for each pitch.
//	unsigned nbr_of_cells_fst_dim = nbr_of_strips * 2 + 1;
//	//unsigned nbr_of_cells_scd_dim = ceil(width_fe/strip_width_fe);
//	unsigned nbr_of_cells_scd_dim = ceil(width/hole_width);
//
//	/*std::cout << "rect_width_fe " << width_fe << std::endl;
//	std::cout << "strip_width_fe " << strip_width_fe << std::endl;
//	std::cout << "nbr_of_cells_fst_dim " << nbr_of_cells_fst_dim << std::endl;
//	std::cout << "nbr_of_cells_scd_dim " << nbr_of_cells_scd_dim << std::endl;*/
//
//	std::vector<Point<2> > points;
//
//	//End of half-pitch on left detector side
//	double left_half_pitch_end_x = pitch_fe / 2.0;
//	double right_half_pitch_start_x = length_fe-left_half_pitch_end_x;
//
//	//for (double y = 0.0; y <= width_fe; y += strip_width_fe) {
//	double y=0.0;
//	for (unsigned j = 0; j <= nbr_of_cells_scd_dim; j++) {
//
//		//First point for half-pitch on left detector side
//		Point<2> left_corner_point(0.0, y);
//		points.push_back(left_corner_point);
//
//		bool strip_left_point = true;
//
//		double x = left_half_pitch_end_x;
//		//for (double x = left_half_pitch_end_x; x <= right_half_pitch_start_x;) {
//		for(unsigned i=0; i < nbr_of_cells_fst_dim-1; i++){
//
//			Point<2> new_point(x, y);
//			points.push_back(new_point);
//
//			if (strip_left_point)
//				x += strip_length_fe;
//			else
//				x += pitch_fe;
//
//			strip_left_point = !strip_left_point;
//		}
//
//		Point<2> end_point(length_fe, y);
//		points.push_back(end_point);
//		y += strip_width_fe;
//	}
//
//	std::vector<CellData<2> > cells;
//
//	std::cout << "nbr of Points: " << points.size() << std::endl;
//
//	cells.resize(nbr_of_cells_fst_dim * nbr_of_cells_scd_dim - nbr_of_strips);
//	unsigned int c = 0;
//	for (unsigned int y = 0; y < nbr_of_cells_scd_dim; ++y){
//		for (unsigned int x = 0; x < nbr_of_cells_fst_dim; ++x) {
//
//			if ((x % 2 == 1) && (y == nbr_of_cells_scd_dim-1))
//				continue;
//			Assert(c < cells.size(), ExcInternalError());
//			cells[c].vertices[0] = y * (nbr_of_cells_fst_dim + 1) + x;
//			cells[c].vertices[1] = y * (nbr_of_cells_fst_dim + 1) + x + 1;
//			cells[c].vertices[2] = (y + 1) * (nbr_of_cells_fst_dim + 1) + x;
//			cells[c].vertices[3] = (y + 1) * (nbr_of_cells_fst_dim + 1) + x + 1;
//			cells[c].material_id = 0;
//
//			//for(int i=0; i<4; i++)
//			//	std::cout << (unsigned)cells[c].vertices[i] << " ";
//			//std::cout << std::endl;
//
//			++c;
//		}
//	}
//
//
//
//	tria.create_triangulation(points, cells, SubCellData());
//}
