/*
 * CustomGridGen.cpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

#include "CustomGridGen.hpp"

void serrated_rectangle(Triangulation<2> &tria, double length_fe, double width_fe,
		unsigned nbr_of_strips, unsigned strip_length, unsigned strip_width,
		unsigned pitch) {

	/*
	 * Pre:
	 * total_length > 0
	 * nbr_of_strips > 0
	 * total_strip_length > 0
	 * strip_length > 0
	 *
	 */

	//Length "domain language" (microm...)
	double total_strip_length = nbr_of_strips * strip_length;
	double total_pitch_length = nbr_of_strips * pitch;
	double total_length = total_strip_length + total_pitch_length;

	//length_fe are length in the finite element domain
	double total_strip_length_fe = total_strip_length / total_length
			* length_fe;
	double strip_length_fe = total_strip_length_fe / nbr_of_strips;
	double strip_width_fe = strip_width/(double)strip_length*strip_length_fe;

	double total_pitch_length_fe = length_fe - total_strip_length_fe;
	double pitch_fe = total_pitch_length_fe / nbr_of_strips;

	//One cell for each strip, one for each of the two half pitch at left and
	//right and one cell for each pitch.
	unsigned nbr_of_cells_fst_dim = nbr_of_strips * 2 + 1;
	unsigned nbr_of_cells_scd_dim = ceil(width_fe/strip_width_fe);

	/*std::cout << "rect_width_fe " << width_fe << std::endl;
	std::cout << "strip_width_fe " << strip_width_fe << std::endl;
	std::cout << "nbr_of_cells_fst_dim " << nbr_of_cells_fst_dim << std::endl;
	std::cout << "nbr_of_cells_scd_dim " << nbr_of_cells_scd_dim << std::endl;*/

	std::vector<Point<2> > points;

	//End of half-pitch on left detector side
	double left_half_pitch_end_x = pitch_fe / 2.0;
	double right_half_pitch_start_x = length_fe-left_half_pitch_end_x;

	//for (double y = 0.0; y <= width_fe; y += strip_width_fe) {
	double y=0.0;
	for (unsigned j = 0; j <= nbr_of_cells_scd_dim; j++) {

		//First point for half-pitch on left detector side
		Point<2> left_corner_point(0.0, y);
		points.push_back(left_corner_point);

		bool strip_left_point = true;

		double x = left_half_pitch_end_x;
		//for (double x = left_half_pitch_end_x; x <= right_half_pitch_start_x;) {
		for(unsigned i=0; i < nbr_of_cells_fst_dim-1; i++){

			Point<2> new_point(x, y);
			points.push_back(new_point);

			if (strip_left_point)
				x += strip_length_fe;
			else
				x += pitch_fe;

			strip_left_point = !strip_left_point;
		}

		Point<2> end_point(length_fe, y);
		points.push_back(end_point);
		y += strip_width_fe;
	}

	std::vector<CellData<2> > cells;

	std::cout << "nbr of Points: " << points.size() << std::endl;

	cells.resize(nbr_of_cells_fst_dim * nbr_of_cells_scd_dim - nbr_of_strips);
	unsigned int c = 0;
	for (unsigned int y = 0; y < nbr_of_cells_scd_dim; ++y){
		for (unsigned int x = 0; x < nbr_of_cells_fst_dim; ++x) {

			if ((x % 2 == 1) && (y == nbr_of_cells_scd_dim-1))
				continue;
			Assert(c < cells.size(), ExcInternalError());
			cells[c].vertices[0] = y * (nbr_of_cells_fst_dim + 1) + x;
			cells[c].vertices[1] = y * (nbr_of_cells_fst_dim + 1) + x + 1;
			cells[c].vertices[2] = (y + 1) * (nbr_of_cells_fst_dim + 1) + x;
			cells[c].vertices[3] = (y + 1) * (nbr_of_cells_fst_dim + 1) + x + 1;
			cells[c].material_id = 0;

			//for(int i=0; i<4; i++)
			//	std::cout << (unsigned)cells[c].vertices[i] << " ";
			//std::cout << std::endl;

			++c;
		}
	}



	tria.create_triangulation(points, cells, SubCellData());
}

/*
 * CustomGridGen.cpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

//#include "CustomGridGen.hpp"
//
//void serrated_rectangle(Triangulation<2> &tria, unsigned pitch_length,
//		unsigned nbr_of_strips, unsigned stip_length, unsigned strip_width) {
//
//	const unsigned dim = 2;
//
//	AssertDimension(holes.size(), dim);
//
//	// The corner points of the first cell. If there is a desire at
//	// some point to change the geometry of the cells, they can be
//	// made an argument to the function.
//
//	Point<dim> p1;
//	Point<dim> p2;
//
//	for (unsigned int d = 0; d < dim; ++d)
//		p2(d) = 1.;
//
//	// then check that all repetitions
//	// are >= 1, and calculate deltas
//	// convert repetitions from double
//	// to int by taking the ceiling.
//	std::vector<Point<dim> > delta(dim);
//	unsigned int repetitions[dim];
//	for (unsigned int i = 0; i < dim; ++i) {
//		//Assert(holes[i] >= 1,
//		//		ExcMessage("At least one hole needed in each direction"));
//		repetitions[i] = 2 * holes[i] + 1;
//		delta[i][i] = (p2[i] - p1[i]);
//	}
//
//	// then generate the necessary
//	// points
//	std::vector<Point<dim> > points;
//	switch (dim) {
//	case 2:
//		for (unsigned int y = 0; y <= repetitions[1]; ++y)
//			for (unsigned int x = 0; x <= repetitions[0]; ++x)
//				points.push_back(
//						p1 + (double) x * delta[0] + (double) y * delta[1]);
//		break;
//
//	default:
//		Assert(false, ExcNotImplemented());
//	}
//
//	// next create the cells
//	// Prepare cell data
//	std::vector<CellData<dim> > cells;
//	switch (dim) {
//	case 2: {
//		cells.resize(repetitions[1] * repetitions[0]-2);// - holes[1] * holes[0]);
//		unsigned int c = 0;
//		for (unsigned int y = 0; y < repetitions[1]; ++y)
//			for (unsigned int x = 0; x < repetitions[0]; ++x) {
//				//if ((x % 2 == 1) && (y % 2 == 1))
//				if ((x % 2 == 1) && (y == 4))
//					continue;
//				Assert(c < cells.size(), ExcInternalError());
//				cells[c].vertices[0] = y * (repetitions[0] + 1) + x;
//				cells[c].vertices[1] = y * (repetitions[0] + 1) + x + 1;
//				cells[c].vertices[2] = (y + 1) * (repetitions[0] + 1) + x;
//				cells[c].vertices[3] = (y + 1) * (repetitions[0] + 1) + x + 1;
//				cells[c].material_id = 0;
//				++c;
//			}
//		break;
//	}
//
//	default:
//		Assert(false, ExcNotImplemented());
//	}
//
//	tria.create_triangulation(points, cells, SubCellData());
//}
//
