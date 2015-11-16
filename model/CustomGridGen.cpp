/*
 * CustomGridGen.cpp
 *
 *  Created on: 13 nov. 2015
 *      Author: aschils
 */

#include "CustomGridGen.hpp"

void serrated_rectangle(Triangulation<2> &tria,
		unsigned nbr_of_strips) {

	std::vector<Point<2> > points;

	unsigned nbr_of_cells_fst_dim = nbr_of_strips*2+1;
	unsigned nbr_of_cells_scd_dim = 5;

	for (unsigned y = 0; y <= nbr_of_cells_scd_dim; ++y){
			for (unsigned x = 0; x <= nbr_of_cells_fst_dim; ++x){
				Point<2> new_point(x,y);
				points.push_back(new_point);
			}
	}

	std::vector<CellData<2> > cells;

	cells.resize(nbr_of_cells_fst_dim * nbr_of_cells_scd_dim - nbr_of_strips);
	unsigned int c = 0;
	for (unsigned int y = 0; y < nbr_of_cells_scd_dim; ++y)
		for (unsigned int x = 0; x < nbr_of_cells_fst_dim; ++x) {
			//if ((x % 2 == 1) && (y % 2 == 1))
			if ((x % 2 == 1) && (y == 4))
				continue;
			Assert(c < cells.size(), ExcInternalError());
			cells[c].vertices[0] = y * (nbr_of_cells_fst_dim + 1) + x;
			cells[c].vertices[1] = y * (nbr_of_cells_fst_dim + 1) + x + 1;
			cells[c].vertices[2] = (y + 1) * (nbr_of_cells_fst_dim + 1) + x;
			cells[c].vertices[3] = (y + 1) * (nbr_of_cells_fst_dim + 1) + x + 1;
			cells[c].material_id = 0;
			++c;
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
