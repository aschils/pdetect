//============================================================================
// Name        : pdetect.cpp
// Author      : Arnaud Schils
//============================================================================

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace dealii;

void first_grid(){
  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  unsigned level = 0;

  Triangulation<2>::active_cell_iterator cell, endc;
  cell = triangulation.begin_active(level),
  endc = triangulation.end();

  const Point<2> bottom_left(0,0);

  for (; cell!=endc; cell++){
	  for(unsigned j=0; j<GeometryInfo<2>::vertices_per_cell; j++){

          const double distance_from_bottom_left
            = (cell->vertex(j)).distance(bottom_left);

		  if(distance_from_bottom_left < triangulation.space_dimension/4.0){
			  cell->set_refine_flag ();
			  break;
		  }
	  }
  }


  triangulation.execute_coarsening_and_refinement();


  std::string fileName = "grid-";
  fileName = std::to_string(level)+".eps";

  std::ofstream out(fileName);
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);
  std::cout << "Grid written to " << fileName << std::endl;
}


int main() {
	first_grid();
	return 0;
}
