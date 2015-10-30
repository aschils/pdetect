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


  Triangulation<2>::active_cell_iterator cell, endc;

  unsigned level = 0;
  unsigned nbr_of_cell_to_refine = 50;
  cell = triangulation.begin_active(level),
  endc = triangulation.end();

  for (unsigned i = 0; cell!=endc && i!=nbr_of_cell_to_refine; cell++, i++)
	  cell->set_refine_flag ();

  triangulation.execute_coarsening_and_refinement ();


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
