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
using namespace dealii;

void first_grid (){
  Triangulation<2> triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (4);
  std::ofstream out ("grid-1.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);
  std::cout << "Grid written to grid-1.eps" << std::endl;
}


int main() {
	first_grid ();
	return 0;
}
