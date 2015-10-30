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
<<<<<<< HEAD
	first_grid ();
=======
	cout << "Third" << endl;
>>>>>>> f6c5d69d56851a39aff393fd57ac6698e3fc48fc
	return 0;
}
