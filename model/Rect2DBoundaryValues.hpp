/*
 * BoundaryValues.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#ifndef __BOUNDARY_VALUES__HPP__
#define __BOUNDARY_VALUES__HPP__

#include <deal.II/grid/tria.h>
#include <deal.II/base/function.h>

using namespace dealii;

template <int dim>
class Rect2DBoundaryValues : public Function<dim>
{
public:
	Rect2DBoundaryValues () : Function<dim>() {}
  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

#include "Rect2DBoundaryValues.cpp"

#endif
