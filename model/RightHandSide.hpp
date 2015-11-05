/*
 * RightHandSide.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#ifndef __RIGHT__HAND_SIDE__HPP__
#define __RIGHT__HAND_SIDE__HPP__

#include <deal.II/grid/tria.h>
#include <deal.II/base/function.h>

using namespace dealii;

template <int dim>
class RightHandSide : public Function<dim>
{

public:
  RightHandSide () : Function<dim>() {}
  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

#include "RightHandSide.cpp"

#endif
