/*
 * RightHandSide.hpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

#pragma once

#include <deal.II/grid/tria.h>
#include <deal.II/base/function.h>

using namespace dealii;

template <unsigned dim>
class ZeroRightHandSide : public Function<dim>
{

public:
  ZeroRightHandSide () : Function<dim>() {}
  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <unsigned dim>
double ZeroRightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
	return 0.0;
}
