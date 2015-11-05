/*
 * BoundaryValues.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
  return p.square();
}
