/*
 * RightHandSide.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
  double return_value = 0.0;
  for (unsigned int i=0; i<dim; ++i)
    return_value += 4.0 * std::pow(p(i), 4.0);
  return return_value;
}
