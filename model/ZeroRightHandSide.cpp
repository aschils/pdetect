/*
 * RightHandSide.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template <unsigned dim>
double ZeroRightHandSide<dim>::value (const Point<dim> &p,
                                  const unsigned int /*component*/) const
{
	return 0.0;
}
