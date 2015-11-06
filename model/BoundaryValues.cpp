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
	double stripWidth = 1.95; //Strip size is 1/3 of domain size
	return p[1] == -1 && -stripWidth/2 <= p[0] && p[0] <= stripWidth/2;
}
