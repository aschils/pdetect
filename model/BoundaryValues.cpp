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
	//2D case
	//double stripWidth = 0.5; //Strip size is 1/3 of domain size
	//return p[1] == -1 && -stripWidth/2 <= p[0] && p[0] <= stripWidth/2;
	double stripWidth = 1.0/3.0; //Strip size is 1/3 of domain size

	double middle = -stripWidth/2 <= p[0] && p[0] <= stripWidth/2;
	double left = p[0] >= -1+stripWidth/2 && p[0] <= -1+3*stripWidth/2;
	double right = p[0] <= 1-stripWidth/2 && p[0] >= 1-3*stripWidth/2;
	return p[1] == -1 && (middle || right || left);
}
