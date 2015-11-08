/*
 * BoundaryValues.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template <int dim>
double Rect2DBoundaryValues<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
	//2D case
	//double stripWidth = 0.5; //Strip size is 1/3 of domain size
	//return p[1] == -1 && -stripWidth/2 <= p[0] && p[0] <= stripWidth/2;
	double stripWidth = 2.0/5.0; //Strip size is 1/3 of domain size

	double middle  = p[0] >= -2+9*stripWidth/2 && p[0] <= 2-9*stripWidth/2;
	double m_left  = p[0] >= -2+5*stripWidth/2 && p[0] <= -2+7*stripWidth/2;
	double m_right = p[0] <= 2-5*stripWidth/2 && p[0] >= 2-7*stripWidth/2;
	double o_left  = p[0] >= -2+stripWidth/2 && p[0] <= -2+3*stripWidth/2;
	double o_right = p[0] <= 2-stripWidth/2 && p[0] >= 2-3*stripWidth/2;
	
	return p[1] == -1 && (middle || m_right || m_left || o_right || o_left);
}
