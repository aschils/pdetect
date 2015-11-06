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


	double stripWidth = 2.0/3.0; //Strip size is 1/3 of domain size

	//stripCenter.unrolled_to_component_indices(1)

	if(p(1) == -1 && -1+stripWidth <= p(2) && p(2) <= 1-stripWidth)
		return 1; //1V
	else
		return 0;

  //return p.square();
}
