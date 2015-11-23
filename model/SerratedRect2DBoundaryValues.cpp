/*
 * SerratedRect2DBoundaryValues.cpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

template<int dim>
SerratedRect2DBoundaryValues<dim>::SerratedRect2DBoundaryValues(unsigned nbr_of_strips,
		double rect_length_fe,	double rect_width_fe, double strip_potential){
	this->rect_width_fe = rect_width_fe;
	this->rect_length_fe = rect_length_fe;
	this->strip_potential = strip_potential;
	this->nbr_of_strips = nbr_of_strips;
}

template<int dim>
bool SerratedRect2DBoundaryValues<dim>::is_strip(const Point<dim> &p) const {

	double x = p[0];
	double y = p[1];
	double epsilon = 0.00001;

	return !(Utils::is_same_double(y, rect_width_fe, epsilon) || x == 0 ||
			x == rect_length_fe || y == 0);
}

template<int dim>
double SerratedRect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {
	if(nbr_of_strips > 0 && is_strip(p))
		return strip_potential;
	else
		return 0.0;
}
