/*
 * SerratedRect2DBoundaryValues.cpp
 *
 *  Created on: 22 nov. 2015
 *      Author: aschils
 */

template<int dim>
SerratedRect2DBoundaryValues<dim>::SerratedRect2DBoundaryValues(unsigned nbr_of_strip,
		unsigned strip_length, unsigned strip_width, unsigned pitch, double rect_length_fe,
		double rect_width_fe, double strip_potential){

	if(rect_length_fe <= 0)
		throw ZERO_OR_NEGATIVE_DOMAIN_LENGTH_ERROR;

	this->rect_length_fe = rect_length_fe;
	this->strip_potential = strip_potential;

	if(nbr_of_strip == 0 || strip_length == 0){
		strip_length_fe = 0;
		pitch_length_fe = 0;
		strip_pitch_pair_length_fe = 0;
		return;
	}

	double total_strips_length = nbr_of_strip*strip_length;
	double total_pitches_length = nbr_of_strip*pitch;
	double total_length = total_strips_length + total_pitches_length;

	double total_strips_length_fe =
			total_strips_length/total_length*rect_length_fe;
	double total_pitches_length_fe = rect_length_fe - total_strips_length_fe;

	strip_length_fe = total_strips_length_fe/nbr_of_strip;
	this->strip_width_fe = strip_width/strip_length*strip_length_fe;
	pitch_length_fe = total_pitches_length_fe/nbr_of_strip;

	//Computing these values here for performance reason
	strip_pitch_pair_length_fe = strip_length_fe + pitch_length_fe;
	half_pitch_length_fe = pitch_length_fe/2.0;
	this->rect_width_fe = rect_width_fe;
}

/**
 *
 * We define a periodic structure as a strip surrounded by half a pitch on
 * each side. It is the structure repeated on the detector upper surface to
 * implement a grid of strips/pixels. Examples:
 *
 * Detector with one periodic structure:
 * |  ---  |
 *
 * Detector with two periodic structures:
 * |  ---    ---  |
 *
 * @pre:
 * - if no strip or strip length is 0, half_pitch_length_fe must have been set
 * to 0 in the constructor.
 *
 * @return: true if the there is a strip at position p, false otherwise.
 */
template<int dim>
bool SerratedRect2DBoundaryValues<dim>::is_strip(const Point<dim> &p) const {

	double x = p[0];
	double y = p[1];

	return !(y == rect_width_fe || x == 0 || x == rect_length_fe || y == 0);


/*
	if(p[1] < rect_width_fe-strip_width_fe)
		return false;

	double component = p[0];
	//points start at component 0 instead of -rect_length_fe/2
	//double translated_component = rect_length_fe/2 + component;
	unsigned nbr_of_prev_periodic_str = (strip_pitch_pair_length_fe > 0)?
			component/strip_pitch_pair_length_fe : 0;
	double delta_from_prev_periodic_str = component -
			nbr_of_prev_periodic_str*strip_pitch_pair_length_fe;

	if(delta_from_prev_periodic_str < half_pitch_length_fe)
		return false;
	else if(delta_from_prev_periodic_str <= strip_length_fe+half_pitch_length_fe)
		return true;
	else
		return false;*/
}

template<int dim>
double SerratedRect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {
	if(is_strip(p))
		return strip_potential;
	else
		return 0.0;
}
