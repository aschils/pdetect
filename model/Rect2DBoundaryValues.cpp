/*
 * BoundaryValues.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */

template<int dim>
Rect2DBoundaryValues<dim>::Rect2DBoundaryValues(unsigned nbr_of_strip,
		unsigned strip_width, unsigned pitch, double rect_length_fe,
		double rect_width_fe, double strip_potential){

	if(rect_length_fe <= 0)
		throw ZERO_OR_NEGATIVE_DOMAIN_LENGTH_ERROR;

	this->rect_length_fe = rect_length_fe;
	this->strip_potential = strip_potential;

	if(nbr_of_strip == 0 || strip_width == 0){
		strip_length_fe = 0;
		pitch_length_fe = 0;
		strip_pitch_pair_length_fe = 0;
		return;
	}

	double total_strips_length = nbr_of_strip*strip_width;
	double total_pitches_length = nbr_of_strip*pitch;
	double total_length = total_strips_length + total_pitches_length;

	double total_strips_length_fe =
			total_strips_length/total_length*rect_length_fe;
	double total_pitches_length_fe = rect_length_fe - total_strips_length_fe;

	strip_length_fe = total_strips_length_fe/nbr_of_strip;
	pitch_length_fe = total_pitches_length_fe/nbr_of_strip;

	//Computing these values here for performance reason
	strip_pitch_pair_length_fe = strip_length_fe + pitch_length_fe;
	half_pitch_length_fe = pitch_length_fe/2.0;
	half_rect_width_fe = rect_width_fe/2.0;
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
bool Rect2DBoundaryValues<dim>::is_strip(const Point<dim> &p) const {

	if(p[1] != half_rect_width_fe)
		return false;

	double component = p[0];
	//points start at component 0 instead of -rect_length_fe/2
	double translated_component = rect_length_fe/2 + component;
	unsigned nbr_of_prev_periodic_str = (strip_pitch_pair_length_fe > 0)?
			translated_component/strip_pitch_pair_length_fe : 0;
	double delta_from_prev_periodic_str = translated_component -
			nbr_of_prev_periodic_str*strip_pitch_pair_length_fe;

	if(delta_from_prev_periodic_str <= half_pitch_length_fe)
		return false;
	else if(delta_from_prev_periodic_str < strip_length_fe+half_pitch_length_fe)
		return true;
	else
		return false;
}

template<int dim>
double Rect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {
	if(is_strip(p)){
		return strip_potential;
		std::cout << "potential setted"<< std::endl;
	}
	else
		return 0.0;
}
