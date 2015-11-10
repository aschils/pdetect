/*
 * BoundaryValues.cpp
 *
 *  Created on: 5 nov. 2015
 *      Author: aschils
 */
template<int dim>
Rect2DBoundaryValues<dim>::Rect2DBoundaryValues(unsigned nbr_of_strip,
		unsigned strip_width, unsigned pitch, double rect_length_fe){

	if(rect_length_fe <= 0)
		throw ZERO_OR_NEGATIVE_DOMAIN_LENGTH_ERROR;

	this->rect_length_fe = rect_length_fe;

	if(nbr_of_strip == 0 || strip_width == 0){
		this->strip_length_fe = 0;
		this->pitch_length_fe = 0;
		this->strip_pitch_pair_length_fe = 0;
		return;
	}


	double total_strips_length = nbr_of_strip*strip_width;
	double total_pitches_length = (nbr_of_strip+1)*pitch;
	double total_length = total_strips_length + total_pitches_length;

	double total_strips_length_fe =
			total_strips_length/total_length*rect_length_fe;
	double total_pitches_length_fe = rect_length_fe - total_strips_length_fe;

	this->strip_length_fe = total_strips_length_fe/nbr_of_strip;
	this->pitch_length_fe = (nbr_of_strip >= 1)?
			total_pitches_length_fe/(nbr_of_strip+1) : 0;
	this->strip_pitch_pair_length_fe = this->strip_length_fe + this->pitch_length_fe;
}

template<int dim>
bool Rect2DBoundaryValues<dim>::is_strip(double pos) const {

	//points start at component 0 instead of -rect_length_fe/2
	/*double translated_component = this->rect_length_fe/2 + component;
	unsigned strip_picth_pair_before = (this->strip_pitch_pair_length_fe > 0)?
			translated_component/this->strip_pitch_pair_length_fe : 0;
	double delta_from_strip_pitch_pair =
			translated_component - strip_picth_pair_before*this->strip_pitch_pair_length_fe;
	return delta_from_strip_pitch_pair < this->strip_length_fe;*/

	if(this->strip_length_fe == 0)
		return 0;

	double i = -this->rect_length_fe/2 + this->pitch_length_fe;
	double j = -i;

	for(; i <= j; i = i+this->strip_pitch_pair_length_fe){
		double in_strip = pos >= i && pos <= i+this->strip_length_fe;
		if(in_strip)
			return 1;
	}
	return 0;
}

template<int dim>
double Rect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {

	return p[1] == 1 && is_strip(p[0]);
}
