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
		strip_length_fe = 0;
		pitch_length_fe = 0;
		strip_pitch_pair_length_fe = 0;
		return;
	}


	double total_strips_length = nbr_of_strip*strip_width;
	double total_pitches_length = (nbr_of_strip-1)*pitch;
	double total_length = total_strips_length+total_pitches_length;

	double total_strips_length_fe =
			total_strips_length/total_length*rect_length_fe;
	double total_pitches_length_fe = rect_length_fe-total_strips_length_fe;

	strip_length_fe = total_strips_length_fe/nbr_of_strip;
	pitch_length_fe = (nbr_of_strip > 1)?
			total_pitches_length_fe/(nbr_of_strip-1): 0;
	strip_pitch_pair_length_fe = strip_length_fe+pitch_length_fe;
}

template<int dim>
bool Rect2DBoundaryValues<dim>::is_strip(double component) const {

	//points start at component 0 instead of -rect_length_fe/2
	double translated_component = rect_length_fe/2+component;
	unsigned strip_picth_pair_before = (strip_pitch_pair_length_fe > 0)?
			translated_component/strip_pitch_pair_length_fe: 0;
	double delta_from_strip_pitch_pair =
			translated_component - strip_picth_pair_before*strip_pitch_pair_length_fe;
	return delta_from_strip_pitch_pair < strip_length_fe;
}

template<int dim>
double Rect2DBoundaryValues<dim>::value(const Point<dim> &p,
		const unsigned int /*component*/) const {

	return p[1] == -1 && is_strip(p[0]);

	//2D case
	//double stripWidth = 0.5; //Strip size is 1/3 of domain size
	//return p[1] == -1 && -stripWidth/2 <= p[0] && p[0] <= stripWidth/2;
	/*double stripWidth = 2.0 / 5.0; //Strip size is 1/3 of domain size

	double middle = p[0] >= -2 + 9 * stripWidth / 2
			&& p[0] <= 2 - 9 * stripWidth / 2;
	double m_left = p[0] >= -2 + 5 * stripWidth / 2
			&& p[0] <= -2 + 7 * stripWidth / 2;
	double m_right = p[0] <= 2 - 5 * stripWidth / 2
			&& p[0] >= 2 - 7 * stripWidth / 2;
	double o_left = p[0] >= -2 + stripWidth / 2
			&& p[0] <= -2 + 3 * stripWidth / 2;
	double o_right = p[0] <= 2 - stripWidth / 2
			&& p[0] >= 2 - 3 * stripWidth / 2;

	return p[1] == -1 && (middle || m_right || m_left || o_right || o_left);*/
}
