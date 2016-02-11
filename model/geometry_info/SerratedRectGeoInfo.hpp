/*
 * SerratedRect2DGeoInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "MyGeometryInfo.hpp"
#include "../errors.hpp"

class SerratedRectGeoInfo: public MyGeometryInfo {

public:

	SerratedRectGeoInfo(unsigned dim, unsigned nbr_of_strips, unsigned width,
			unsigned strip_length, unsigned strip_width, unsigned half_pitch) {
		this->nbr_of_strips = nbr_of_strips;
		this->width = width;
		this->strip_length = strip_length;
		this->strip_width = strip_width;
		this->half_pitch = half_pitch;
		this->length = compute_total_length();
		this->dim = dim;
	}

	unsigned get_width() {
		return width;
	}

	unsigned get_length() {
		return length;
	}

	unsigned get_nbr_of_strips() {
		return nbr_of_strips;
	}

	unsigned get_dimension(){
		return dim;
	}

	unsigned get_strip_length(){
		return strip_length;
	}

	unsigned get_strip_width(){
		return strip_width;
	}

	unsigned get_half_pitch(){
		return half_pitch;
	}

	bool is_point_inside_geometry(std::vector<double> point_coord) {

		switch (dim) {
		case 2:
			return is_point_inside_detector_2D(point_coord);
		default:
			throw NOT_YET_IMPLEMENTED;
		}

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
	 */
	template<unsigned dim>
	bool is_strip(const Point<dim> &p) {

		unsigned pitch = 2 * half_pitch;
		unsigned periodic_str_length = pitch + strip_length;

		double epsilon = 0.00000001;

		double x = p[0];
		double y = p[1];

		double bottom_of_strip = width - strip_width;

		if (!Utils::greater_than_or_equals_double(y, bottom_of_strip, epsilon)
				|| periodic_str_length == 0.0 || strip_length == 0.0)
			return false;

		unsigned nbr_of_prev_periodic_str = x / periodic_str_length;
		double delta_from_prev_periodic_str = x
				- nbr_of_prev_periodic_str * periodic_str_length;

		double strip_border_right = half_pitch + strip_length;

		return Utils::greater_than_or_equals_double(
				delta_from_prev_periodic_str, half_pitch, epsilon)
				&& Utils::less_than_or_equals_double(
						delta_from_prev_periodic_str, strip_border_right,
						epsilon);
	}

	template<unsigned dim>
	bool is_middle_strip(const Point<dim> &p) {

		unsigned periodic_str_length = 2 * half_pitch + strip_length;
		unsigned nbr_of_prev_periodic_str = p[0] / periodic_str_length;

		return this->nbr_of_strips > 0 && is_strip<dim>(p)
				&& nbr_of_prev_periodic_str == this->nbr_of_strips / 2;
	}

	std::vector<Utils::Segment<2>> segments_at_intersec(
				Utils::Line<2> line){
		//TODO
	}

private:
	unsigned nbr_of_strips, width, strip_length, strip_width, half_pitch,
			length;
	unsigned dim;

	bool is_point_inside_detector_2D(std::vector<double> point_coord) {

		double x = point_coord[0];
		double y = point_coord[1];

		//Quick check to exclude points obviously not inside detector domain
		if (x < 0 || x > length || y < 0 || y > width)
			return false;

		if (y > width - strip_width) {
			dealii::Point<2> point(x, y);
			return !is_strip<2>(point);
		} else {
			return true;
		}
	}

	unsigned compute_total_length() {
		unsigned pitch = 2 * half_pitch;
		if (nbr_of_strips == 0)
			return pitch;
		else
			return nbr_of_strips * (strip_length + pitch);
	}

};
