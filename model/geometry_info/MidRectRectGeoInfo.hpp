/*
 * MidRectRectGeoInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "MyGeometryInfo.hpp"
#include "../Utils.hpp"

class MidRectRectGeoInfo: public MyGeometryInfo {

public:

	MidRectRectGeoInfo(unsigned half_width, unsigned half_strip_width,
			unsigned strip_length, unsigned nbr_of_strips,
			unsigned half_inter_strip_dist) {
		this->half_width = half_width;
		this->length = nbr_of_strips
				* (strip_length + half_inter_strip_dist * 2);
		this->width = 2 * half_width;
		this->strip_length = strip_length;
		this->half_strip_width = half_strip_width;
		this->nbr_of_strips = nbr_of_strips;
		this->half_inter_strip_dist = half_inter_strip_dist;

		//More efficient to compute this once instead of
		//at each "is_middle_strip" call
		nbr_of_strips_before_mid_strip = nbr_of_strips / 2;
		middle_strip_left_boundary_x = nbr_of_strips_before_mid_strip
				* (strip_length + 2 * half_inter_strip_dist)
				+ half_inter_strip_dist;
		middle_strip_right_boundary_x = middle_strip_left_boundary_x
				+ strip_length;
		strip_bottom_boundary_y = half_width - half_strip_width;
		strip_top_boundary_y = half_width + half_strip_width;
	}

	bool is_strip(bpoint p) {

		double y = p.get<1>();

		if (y > strip_top_boundary_y || y < strip_bottom_boundary_y)
			return false;

		double x = p.get<0>();
		unsigned periodic_struct_lgth = strip_length + 2*half_inter_strip_dist;
		unsigned nbr_of_periodic_struct_before = x / periodic_struct_lgth;

		//Shift the coordinate system origin at the middle of the periodic
		//structure
		x = x - nbr_of_periodic_struct_before * periodic_struct_lgth;
		return x >= half_inter_strip_dist
				&& x <= periodic_struct_lgth - half_inter_strip_dist;
	}

	bool is_middle_strip(const bpoint &p) const {
		double x = p.get<0>();
		double y = p.get<1>();
		return x >= middle_strip_left_boundary_x
				&& x <= middle_strip_right_boundary_x
				&& y <= strip_top_boundary_y && y >= strip_bottom_boundary_y;
	}

	bool is_middle_strip(const Point<2> &p){
		bpoint bp = Utils::dealii_point_to_bpoint<2>(p);
		return is_middle_strip(bp);
	}

	bool is_point_inside_geometry(bpoint p) {
		return !coord_outside_geo_coord_range(p) && !is_strip(p);
	}

	//Return the maximum width of the domain
	unsigned get_width() {
		return width;
	}

	unsigned get_half_width(){
		return half_width;
	}

	//Return the maximum length of the domain
	unsigned get_length() {
		return length;
	}

	unsigned get_nbr_of_strips() {
		return nbr_of_strips;
	}

	unsigned get_dimension() {
		return 2;
	}

	unsigned get_strip_length() {
		return strip_length;
	}

	unsigned get_strip_width() {
		return 2 * half_strip_width;
	}

	unsigned get_half_strip_width() {
		return half_strip_width;
	}

	unsigned get_half_inter_strip_dist() {
		return half_inter_strip_dist;
	}

	Line get_mid_length_vertical_line() {
		double mid_length = length / 2.0;
		bpoint p1(mid_length, 0);
		bpoint p2(mid_length, 1);
		Line l(p1, p2);
		return l;
	}

	std::vector<bg::model::segment<bpoint> > get_geometry_segments() {
		std::vector<bg::model::segment<bpoint> > segments;

		//Add external rectangle borders
		bpoint top_left(0.0, width);
		bpoint top_right(length, width);
		bpoint bot_left(0.0, 0.0);
		bpoint bot_right(length, 0.0);

		bg::model::segment<bpoint> bot(bot_left, bot_right);
		bg::model::segment<bpoint> top(top_left, top_right);
		bg::model::segment<bpoint> left(bot_left, top_left);
		bg::model::segment<bpoint> right(bot_right, top_right);

		segments.push_back(bot);
		segments.push_back(top);
		segments.push_back(left);
		segments.push_back(right);

		double left_x_strip = half_inter_strip_dist;
		double right_x_strip = half_inter_strip_dist + strip_length;
		double top_y_strip = half_width + half_strip_width;
		double bot_y_strip = half_width - half_strip_width;

		//Add strip borders
		double x_shift = 0;
		double periodic_struct_length = strip_length
				+ 2 * half_inter_strip_dist;

		for (unsigned i = 0; i < nbr_of_strips; i++) {

			x_shift += i * periodic_struct_length;

			bpoint top_left_strip(left_x_strip + x_shift, top_y_strip);
			bpoint bot_left_strip(left_x_strip + x_shift, bot_y_strip);
			bpoint top_right_strip(right_x_strip + x_shift, top_y_strip);
			bpoint bot_right_strip(right_x_strip + x_shift, bot_y_strip);

			bg::model::segment<bpoint> bot_strip(bot_left_strip, bot_right_strip);
			bg::model::segment<bpoint> top_strip(top_left_strip, top_right_strip);
			bg::model::segment<bpoint> left_strip(bot_left_strip, top_left_strip);
			bg::model::segment<bpoint> right_strip(bot_right_strip, top_right_strip);

			segments.push_back(bot_strip);
			segments.push_back(top_strip);
			segments.push_back(left_strip);
			segments.push_back(right_strip);
		}

		return segments;
	}

private:
	unsigned half_width, length, half_strip_width, strip_length, nbr_of_strips;
	unsigned half_inter_strip_dist, width;

	unsigned nbr_of_strips_before_mid_strip, middle_strip_left_boundary_x;
	unsigned middle_strip_right_boundary_x, strip_bottom_boundary_y;
	unsigned strip_top_boundary_y;

	bool coord_outside_geo_coord_range(bpoint p) {
		double x = p.get<0>();
		double y = p.get<1>();
		return x < 0 || x > length || y < 0 || y > width;
	}

};
