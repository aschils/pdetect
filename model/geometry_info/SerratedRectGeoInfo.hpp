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

	unsigned get_dimension() {
		return dim;
	}

	unsigned get_strip_length() {
		return strip_length;
	}

	unsigned get_strip_width() {
		return strip_width;
	}

	unsigned get_half_pitch() {
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

	/**
	 * intersections are placed in increasing x position
	 */
	std::vector<Point<2>> segments_at_intersec(Line line) {

		std::vector<Segment> segments;

		//left side
		Point<2> left_side_bot(0.0, 0.0);
		Point<2> left_side_top(0.0, width);
		Segment left_side(left_side_bot, left_side_top);
		segments.push_back(left_side);

		//half pitch at left
		Point<2> left_half_pitch_left(0.0, width);
		Point<2> left_half_pitch_right(half_pitch, width);
		Segment left_half_pitch(left_half_pitch_left, left_half_pitch_right);
		segments.push_back(left_half_pitch);

		unsigned x = half_pitch;
		unsigned strip_y = width - strip_width;
		unsigned pitch = 2 * half_pitch;
		//strips and pitchs
		for (unsigned i = 0; i < nbr_of_strips; i++) {

			//Add vertical segment "between pitch and strip"
			Point<2> vert_bot(x, strip_y);
			Point<2> vert_top(x, width);
			Segment vert_seg(vert_bot, vert_top);
			segments.push_back(vert_seg);

			if (i % 2 == 0) { //strip
				Point<2> left(x, strip_y);
				Point<2> right(x + strip_length, strip_y);
				Segment strip_seg(left, right);
				segments.push_back(strip_seg);
				x += strip_length;
			} else { //pitch
				Point<2> left(x, width);
				Point<2> right(x + pitch, width);
				Segment pitch_seg(left, right);
				segments.push_back(pitch_seg);
				x += pitch;
			}
		}

		//Add last vert segment
		if (nbr_of_strips > 0) {
			Point<2> vert_bot(x, strip_y);
			Point<2> vert_top(x, width);
			Segment vert_seg(vert_bot, vert_top);
			segments.push_back(vert_seg);
		}

		//Right side
		Point<2> right_side_bot(length, 0.0);
		Point<2> right_side_top(length, width);
		Segment right_side(right_side_bot, right_side_top);
		segments.push_back(right_side);

		//Bottom side
		Point<2> bot_side_left(0.0, 0.0);
		Point<2> bot_side_right(length, 0.0);
		Segment bot_side_seg(bot_side_left, bot_side_right);
		segments.push_back(bot_side_seg);

		//Compute intersections with all segments of detector boundary
		std::vector<Point<2>> intersections;

		for(unsigned i=0; i<segments.size(); i++)
			add_intersection(line, segments[i], intersections);

		Utils::sort_points_by_coord<2>(&intersections);
		return intersections;
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

	bool intersect_segment(Line line, Segment seg, Point<2> &intersect_pt) {

		//Build Line passing through segment seg
		Line seg_line(seg);
		try {
			intersect_pt = line.intersection_point(seg_line);
		} catch (int e) {
			return false;
		}

		double seg_low_y = seg.get_low_y();
		double seg_high_y = seg.get_high_y();
		double seg_low_x = seg.get_low_x();
		double seg_high_x = seg.get_high_x();

		if (seg.is_vertical())
			return seg_low_y <= intersect_pt[1] && intersect_pt[1] <= seg_high_y;
		else
			return seg_low_x <= intersect_pt[0] && intersect_pt[0] <= seg_high_x;
	}

	void add_intersection(Line line, Segment seg,
			std::vector<Point<2>> &intersections) {
		Point<2> intersect;
		if (intersect_segment(line, seg, intersect))
			intersections.push_back(intersect);
	}

};
