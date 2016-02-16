/*
 * SerratedRect2DGeoInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "MyGeometryInfo.hpp"
#include "../errors.hpp"
#include "../Utils.hpp"

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

	bool is_point_inside_geometry(Point<2> p) {

		switch (dim) {
		case 2:
			return is_point_inside_detector_2D(p);
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

		if (coord_outside_geo_coord_range(p))
			return false;

		unsigned pitch = 2 * half_pitch;
		unsigned periodic_str_length = pitch + strip_length;

		double epsilon = 0.00000001;

		double x = p[0];
		double y = p[1];

		double bottom_of_strip = width - strip_width;

		if (periodic_str_length == 0.0 || strip_length == 0.0
				|| !Utils::greater_than_or_equals_double(y, bottom_of_strip,
						epsilon)
				|| !Utils::less_than_or_equals_double(y, width, epsilon))
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
		unsigned periodic_str_before_mid_strip = this->nbr_of_strips / 2;
		unsigned nbr_of_prev_periodic_str = p[0] / periodic_str_length;

		//Fixed bug when half_pitch == 0, we must take the last point at the right
		//side of the strip
		if (half_pitch == 0
				&& p[0]
						== (periodic_str_before_mid_strip + 1)
								* periodic_str_length) {
			nbr_of_prev_periodic_str -= 1;
		}

		return this->nbr_of_strips > 0 && is_strip<dim>(p)
				&& nbr_of_prev_periodic_str == periodic_str_before_mid_strip;
	}

	/**
	 * intersections are placed in increasing x position
	 */
	std::vector<Point<2>> boundaries_intersections(Line line) {

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

		//half pitch at right
		Point<2> right_half_pitch_left(length - half_pitch, width);
		Point<2> right_half_pitch_right(length, width);
		Segment right_half_pitch(right_half_pitch_left, right_half_pitch_right);
		segments.push_back(right_half_pitch);

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

		for (unsigned i = 0; i < segments.size(); i++) {
			add_intersection(line, segments[i], intersections);
			//std::cout << "[(" << segments[i].p1[0] << "," <<
			//		segments[i].p1[1] << "),(" <<  segments[i].p2[0] << "," <<
			//		segments[i].p2[1] << ")]" << std::endl;
		}

		Utils::sort_points_by_coord<2>(&intersections);
		intersections.erase(unique(intersections.begin(), intersections.end()),
				intersections.end());
		intersections.shrink_to_fit();
		return intersections;
	}

	bool is_point_inside_detector_and_not_strip_2D(Point<2> p) {
		return is_point_inside_detector_2D(p) && !is_strip<2>(p);
	}

	Line get_mid_length_vertical_line() {
		double mid_length = length/2.0;
		Point<2> p1(mid_length, 0);
		Point<2> p2(mid_length, 1);
		Segment seg(p1,p2);
		return Line(seg);
	}

private:
	unsigned nbr_of_strips, width, strip_length, strip_width, half_pitch,
			length;
	unsigned dim;

	bool coord_outside_geo_coord_range(Point<2> p) {
		double x = p[0];
		double y = p[1];
		return x < 0 || x > length || y < 0 || y > width;
	}

	bool is_point_inside_detector_2D(Point<2> point) {

		//Quick check to exclude points obviously not inside detector domain
		if (coord_outside_geo_coord_range(point))
			return false;

		double y = point[1];

		if (y > width - strip_width) {
			return is_strip<2>(point);
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
		return line.intersection_point(seg, intersect_pt);
	}

	void add_intersection(Line line, Segment seg,
			std::vector<Point<2>> &intersections) {
		Point<2> intersect;
		if (intersect_segment(line, seg, intersect)) {

			//We must remove the intersection point
			//if at + and -epsilon not in detector OR
			// if at + and - epsilon in detector

			double epsilon = 0.001;
			Point<2> eps_direction_vec = epsilon * line.get_direction_vector();
			Point<2> before(intersect[0] - eps_direction_vec[0],
					intersect[1] - eps_direction_vec[1]);
			Point<2> after = intersect + eps_direction_vec;

			bool before_in_det = is_point_inside_detector_and_not_strip_2D(
					before);
			bool after_in_det = is_point_inside_detector_and_not_strip_2D(
					after);
			bool is_a_cross_boundary_point = !((before_in_det && after_in_det)
					|| (!before_in_det && !after_in_det));

			if (is_a_cross_boundary_point) {
				intersections.push_back(intersect);
			}
		}
	}

};
