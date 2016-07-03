/*
 * BasicGeometry.hpp
 *
 *  Created on: 12 févr. 2016
 *      Author: aschils
 */

#pragma once

using namespace dealii;

#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/algorithms/length.hpp>

#include "../errors.hpp"

namespace bg = boost::geometry;

/*class Segment {
 public:
 Point<2> p1, p2;

 Segment(Point<2> p1, Point<2> p2) {
 this->p1 = p1;
 this->p2 = p2;
 double x_diff = p1[0] - p2[0];
 vertical = (x_diff == 0.0);

 if (!vertical)
 slope = (p1[1] - p2[1]) / x_diff;
 }

 bool is_vertical() {
 return vertical;
 }

 //Call this only if segment is not vertical
 double get_slope() {
 if (vertical)
 throw PRECONDITIONS_VIOLATED;
 return slope;
 }

 double get_low_x() {
 if (p1[0] < p2[0])
 return p1[0];
 else
 return p2[0];
 }

 double get_low_y() {
 if (p1[1] < p2[1])
 return p1[1];
 else
 return p2[1];
 }

 double get_high_x() {
 if (p1[0] > p2[0])
 return p1[0];
 else
 return p2[0];
 }

 double get_high_y() {
 if (p1[1] > p2[1])
 return p1[1];
 else
 return p2[1];
 }

 double compute_length() {
 double diff_x = p2[0] - p1[0];
 double diff_y = p2[1] - p1[1];
 return std::sqrt(diff_x * diff_x + diff_y * diff_y);
 }

 private:
 bool vertical;
 double slope;
 };*/

class Line {

public:

	typedef bg::model::point<double, 2, bg::cs::cartesian> bpoint;

	Line() {
	}

	Line(double slope, double intercept) {
		bpoint p1(0.0, intercept);
		double y = slope + intercept;
		bpoint p2(1, y);
		bg::model::segment<bpoint> seg(p1, p2);
		constructor_from_segment(seg);
		//std::cout << "Point (" << point[0] << "," << point[1] << ")"
		//		<< "Vec (" << vec[0] << "," << vec[1] << ")" << std::endl;
	}

	Line(bg::model::segment<bpoint> seg) {
		constructor_from_segment(seg);
	}

	Line(bpoint p1, bpoint p2) {
		bg::model::segment<bpoint> seg(p1, p2);
		constructor_from_segment(seg);
	}

	void constructor_from_segment(bg::model::segment<bpoint> seg) {

		point = seg.first;
		bpoint point2 = seg.second;
		double norm = bg::length(seg);
		double vec_x = (point2.get<0>() - point.get<0>()) / norm;
		double vec_y = (point2.get<1>() - point.get<1>()) / norm;
		bpoint vec(vec_x, vec_y);
		this->vec = vec;

		if (!is_vertical()) {
			slope = vec_y / vec_x;
			intercept = point.get<1>() - point.get<0>() * vec_y / vec_x;
		}
	}

	bool is_parallel(Line l) {
		double epsilon = 0.00001;
		bool same_vec_x = Utils::equals_double(vec.get<0>(), l.vec.get<0>(),
				epsilon);
		bool same_vec_y = Utils::equals_double(vec.get<1>(), l.vec.get<1>(),
				epsilon);
		return same_vec_x && same_vec_y;
	}

	bool intersection_point(Line l, bpoint &intersect_pt) {

		if (is_parallel(l)) {
			return false;
		}

		//Well it's geometry: intersection between two lines described
		//by parametric equations
		double v1 = vec.get<0>();
		double v2 = vec.get<1>();
		double w1 = l.vec.get<0>();
		double w2 = l.vec.get<1>();
		double x1 = point.get<0>();
		double x2 = l.point.get<0>();
		double y1 = point.get<1>();
		double y2 = l.point.get<1>();

		double alpha = 1 / (w1 * v2 - w2 * v1)
				* (w2 * (x1 - x2) + w1 * (y2 - y1));

		double x = alpha * v1 + x1;
		double y = alpha * v2 + y1;
		bpoint intersection_point(x, y);
		intersect_pt = intersection_point;
		return true;
	}

	bool intersection_point(bg::model::segment<bpoint> seg,
			bpoint &intersect_pt) {

		Line l(seg);
		bool found = intersection_point(l, intersect_pt);

		if (!found)
			return false;

		//std::cout << "(" << intersect_pt[0] << "," << intersect_pt[1] << ")";

		double pt1_x = seg.first.get<0>();
		double pt1_y = seg.first.get<1>();
		double pt2_x = seg.second.get<0>();
		double pt2_y = seg.second.get<1>();

		double seg_low_y = (pt1_y < pt2_y) ? pt1_y : pt2_y;
		double seg_high_y = (pt1_y > pt2_y) ? pt1_y : pt2_y;
		double seg_low_x = (pt1_x < pt2_x) ? pt1_x : pt2_x;
		double seg_high_x = (pt1_x > pt2_x) ? pt1_x : pt2_x;

		if (is_segment_vertical(seg)){
			return seg_low_y <= intersect_pt.get<1>()
					&& intersect_pt.get<1>() <= seg_high_y;
		}
		else{
			return seg_low_x <= intersect_pt.get<0>()
					&& intersect_pt.get<0>() <= seg_high_x;
		}
	}

	bpoint get_direction_vector() {
		return vec;
	}

	bool is_vertical() {
		bpoint direction_vec(0, 1);
		Line vert_line;
		vert_line.vec = direction_vec;
		return is_parallel(vert_line);
	}

	int eval(double x, double *y) {
		if (is_vertical())
			return RETURN_FAILURE;
		*y = slope * x + intercept;
		return RETURN_OK;
	}

	int get_segment(double x1, double x2, bg::model::segment<bpoint> *seg) {

		if (is_vertical())
			return RETURN_FAILURE;

		double y1, y2;
		eval(x1, &y1);
		eval(x2, &y2);
		bpoint p1(x1, y1);
		bpoint p2(x2, y2);
		std::pair<bpoint, bpoint> seg_points(p1, p2);
		seg->swap(seg_points);
		return RETURN_OK;
	}

	int get_vertical_segment(double y1, double y2,
			bg::model::segment<bpoint> *seg) {

		if (!is_vertical())
			return RETURN_FAILURE;

		double x = point.get<0>();
		bpoint p1(x, y1);
		bpoint p2(x, y2);
		std::pair<bpoint, bpoint> seg_points(p1, p2);
		seg->swap(seg_points);
		return RETURN_OK;
	}

private:
	bpoint point;
	bpoint vec;

	//Undefined if Line is vertical!
	double slope;
	double intercept;

	static bool is_segment_vertical(bg::model::segment<bpoint> &seg) {
		bpoint p1 = seg.first;
		bpoint p2 = seg.second;
		double x_diff = p1.get<0>() - p2.get<0>();
		return (x_diff == 0.0);
	}
};
