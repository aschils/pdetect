/*
 * BasicGeometry.hpp
 *
 *  Created on: 12 févr. 2016
 *      Author: aschils
 */

#pragma once

using namespace dealii;

#include "../errors.hpp"

class Segment {
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
};

class Line {

public:

	Line(){}

	Line(double slope, double intercept) {
		Point<2> p1(0.0, intercept);
		double y = slope + intercept;
		Point<2> p2(1, y);
		Segment seg(p1, p2);
		constructor_from_segment(seg);
		//std::cout << "Point (" << point[0] << "," << point[1] << ")"
		//		<< "Vec (" << vec[0] << "," << vec[1] << ")" << std::endl;
	}

	Line(Segment seg) {
		constructor_from_segment(seg);
	}

	Line(Point<2> p1, Point<2> p2){
		Segment seg(p1,p2);
		constructor_from_segment(seg);
	}

	void constructor_from_segment(Segment seg) {
		point = seg.p1;
		double norm = seg.compute_length();
		double vec_x = (seg.p2[0] - seg.p1[0]) / norm;
		double vec_y = (seg.p2[1] - seg.p1[1]) / norm;
		Point<2> vec(vec_x, vec_y);
		this->vec = vec;
	}

	bool is_parallel(Line l) {
		double epsilon = 0.00001;
		bool same_vec_x = Utils::equals_double(vec[0], l.vec[0], epsilon);
		bool same_vec_y = Utils::equals_double(vec[1], l.vec[1], epsilon);
		return same_vec_x && same_vec_y;
	}

	bool intersection_point(Line l, Point<2> &intersect_pt) {

		if (is_parallel(l)) {
			return false;
		}

		//Well it's geometry: intersection between two lines described
		//by parametric equations
		double v1 = vec[0];
		double v2 = vec[1];
		double w1 = l.vec[0];
		double w2 = l.vec[1];
		double x1 = point[0];
		double x2 = l.point[0];
		double y1 = point[1];
		double y2 = l.point[1];

		double alpha = 1 / (w1 * v2 - w2 * v1)
				* (w2 * (x1 - x2) + w1 * (y2 - y1));

		double x = alpha * v1 + x1;
		double y = alpha * v2 + y1;
		Point<2> intersection_point(x, y);
		intersect_pt = intersection_point;
		return true;
	}

	bool intersection_point(Segment seg, Point<2> &intersect_pt) {

		Line l(seg);
		bool found = intersection_point(l, intersect_pt);

		if (!found)
			return false;

		//std::cout << "(" << intersect_pt[0] << "," << intersect_pt[1] << ")";

		double seg_low_y = seg.get_low_y();
		double seg_high_y = seg.get_high_y();
		double seg_low_x = seg.get_low_x();
		double seg_high_x = seg.get_high_x();

		if (seg.is_vertical())
			return seg_low_y <= intersect_pt[1] && intersect_pt[1] <= seg_high_y;
		else
			return seg_low_x <= intersect_pt[0] && intersect_pt[0] <= seg_high_x;
	}

	Point<2> get_direction_vector(){
		return vec;
	}

	bool is_vertical(){
		Point<2> p1(0,0);
		Point<2> p2(0,1);
		Line vert_line(p1,p2);
		return is_parallel(vert_line);
	}

private:
	Point<2> point;
	Point<2> vec;
};
