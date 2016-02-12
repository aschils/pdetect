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
		vertical = x_diff == 0.0;

		if (!vertical)
			slope = p1[1] - p2[1] / x_diff;
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

private:
	bool vertical;
	double slope;
};

class NonVerticalLine {

public:

	NonVerticalLine() {
	}

	NonVerticalLine(double slope, double intercept) {
		this->slope = slope;
		this->intercept = intercept;
	}

	double eval(double x) {
		return intercept + slope * x;
	}

	double inverse(double y) {
		return (y - intercept) / slope;
	}

	double get_intercept() {
		return intercept;
	}

	double get_slope() {
		return slope;
	}

	void set_intercept(double intercept) {
		this->intercept = intercept;
	}

	void set_slope(double slope) {
		this->slope = slope;
	}

private:
	double slope = 0.0;
	double intercept = 0.0;

};

class VerticalLine {

public:

	VerticalLine() {
	}

	VerticalLine(double x) {
		this->x = x;
	}

	double inverse(double y) {
		return x;
	}

	double get_x() {
		return x;
	}

	void set_x(double x) {
		this->x = x;
	}

private:
	double x;
};

class Line {

public:

	Line(double slope, double intercept) {
		non_vert_line.set_slope(slope);
		non_vert_line.set_intercept(intercept);
		vertical = false;
	}

	/**
	 * segment (p1,p2)

	 vecteur directeur = p2-p1

	 eq param: x = x_1 + alpha (x_2 - x_1)
	 y = y_1 + alpha(y_2 - y_1)

	 alpha = (x - x_1)/(x_2-x_1)

	 y(x) = y_1 + (x-x_1)/(x_2-x_1) (y_2-y_1) = y_1 - x_1 (y_2 - y_1) / (x_2-x_1)
	 + x (y_2 - y_1) / (x_2-x_1)

	 => intercept = y_1 - x_1 (y_2 - y_1) / (x_2-x_1)
	 = (y_1 x_2 - y_1 x_1 - x_1 y_2 + x_1 y_1)/(x_2 - x_1)
	 = (y_1 x_2 - x_1 y_2) / (x_2 - x_1)

	 *
	 */
	Line(Segment seg) {
		vertical = seg.is_vertical();
		if (vertical) {
			vert_line.set_x(seg.p1[0]);
		} else {
			non_vert_line.set_slope(seg.get_slope());
			double x1 = seg.p1[0];
			double x2 = seg.p2[0];
			double y1 = seg.p1[1];
			double y2 = seg.p2[1];
			double intercept = (y1 * x2 - x1 * y2) / (x2 - x1);
			non_vert_line.set_intercept(intercept);
		}
	}

	Point<2> intersection_point(Line l) {

		if (l.vertical && vertical)
			throw PRECONDITIONS_VIOLATED;

		if (vertical) {
			return intersection_point(vert_line, l.non_vert_line);
		} else {
			return intersection_point(l.vert_line, non_vert_line);
		}
	}

private:
	bool vertical;
	VerticalLine vert_line;
	NonVerticalLine non_vert_line;

	Point<2> intersection_point(VerticalLine v, NonVerticalLine l) {
		Point<2> intersect(v.get_x(), l.eval(v.get_x()));
		return intersect;
	}

	Point<2> intersection_point(NonVerticalLine l1, NonVerticalLine l2) {
		double x = (l1.get_intercept() - l2.get_intercept())
				/ (l1.get_slope() - l2.get_slope());
		double y = l1.eval(x);
		Point<2> intersection_point(x, y);
		return intersection_point;
	}
};
