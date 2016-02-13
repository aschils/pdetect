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

	double compute_length() {
		double diff_x = p2[0] - p1[0];
		double diff_y = p2[1] - p1[1];
		return std::sqrt(diff_x * diff_x + diff_y * diff_y);
	}

private:
	bool vertical;
	double slope;
};

/*
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
 */

//class Line {
//
//public:
//
//	Line(double slope, double intercept) {
//		non_vert_line.set_slope(slope);
//		non_vert_line.set_intercept(intercept);
//		vertical = false;
//	}
//
//	/**
//	 * segment (p1,p2)
//
//	 vecteur directeur = p2-p1
//
//	 eq param: x = x_1 + alpha (x_2 - x_1)
//	 y = y_1 + alpha(y_2 - y_1)
//
//	 alpha = (x - x_1)/(x_2-x_1)
//
//	 y(x) = y_1 + (x-x_1)/(x_2-x_1) (y_2-y_1) = y_1 - x_1 (y_2 - y_1) / (x_2-x_1)
//	 + x (y_2 - y_1) / (x_2-x_1)
//
//	 => intercept = y_1 - x_1 (y_2 - y_1) / (x_2-x_1)
//	 = (y_1 x_2 - y_1 x_1 - x_1 y_2 + x_1 y_1)/(x_2 - x_1)
//	 = (y_1 x_2 - x_1 y_2) / (x_2 - x_1)
//
//	 *
//	 */
//	Line(Segment seg) {
//		vertical = seg.is_vertical();
//
//		if (vertical) {
//			vert_line.set_x(seg.p1[0]);
//		} else {
//			non_vert_line.set_slope(seg.get_slope());
//			double x1 = seg.p1[0];
//			double x2 = seg.p2[0];
//			double y1 = seg.p1[1];
//			double y2 = seg.p2[1];
//			double intercept = (y1 * x2 - x1 * y2) / (x2 - x1);
//			non_vert_line.set_intercept(intercept);
//		}
//	}
//
//	bool intersection_point(Line l, Segment seg, Point<2> &intersect){
//
//	}
//
//	bool intersection_point(Line l, Point<2> &intersect) {
//
//		if (l.vertical && vertical)
//			return false;
//
//		if (vertical) {
//			return intersection_point(vert_line, l.non_vert_line);
//		} else if(l.vertical) {
//			return intersection_point(l.vert_line, non_vert_line);
//		}
//		else{
//			return intersection_point(non_vert_line, l.non_vert_line);
//		}
//	}
//
//private:
//	bool vertical;
//	VerticalLine vert_line;
//	NonVerticalLine non_vert_line;
//
//	Point<2> intersection_point(VerticalLine v, NonVerticalLine l) {
//		Point<2> intersect(v.get_x(), l.eval(v.get_x()));
//		//std::cout << "intersection_point(VerticalLine v, NonVerticalLine l)" << std::endl;
//		//std::cout << "x=" << intersect[0] << " y=" << intersect[1] << std::endl;
//		return intersect;
//	}
//
//	Point<2> intersection_point(NonVerticalLine l1, NonVerticalLine l2) {
//		double x = (l1.get_intercept() - l2.get_intercept())
//				/ (l1.get_slope() - l2.get_slope());
//		double y = l1.eval(x);
//		Point<2> intersection_point(x, y);
//		return intersection_point;
//	}
//};
class Line {

public:

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
			std::cout << "parallel" << std::endl;
			return false;
		}

		//Well it's math
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

		//std::cout << "alpha " << alpha << " v1 " << v1 << " x1 " << x1 << std::endl;

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

private:
	Point<2> point;
	Point<2> vec;
};
