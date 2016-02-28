/*
 * MyGeometryInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include "../geometry_info/BasicGeometry.hpp"

class MyGeometryInfo {

public:

	/**
	 * Note: a point inside a strip is not considered as being
	 * part of the geometry, since a strip is a hole in the
	 * geometry
	 */
	virtual bool is_point_inside_geometry(Point<2> p) = 0;

	/* Returns the intersection points between the particle trajectory
	 * and the domain boundaries
	 * intersections are placed in increasing x position
	 */
	std::vector<Point<2>> boundaries_intersections(Line line) {

		std::vector<Segment> segments = get_geometry_segments();

		//Compute intersections with all segments of detector boundary
		std::vector<Point<2>> intersections;

		for (unsigned i = 0; i < segments.size(); i++) {
			add_intersection(line, segments[i], intersections);
		}

		Utils::sort_points_by_coord<2>(&intersections);
		intersections.erase(unique(intersections.begin(), intersections.end()),
				intersections.end());
		intersections.shrink_to_fit();
		return intersections;
	}

	//Return the maximum width of the domain
	virtual unsigned get_width() = 0;

	//Return the maximum length of the domain
	virtual unsigned get_length() = 0;

	virtual unsigned get_nbr_of_strips() = 0;

	virtual unsigned get_dimension() = 0;

	virtual Line get_mid_length_vertical_line() = 0;

	virtual unsigned get_strip_width() = 0;

	~MyGeometryInfo() {
	}

protected:

	virtual std::vector<Segment> get_geometry_segments() = 0;

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

			bool before_in_det = is_point_inside_geometry(before);
			bool after_in_det = is_point_inside_geometry(after);

			bool is_a_cross_boundary_point = !((before_in_det && after_in_det)
					|| (!before_in_det && !after_in_det));

			if (is_a_cross_boundary_point) {
				intersections.push_back(intersect);
			}
		}
	}
};
