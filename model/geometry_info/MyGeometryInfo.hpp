/*
 * MyGeometryInfo.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: aschils
 */

#pragma once

#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>

#include "../geometry_info/BasicGeometry.hpp"

namespace bg = boost::geometry;

class MyGeometryInfo {

public:

	typedef bg::model::point<double, 2, bg::cs::cartesian> bpoint;

	/**
	 * Note: a point inside a strip is not considered as being
	 * part of the geometry, since a strip is a hole in the
	 * geometry
	 */
	virtual bool is_point_inside_geometry(bpoint p) = 0;

	/* Returns the intersection points between the particle trajectory
	 * and the domain boundaries
	 * intersections are placed in increasing x position
	 */
	std::vector<bpoint> boundaries_intersections(Line line) {

		std::vector<bg::model::segment<bpoint> > segments = get_geometry_segments();
		//std::vector<bg::model::segment<bpoint> > segments = get_geometry_segments();

		//Compute intersections with all segments of detector boundary
		std::vector<bpoint> intersections;

		for (unsigned i = 0; i < segments.size(); i++) {
			add_intersection(line, segments[i], intersections);
		}

		Utils::sort_points_by_coord(&intersections);
		intersections = Utils::unique_bpoints(intersections);
		//unique(intersections.begin(), intersections.end());
		//intersections.erase(unique(intersections.begin(), intersections.end()),
		//		intersections.end());
		//intersections.shrink_to_fit();
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

	virtual std::vector<bg::model::segment<bpoint> > get_geometry_segments() = 0;

	bool intersect_segment(Line line, bg::model::segment<bpoint> seg,
			bpoint &intersect_pt) {
		return line.intersection_point(seg, intersect_pt);
	}

	void add_intersection(Line line, bg::model::segment<bpoint> seg,
			std::vector<bpoint> &intersections) {
		bpoint intersect;
		if (intersect_segment(line, seg, intersect)) {

			//We must remove the intersection point
			//if at + and -epsilon not in detector OR
			// if at + and - epsilon in detector

			const double epsilon = 0.001;

			//bpoint temp = line.get_direction_vector();
			//bpoint eps_direction_vec = Utils::mul_coords_by_const(
			//		temp, epsilon);
			bpoint eps_direction_vec = line.get_direction_vector();
			bg::multiply_value(eps_direction_vec, epsilon);
			bpoint before(intersect.get<0>() - eps_direction_vec.get<0>(),
					intersect.get<1>() - eps_direction_vec.get<1>());
			bpoint after = intersect;
			bg::add_point(after, eps_direction_vec);

			bool before_in_det = is_point_inside_geometry(before);
			bool after_in_det = is_point_inside_geometry(after);

			bool is_a_cross_boundary_point = !((before_in_det && after_in_det)
					|| (!before_in_det && !after_in_det));

			if (is_a_cross_boundary_point) {
				intersections.push_back(intersect);
			}
		}
	}

	/**
	 * Get segment of particle trajectory line that may intersects the
	 * detector boundaries. Required because boost geometry library
	 * does not compute intersections between a line and a segment
	 * (only between segments).
	 */
	bg::model::segment<bpoint> segment_from_line(Line &line){

		bg::model::segment<bpoint> seg;

		if(line.is_vertical())
			line.get_vertical_segment(-1.0, get_width()+1.0, &seg);
		else
			line.get_segment(-1.0, get_length()+1.0, &seg);
		return seg;
	}

};
