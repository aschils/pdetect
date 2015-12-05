/*
 * VectorUtils.hpp
 *
 *  Created on: 4 déc. 2015
 *      Author: aschils
 */

#ifndef __VECTOR_UTILS_HPP__
#define __VECTOR_UTILS_HPP__

#include "Utils.hpp"

class VectorUtils {

public:

	/**
	 * First element of the pair is supposed to be vector with dim elements.
	 * v is sorted in ascending order regarding the std::pair(coord, val)
	 *
	 * coordx = (x_1,...,x_n), coordy = (y_1,...,y_n)
	 *
	 * coordx < coordy if x_n < y_n
	 *  or x_n = y_n and x_(n-1) < y_(n-1)
	 *  or x_n = y_n and x_(n-1) = y_(n-1) and x_(n-2) < y_(n-2)
	 *  and so on.
	 *
	 *
	 */
	template <unsigned dim, typename T>
	static void sort_by_coord(
			std::vector<std::pair<std::vector<double>, T> > &v) {

		auto cmp =
				[](std::pair<std::vector<double>,
						T> const & a,
						std::pair<std::vector<double>,
						T> const & b) {
					double epsilon = 0.00001;
					for(int i=dim-1; i>=0; i--) {
						if(!Utils::greater_than_or_equals_double((a.first)[i], (b.first)[i], epsilon))
						return true;
						else if(!Utils::less_than_or_equals_double(a.first[i], b.first[i], epsilon))
						return false;
					}
					return false;
				};

		std::sort(v.begin(), v.end(), cmp);
	}

	static void print_vec_components(std::vector<double> vec) {
		for (unsigned j = 0; j < vec.size(); j++) {
			std::cout << vec[j];

			if (j != vec.size() - 1)
				std::cout << ",";
		}
	}

	static void print_vec_of_pair_of_vec(
			std::vector<std::pair<std::vector<double>, std::vector<double> > > vec) {
		for (unsigned i = 0; i < vec.size(); i++) {

			std::cout << "[(";
			print_vec_components(vec[i].first);
			std::cout << ") (";
			print_vec_components(vec[i].second);
			std::cout << ")] ";
		}
	}

	static std::vector<double> multiply_vector_by_scalar(std::vector<double> v,
			double lambda) {
		std::vector<double> res(v.size());
		std::transform(v.begin(), v.end(), res.begin(),
				std::bind1st(std::multiplies<double>(), lambda));
		return res;
	}

	static std::vector<double> opposite_vector(std::vector<double> v) {
		return multiply_vector_by_scalar(v, -1.0);
	}
};
#endif
