/*
 * VectorUtils.hpp
 *
 *  Created on: 4 déc. 2015
 *      Author: aschils
 */

#ifndef __VECTOR_UTILS_HPP__
#define __VECTOR_UTILS_HPP__

class VectorUtils {

public:
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
			double lambda){
		std::vector<double> res(v.size());
		std::transform(v.begin(), v.end(), res.begin(),
		               std::bind1st(std::multiplies<double>(),lambda));
		return res;
	}

	static std::vector<double> opposite_vector(std::vector<double> v){
		return multiply_vector_by_scalar(v, -1.0);
	}
};
#endif
