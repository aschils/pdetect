/*
 * utils.hpp

 *
 *  Created on: 23 nov. 2015
 *      Author: aschils
 */

#include <stdlib.h>

class Utils {

public:
	static bool is_same_double(double a, double b, double epsilon){
		return abs(a-b) <= epsilon;
	}
};
