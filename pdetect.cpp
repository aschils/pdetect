/*
 * Author: Arnaud Schils, Simon Lardinois, Catholic University of Louvain, 2015
 */

#include "model/LaplaceSolver.hpp"

int main() {
    LaplaceSolver<2> laplace_problem_2d;
    laplace_problem_2d.run ();
	return 0;
}
