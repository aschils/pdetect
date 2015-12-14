/*
 * StraightLine.cpp
 *
 *  Created on: 14 déc. 2015
 *      Author: slardinois
 */

template<unsigned dim>
StraightLine<dim>::StraightLine(double alpha, Point<dim> const &pass, 
								Solution<dim> *sol,
								double precision) {

	this->precision = precision;
	this->sol = sol;
	int l = sol->values_at_cells.size() - 1;
	rect_length_fe = sol->values_at_cells[l].first->vertex(3)[0];
	rect_width_fe = sol->values_at_cells[l].first->vertex(3)[1];

	construct_line(alpha, pass);
}

template<unsigned dim>
std::pair<Point<dim>, Point<dim>> StraightLine<dim>::get_extremeties(double alpha, 
														Point<dim> const &pass) {
	Point<dim> begin = pass;
	Point<dim> end = pass;
	if(alpha == 0) { //Horizontal line
		begin[0] = 0;
		end[0] = rect_length_fe;
	}
	else if(alpha == 90) { //Vertical line
		begin[1] = 0;
		end[1] = rect_width_fe;
	}

	std::pair<Point<dim>, Point<dim>> extrem;
	extrem.first = begin;
	extrem.second = end;

	return extrem;
}

template <unsigned dim>
void StraightLine<dim>::construct_line(double alpha, Point<dim> const &pass) { 

	std::pair<Point<dim>, Point<dim>> extrem = get_extremeties(alpha, pass);

	Point<dim> point = extrem.first;
	while(point != extrem.second) {
		if(alpha == 0) {
			point[0] = point[0] + precision;
			if(point[0] >= extrem.second[0])
				point = extrem.second;
		}
		else if(alpha == 90) {
			point[1] = point[1] + precision;
			if(point[1] >= extrem.second[1])
				point = extrem.second;
		}

		std::cout << point[0] << std::endl
				<< point[1] << std::endl;
		ValuesAtPoint<dim> value = sol->get_values(point);

		std::pair<ValuesAtPoint<dim>, Point<dim>> values_at_point;

		values_at_point.first = value;
		values_at_point.second = point;

		values_on_line.push_back(values_at_point);
	}
}