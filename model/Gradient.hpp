#ifndef __GRADIENT_HPP__
#define __GRADIENT_HPP__

#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

template<int dim>
class Gradient: public DataPostprocessor<dim> {

public:

	Gradient();

	virtual void compute_derived_quantities_scalar(
			const std::vector<double> &uh,
			const std::vector<Tensor<1, dim> > &duh,
			const std::vector<Tensor<2, dim> > &dduh,
			const std::vector<Point<dim> > &normals,
			const std::vector<Point<dim> > &evaluation_points,
			std::vector<Vector<double> > &computed_quantities) const;

	virtual UpdateFlags get_needed_update_flags() const;

	virtual std::vector<std::string> get_names() const;

	virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
	get_data_component_interpretation() const;

};

template<int dim>
Gradient<dim>::Gradient(){}

template<int dim>
UpdateFlags Gradient<dim>::get_needed_update_flags() const {
	return update_values | update_gradients | update_q_points;
}

template<int dim>
void Gradient<dim>::compute_derived_quantities_scalar(
		const std::vector<double> &uh, const std::vector<Tensor<1, dim> > &duh,
		const std::vector<Tensor<2, dim> > &dduh,
		const std::vector<Point<dim> > &normals,
		const std::vector<Point<dim> > &evaluation_points,
		std::vector<Vector<double> > &computed_quantities) const {

	const unsigned int n_quadrature_points = uh.size();

	for (unsigned int q = 0; q < n_quadrature_points; ++q) {
		for (unsigned int d = 0; d < dim; ++d)
			computed_quantities[q](d) = duh[q][d];
		//sqrt(
		//		duh[q][0] * duh[q][0] + duh[q][1] * duh[q][1]);
	}
}

template<int dim>
std::vector<std::string> Gradient<dim>::get_names() const {
	std::vector<std::string> solution_names(dim, "Electric_field");
	return solution_names;
}

template<int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Gradient<
		dim>::get_data_component_interpretation() const {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation (dim,
                    DataComponentInterpretation::component_is_part_of_vector);
	//std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(
	//		1, DataComponentInterpretation::component_is_scalar);
	return interpretation;
}

#endif
