#ifndef __GRADIENT_HPP__
#define __GRADIENT_HPP__

#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

template<unsigned dim>
class Derivatives: public DataPostprocessor<dim> {

public:

	Derivatives();

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

template<unsigned dim>
Derivatives<dim>::Derivatives(){}

template<unsigned dim>
UpdateFlags Derivatives<dim>::get_needed_update_flags() const {
	return update_values | update_gradients | update_q_points |
			update_second_derivatives;
}

template<unsigned dim>
void Derivatives<dim>::compute_derived_quantities_scalar(
		const std::vector<double> &uh, const std::vector<Tensor<1, dim> > &duh,
		const std::vector<Tensor<2, dim> > &dduh,
		const std::vector<Point<dim> > &normals,
		const std::vector<Point<dim> > &evaluation_points,
		std::vector<Vector<double> > &computed_quantities) const {

	const unsigned int n_quadrature_points = uh.size();

	for (unsigned int q = 0; q < n_quadrature_points; ++q) {
		for (unsigned int d = 0; d < dim; ++d)
			computed_quantities[q](d) = duh[q][d];

		for(unsigned int i = 0; i<dim; i++){
			for(unsigned j = 0; j<dim; j++)
				computed_quantities[q](dim+i*dim+j) = dduh[q][i][j];

		}
	}
}

template<unsigned dim>
std::vector<std::string> Derivatives<dim>::get_names() const {
	std::vector<std::string> solution_names(dim, "gradient");
	for(unsigned i=0; i<(dim*dim); i++)
		solution_names.push_back("seconde_derivatives");
	return solution_names;
}

template<unsigned dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Derivatives<
		dim>::get_data_component_interpretation() const {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation (dim+dim*dim,
                    DataComponentInterpretation::component_is_part_of_vector);
	return interpretation;
}

#endif
