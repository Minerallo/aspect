// /*
//   Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

//   This file is part of ASPECT.

//   ASPECT is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2, or (at your option)
//   any later version.

//   ASPECT is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with ASPECT; see the file LICENSE.  If not see
//   <http://www.gnu.org/licenses/>.
// */
// #include <aspect/postprocess/visualization.h>
// #include <aspect/simulator_access.h>
// #include <deal.II/numerics/data_postprocessor.h>

// #include <aspect/postprocess/visualization/stress.h>


// namespace aspect
// {
//   namespace Postprocess
//   {
//     namespace VisualizationPostprocessors
//     {
//       /**
//        * A class derived from DataPostprocessor that takes an output vector
//        * and computes a variable that represents the 3 or 6 independent
//        * components (in 2d and 3d, respectively) of the stress tensor at every
//        * point. The shear stress is defined as $2 \eta (\varepsilon(\mathbf u)
//        * - \tfrac 13 \textrm{trace}\ \varepsilon(\mathbf u) \mathbf 1) +pI =
//        * 2\eta (\varepsilon(\mathbf u) - \frac 13 (\nabla \cdot \mathbf u)
//        * \mathbf I) + pI$.  The second term in the parentheses is zero if the
//        * model is incompressible.
//        *
//        * The member functions are all implementations of those declared in the
//        * base class. See there for their meaning.
//        */
//       template <int dim>
//       class StressNorm
//         : public DataPostprocessor<dim>,
//           public SimulatorAccess<dim>,
//           public Interface<dim>
//       {
//         public:
//           StressNorm ();

//           void
//           evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
//                                 std::vector<Vector<double> > &computed_quantities) const override;

//       };  
//     }
//   }
// }

// namespace aspect
// {
//   namespace Postprocess
//   {
//     namespace VisualizationPostprocessors
//     {
//       template <int dim>
//       StressNorm<dim>::
//       StressNorm ()
//         :
//         DataPostprocessorScalar<dim> ("stress_norm",
//                                       update_gradients | update_quadrature_points)
//       {}

//       template <int dim>
//       void
//       StressNorm<dim>::
//       evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
//                             std::vector<Vector<double> > &computed_quantities) const
//       {
//         const unsigned int n_quadrature_points = input_data.solution_values.size();
//         Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
//         Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
//         Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());
//         Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,          ExcInternalError());

//         MaterialModel::MaterialModelInputs<dim> in(input_data,
//                                                    this->introspection());
//         MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
//                                                      this->n_compositional_fields());

//         // Compute the viscosity...
//         this->get_material_model().evaluate(in, out);

//         // ...and use it to compute the StressNorme
//         for (unsigned int q=0; q<n_quadrature_points; ++q)
//           {
//             const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
//             const SymmetricTensor<2,dim> compressible_strain_rate
//               = (this->get_material_model().is_compressible()
//                  ?
//                  strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
//                  :
//                  strain_rate);

//             const double eta = out.viscosities[q];

//             const SymmetricTensor<2,dim> stress = 2*eta*compressible_strain_rate +
//                                                   in.pressure[q] ;
//                                                   // * unit_symmetric_tensor<dim>();

//        //   Calculation of the StressNorm , similar to Babeyko et al., 2006 Slim2D

//         // stress_norm = sqrt(pow((SymmetricTensor[1][1]-(SymmetricTensor[2][2]),2)+4*pow(SymmetricTensor[1][2],2))

//             computed_quantities[q](0) = stress;
//             // std::sqrt(stress*stress);
//           }

//         // average the values if requested
//         const auto &viz = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Visualization<dim> >();
//         if (!viz.output_pointwise_stress_and_strain())
//           average_quantities(computed_quantities);

//        }

//     }
//   }
// }


// // explicit instantiations
// namespace aspect
// {
//   namespace Postprocess
//   {
//     namespace VisualizationPostprocessors
//     {
//       ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StressNorm,
//                                                   "stress norm",
//                                                   "A visualization output object that generates output "
//                                                   "for the 3 (in 2d) or 6 (in 3d) components of the stress "
//                                                   "tensor, i.e., for the components of the tensor "
//                                                   "$2\\eta\\varepsilon(\\mathbf u)+pI$ "
//                                                   "in the incompressible case and "
//                                                   "$2\\eta\\left[\\varepsilon(\\mathbf u)-"
//                                                   "\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I\\right]+pI$ "
//                                                   "in the compressible case.")
//     }
//   }
// }