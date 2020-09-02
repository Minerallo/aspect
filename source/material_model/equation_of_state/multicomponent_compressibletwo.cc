/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/utilities.h>

#include <aspect/adiabatic_conditions/interface.h>


// #include <aspect/material_model/interface.h>
// #include <aspect/simulator_access.h>
// #include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/equation_of_state/multicomponent_compressibletwo.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      template <int dim>
      void
      MulticomponentCompressibletwo<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               MaterialModel::EquationOfStateOutputs<dim> &out) const
      {

      for (unsigned int q=0; q < in.n_evaluation_points(); ++q)
        {
          const Point<dim> position = in.position[q];
          const double pressure = in.pressure[q];
          const double temperature = std::max(in.temperature[q], 1.); // temperature can't be zero for correct evaluation

          for (unsigned int c=0; c < out.densities.size(); ++c)
            {

              out.specific_heat_capacities[q] = isochoric_specific_heats[c];
              // out.thermal_conductivities[i] = k_value;
              out.thermal_expansion_coefficients[q] = reference_thermal_expansivities[c];

              double rho = reference_densities[c] * std::exp(reference_compressibility * (pressure - this->get_surface_pressure()));
              rho *= (1 - reference_thermal_expansivities[c] * (temperature - this->get_adiabatic_conditions().temperature(position)));

              out.densities[q] = rho;
              out.compressibilities[q] = reference_compressibility; // 1/rho drho/dp
              out.entropy_derivative_pressure[q] = 0.0;
              out.entropy_derivative_temperature[q] = 0.0;
              // Change in composition due to chemical reactions at the
              // given positions. The term reaction_terms[i][c] is the
              // change in compositional field c at point i.
              // for (unsigned int c=0; c<in.composition[i].size(); ++c)
              //   out.reaction_terms[i][c] = 0.0;
          }
        }
      }



      template <int dim>
      bool
      MulticomponentCompressibletwo<dim>::
      is_compressible () const
      {
        return true;
      }



      template <int dim>
      void
      MulticomponentCompressibletwo<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Densities", "3300.",
                           Patterns::Anything(),
                           "List of densities for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value.  Units: $kg / m^3$");
        prm.declare_entry ("Heat capacities", "1250.",
                           Patterns::Anything(),
                           "List of specific heats $C_p$ for background mantle and compositional fields,"
                           "for a total of N+1 values, where N is the number of compositional fields."
                           "If only one value is given, then all use the same value. Units: $J /kg /K$");
        prm.declare_entry ("Thermal expansivities", "2e-5",
                            Patterns::Anything(),
                            "List of thermal expansivities for background mantle and compositional fields,"
                            "for a total of N+1 values, where N is the number of compositional fields."
                            "If only one value is given, then all use the same value. Units: $1/K$");
        prm.declare_entry ("Reference compressibility", "4e-12",
                            Patterns::Double (0.),
                            "The value of the reference compressibility. "
                            "Units: $1/Pa$.");                        
      }



      template <int dim>
      void
      MulticomponentCompressibletwo<dim>::parse_parameters (ParameterHandler &prm,
                                                           const std::shared_ptr<std::vector<unsigned int>> expected_n_phases_per_composition)
      {
        // // Establish that a background field is required here
        const bool has_background_field = true;

        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        reference_densities = Utilities::parse_map_to_double_array (prm.get("Densities"),
                                                          list_of_composition_names,
                                                          has_background_field,
                                                          "Densities",
                                                          true,
                                                          expected_n_phases_per_composition);

          // k_value                    = prm.get_double ("Thermal conductivity");
        isochoric_specific_heats = Utilities::parse_map_to_double_array (prm.get("Heat capacities"),
                                                               list_of_composition_names,
                                                               has_background_field,
                                                               "Heat capacities",
                                                               true,
                                                               expected_n_phases_per_composition);

        reference_thermal_expansivities = Utilities::parse_map_to_double_array (prm.get("Thermal expansivities"),
                                                                      list_of_composition_names,
                                                                      has_background_field,
                                                                      "Thermal expansivities",
                                                                      true,
                                                                      expected_n_phases_per_composition);

          reference_compressibility  = prm.get_double ("Reference compressibility");
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
#define INSTANTIATE(dim) \
  template class MulticomponentCompressibletwo<dim>;
      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}
