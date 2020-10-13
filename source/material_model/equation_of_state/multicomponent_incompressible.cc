/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>
#include <aspect/utilities.h>
// #include <aspect/adiabatic_conditions/interface.h>
// #include <aspect/simulator_access.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      template <int dim>
      void
      MulticomponentIncompressible<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               const unsigned int input_index,
               MaterialModel::EquationOfStateOutputs<dim> &out) const 
      {


        for (unsigned int c=0; c < out.densities.size(); ++c)        
          {        
            // if(use_density_model1)
            //  {
              // std::cout<<"#Phase = "<<c<<" , Ref_compress = "<<reference_compressibilities[c]<<" , Ref temperature = "<<reference_temperatures[c]<<std::endl;

              const double ak = thermal_expansivities[c]/reference_compressibilities[c];
              const double f = 1 + (in.pressure[input_index] - ak*(in.pressure[input_index] - reference_temperatures[c])) *isothermal_bulk_modulus_pressure_derivatives[c] *reference_compressibilities[c];
              const double term1 = 1. + (in.pressure[input_index] - ak*(in.pressure[input_index] - reference_temperatures[c])); 
              out.densities[c] = densities[c]*std::pow(f, 1./isothermal_bulk_modulus_pressure_derivatives[c]);
              out.thermal_expansion_coefficients[c] = thermal_expansivities[c] / f;
              out.specific_heat_capacities[c] = (specific_heats[c] + (in.temperature[input_index] *thermal_expansivities[c] * ak * std::pow(f, -1.-(1./isothermal_bulk_modulus_pressure_derivatives[c]))/ densities[c]));
              out.compressibilities[c] = reference_compressibilities[c]/f;
              out.entropy_derivative_pressure[c] = 0.;
              out.entropy_derivative_temperature[c] = 0.;  

                            //  std::cout<<"#Phase = "<<c<<" , Out compressibility = "<<out.compressibilities[c]<<" , Out density = "<<out.densities[c]<<std::endl;  
                                                        //  std::cout<<"#Phase = "<<c<<" , ak = "<<ak<<" , f = "<<f<<std::endl;
                                                        std::cout<<"#Phase = "<<c<<" , term1 = "<<term1<<" , pressure= "<<in.pressure[input_index]<<" , index = "<<input_index<<std::endl;   

                      
          }          
      }

      template <int dim>
      bool
      MulticomponentIncompressible<dim>::
      is_compressible () const
      {
          return true; 

      }



      template <int dim>
      void
      MulticomponentIncompressible<dim>::declare_parameters (ParameterHandler &prm,
                                                             const double default_thermal_expansion)
      {


        prm.declare_entry ("Reference temperatures", "600",
                          Patterns::Anything(),
                          "List of reference temperatures $T_0$ for background mantle and compositional fields,"
                          "for a total of N+1 values, where N is the number of compositional fields."
                          "If only one value is given, then all use the same value. Units: \\si{\\kelvin}.");
        prm.declare_entry ("Compressibility", "4e-12",
                          Patterns::Anything(),
                          "List of isothermal compressibilities for background mantle and compositional fields,"
                          "for a total of N+1 values, where N is the number of compositional fields."
                          "If only one value is given, then all use the same value. "
                          "Units: \\si{\\per\\pascal}."); 
        prm.declare_entry ("Isothermal bulk modulus pressure derivatives", "4.",
                          Patterns::Anything(),
                          "List of isothermal pressure derivatives of the bulk moduli for background mantle and compositional fields,"
                          "for a total of N+1 values, where N is the number of compositional fields."
                          "If only one value is given, then all use the same value. "
                          "Units: [].");                                                              
                       

        prm.declare_entry ("Densities", "3300.",
                           Patterns::Anything(),
                           "List of densities for background mantle and compositional fields,"
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        prm.declare_entry ("Thermal expansivities", std::to_string(default_thermal_expansion),
                           Patterns::Anything(),
                           "List of thermal expansivities for background mantle and compositional fields,"
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value. Units: \\si{\\per\\kelvin}.");
        prm.declare_entry ("Heat capacities", "1250.",
                           Patterns::Anything(),
                           "List of specific heats $C_p$ for background mantle and compositional fields,"
                           "for a total of N+M+1 values, where N is the number of compositional fields and M is the number of phases. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
        prm.declare_alias ("Heat capacities", "Specific heats");

      }



      template <int dim>
      void
      MulticomponentIncompressible<dim>::parse_parameters (ParameterHandler &prm,
                                                           const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {

        // Establish that a background field is required here
        const bool has_background_field = true;

        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        

        reference_temperatures = Utilities::parse_map_to_double_array (prm.get("Reference temperatures"),
                                                                      list_of_composition_names,
                                                                      has_background_field,
                                                                      "Reference temperatures",
                                                                      true,
                                                                      expected_n_phases_per_composition);
        reference_compressibilities = Utilities::parse_map_to_double_array (prm.get("Compressibility"),
                                                                            list_of_composition_names,
                                                                            has_background_field,
                                                                            "Compressibility",
                                                                            true,
                                                                            expected_n_phases_per_composition);
        isothermal_bulk_modulus_pressure_derivatives = Utilities::parse_map_to_double_array (prm.get("Isothermal bulk modulus pressure derivatives"),
                                                                                            list_of_composition_names,
                                                                                            has_background_field,
                                                                                            "Isothermal bulk modulus pressure derivatives",
                                                                                            true,
                                                                                            expected_n_phases_per_composition);                                                                            
        densities = Utilities::parse_map_to_double_array (prm.get("Densities"),
                                                          list_of_composition_names,
                                                          has_background_field,
                                                          "Densities",
                                                          true,
                                                          expected_n_phases_per_composition);

        thermal_expansivities = Utilities::parse_map_to_double_array (prm.get("Thermal expansivities"),
                                                                      list_of_composition_names,
                                                                      has_background_field,
                                                                      "Thermal expansivities",
                                                                      true,
                                                                      expected_n_phases_per_composition);

        specific_heats = Utilities::parse_map_to_double_array (prm.get("Heat capacities"),
                                                               list_of_composition_names,
                                                               has_background_field,
                                                               "Specific heats",
                                                               true,
                                                               expected_n_phases_per_composition);
  
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
  template class MulticomponentIncompressible<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
