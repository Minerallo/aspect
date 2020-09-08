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
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>


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

//           OPTION 1 : Weak zone constant density, lithosphere reference temperature  in order for the average density to be equal to the density of reference, mantle density constrained by mantle adiabat.
/*    set Densities = background:3300 ,UPM:3300|3736|3875|4263 ,TZ:3736|3875|4263 ,TZ2:3875|4263 ,LM:4263 ,Weak_Layer:2670 ,Oceanic_crust:3000|3450|3567|3900|4112|4256 ,Oceanic_mantle:3280|3737|3938|4382 ,Sediments:2670 ,Upper_Crust:2800 ,Lower_Crust:3000|3450|3567|3900|4112|4256 ,Continental_mantle:3280|3737|3938|4382 ,Craton:3280|3737|3938|4382 ,plastic_strain:3300, Damping:4530
         */ 

//     set Densities = background:3300 ,UPM:3300|3736|3875|4263 ,TZ:3736|3875|4263 ,TZ2:3875|4263 ,LM:4263 ,Oceanic_crust:3000|3450|3567|3900|4112|4256 ,Oceanic_mantle:3280|3737|3938|4382 ,Weak_Zone:2670,Sediments:2670 ,Upper_Crust:2920 ,Lower_Crust:3000|3450|3567|3900|4112|4256 ,Continental_mantle:3280|3737|3938|4382 ,Craton:3260|3737|3938|4382 ,plastic_strain:3300, Damping:4530

//     set Densities = background:0 ,UPM:1|2|3|4 ,TZ:5|6|7 ,TZ2:8|9 ,LM:10 ,Oceanic_crust:11|12|13|14|15|16 ,Oceanic_mantle:17|18|19|20 ,Weak_Zone:21,Sediments:22 ,Upper_Crust:23 ,Lower_Crust:24|25|26|27|28|29 ,Continental_mantle:30|31|32|33 ,Craton:34|35|36|37 ,plastic_strain:38, Damping:39

         std::vector<double> T_ref = {293,293,293,293,293,293,293,293,293,293,293,350,650,800,890,1080,1315,1000,1000,1000,1000,293,326,510,803,803,803,803,803,803,1240,1240,1240,1240,1100,1100,1100,1100,293,293};
          
          
                // Loop through all requested points
        // for (unsigned int c=0; c < out.densities.size(); ++c)
        for (unsigned int c=0; c < out.densities.size(); ++c)        
          {
            //   //for some reason it doesn't work if I just write >=12
            if(c==22 || c==23 || c==24 || c==25 || c==26 || c==27 || c==28 || c==29 || c==30 ||  c==31 ||  c==32 ||  c==33 || c==34 || c==35 || c==36 || c==37 || c==38 || c==39 || c==40)
            {
              out.densities[c] = densities[c] * (1 - thermal_expansivities[c] * (in.temperature[input_index] - T_ref[c]));    
              out.thermal_expansion_coefficients[c] = thermal_expansivities[c];
              out.specific_heat_capacities[c] = specific_heats[c];
              out.compressibilities[c] = 0.0;
              out.entropy_derivative_pressure[c] = 0.0;
              out.entropy_derivative_temperature[c] = 0.0;
            }                    
            //for the rest the reference temperature is based on the adiabatic gradient
            else
             {
              out.densities[c] = densities[c] * (1 - thermal_expansivities[c] * (in.temperature[input_index] - (this->get_adiabatic_conditions().temperature(in.position[input_index]))));
              out.thermal_expansion_coefficients[c] = thermal_expansivities[c];
              out.specific_heat_capacities[c] = specific_heats[c];
              out.compressibilities[c] = 0.0;
              out.entropy_derivative_pressure[c] = 0.0;
              out.entropy_derivative_temperature[c] = 0.0;
             }
          }           
          
// //           OPTION 1 : Weak zone constant density, lithosphere reference temperature  in order for the average density to be equal to the density of reference, mantle density constrained by mantle adiabat.
// // /*    set Densities = background:3300 ,UPM:3300|3736|3875|4263 ,TZ:3736|3875|4263 ,TZ2:3875|4263 ,LM:4263 ,Weak_Layer:2670 ,Oceanic_crust:3000|3450|3567|3900|4112|4256 ,Oceanic_mantle:3280|3737|3938|4382 ,Sediments:2670 ,Upper_Crust:2800 ,Lower_Crust:3000|3450|3567|3900|4112|4256 ,Continental_mantle:3280|3737|3938|4382 ,Craton:3280|3737|3938|4382 ,plastic_strain:3300, Damping:4530
// //          */ 
// //          std::vector<double> T_ref = {293,293,293,293,293,293,293,293,293,293,293,293,350,650,800,890,1080,1315,1000,1000,1000,1000,326,510,803,803,803,803,803,803,1240,1240,1240,1240,1100,1100,1100,1100,293,293};
// //           
// //           
// //                 // Loop through all requested points
// //         // for (unsigned int c=0; c < out.densities.size(); ++c)
// //         for (unsigned int c=0; c < out.densities.size(); ++c)        
// //           {
// //               //for some reason it doesn't work if I just write >=12
// //             if(c==13 || c==14 || c==15 || c==16 || c==17 || c==18 || c==19 || c==20 || c==21 || c==22 || c==23 || c==24 || c==25 || c==26 || c==27 || c==28 || c==29 || c==30 ||  c==31 ||  c==32 ||  c==33 || c==34 || c==35 || c==36 || c==37 || c==38 || c==39 || c==40)
// //             {
// //               out.densities[c] = densities[c] * (1 - thermal_expansivities[c] * (in.temperature[input_index] - T_ref[c]));    
// //               out.thermal_expansion_coefficients[c] = thermal_expansivities[c];
// //               out.specific_heat_capacities[c] = specific_heats[c];
// //               out.compressibilities[c] = 0.0;
// //               out.entropy_derivative_pressure[c] = 0.0;
// //               out.entropy_derivative_temperature[c] = 0.0;
// //             }  
// //             else if(c==12)
// //             {
// //               out.densities[c] = densities[c];    
// //               out.thermal_expansion_coefficients[c] = thermal_expansivities[c];
// //               out.specific_heat_capacities[c] = specific_heats[c];
// //               out.compressibilities[c] = 0.0;
// //               out.entropy_derivative_pressure[c] = 0.0;
// //               out.entropy_derivative_temperature[c] = 0.0;
// //             }                     
// //             //for the rest the reference temperature is based on the adiabatic gradient
// //             else
// //              {
// //               out.densities[c] = densities[c] * (1 - thermal_expansivities[c] * (in.temperature[input_index] - (this->get_adiabatic_conditions().temperature(in.position[input_index]))));
// //               out.thermal_expansion_coefficients[c] = thermal_expansivities[c];
// //               out.specific_heat_capacities[c] = specific_heats[c];
// //               out.compressibilities[c] = 0.0;
// //               out.entropy_derivative_pressure[c] = 0.0;
// //               out.entropy_derivative_temperature[c] = 0.0;
// //              }
// //           } 
      }

      template <int dim>
      bool
      MulticomponentIncompressible<dim>::
      is_compressible () const
      {
        return false;
      }



      template <int dim>
      void
      MulticomponentIncompressible<dim>::declare_parameters (ParameterHandler &prm,
                                                             const double default_thermal_expansion)
      {
        prm.declare_entry ("Reference temperature", "293.",
                           Patterns::Double (0.),
                           "The reference temperature $T_0$. Units: \\si{\\kelvin}.");
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
        reference_T = prm.get_double ("Reference temperature");

        // Establish that a background field is required here
        const bool has_background_field = true;

        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // Parse multicomponent properties
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
