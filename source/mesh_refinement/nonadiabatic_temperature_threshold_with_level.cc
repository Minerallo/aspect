/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/nonadiabatic_temperature_threshold_with_level.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    NonadiabaticTemperatureThresholdwithLevel<dim>::tag_additional_cells () const
    {
      // tag_additional_cells is executed before the equations are solved
      // for the very first time. If we do not have the finite element, we
      // do not have the temperature and just do nothing in this plugin.
      if (this->get_dof_handler().n_locally_owned_dofs() == 0)
        return;

      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.temperature).get_unit_support_points());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values);

      std::vector<double> temperature_values (quadrature.size());
      const unsigned int n_dofs_per_cell = this->get_fe().base_element(this->introspection().base_elements.temperature).dofs_per_cell;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool refine = false;
            bool coarsen = false;
            bool clear_refine = false;
            bool clear_coarsen = false;
            bool valid_for_refinement =  false;

            fe_values.reinit(cell);
            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         temperature_values);

            // if the nonadiabatic temperature exceeds the threshold, cell is marked for refinement
            for (unsigned int j=0; j<n_dofs_per_cell; ++j)
              {
                const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(fe_values.quadrature_point(j));

                double nonadiabatic_temperature = 0;
                if (temperature_anomaly_type == absolute_value)
                  nonadiabatic_temperature = std::abs(temperature_values[j] - adiabatic_temperature);
                else if (temperature_anomaly_type == positive_only)
                  nonadiabatic_temperature = temperature_values[j] - adiabatic_temperature;
                else if (temperature_anomaly_type == negative_only)
                  nonadiabatic_temperature = adiabatic_temperature - temperature_values[j];
                else
                  AssertThrow (false, ExcNotImplemented());

                if (nonadiabatic_temperature > threshold)
                  {
                    valid_for_refinement = true;
                    break;
                  }
              }

          int max_refinement_level = max_level;

           if (valid_for_refinement){
                max_refinement_level = max_non_adiabatic_temperature_level;
            }            

            const int cell_level = cell->level();
            if (cell_level >= max_refinement_level)
                {
                clear_refine = true;
                }
            if (cell_level >  max_refinement_level)
                {
                coarsen = true;
                }
            if (cell_level <= max_refinement_level)
                {
                clear_coarsen = false;
                }
            if (cell_level < max_refinement_level)
                {
                refine = true;
                }    
                
            // if both coarsen and refine are true, give preference to refinement
            if (coarsen == true && refine == true)
            {
                coarsen = false;
                clear_refine = false;
            }

            // Perform the actual placement of the coarsening and refinement flags
            // We want to make sure that the refiment never goes below the minimum
            // or above the maximum, so we first check/set the coarsen/refine flag,
            // and then check/set the clear coarsen/refine flag.
            if (coarsen == true)
            {
                cell->set_coarsen_flag ();
            }
            if (clear_coarsen == true)
            {
                cell->clear_coarsen_flag ();
            }
            if (refine == true)
            {
                cell->set_refine_flag ();
            }
            if (clear_refine == true)
            {
                cell->clear_refine_flag ();
            }
          


            // if (refine)
            //   {
            //     cell->clear_coarsen_flag ();
            //     cell->set_refine_flag ();
            //   }
          }
    }

    template <int dim>
    void
    NonadiabaticTemperatureThresholdwithLevel<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Nonadiabatic temperature threshold");
        {
          prm.declare_entry ("Threshold",
                             "100",
                             Patterns::Double (0.),
                             "A threshold that the nonadiabatic temperature "
                             "will be evaluated against. "
                             "Units: \\si{\\kelvin}");
          prm.declare_entry ("Temperature anomaly type",
                             "absolute value",
                             Patterns::Selection ("negative only|positive only|absolute value"),
                             "What type of temperature anomaly should be considered when "
                             "evaluating against the threshold: Only negative anomalies "
                             "(negative only), only positive anomalies (positive only) "
                             "or the absolute value of the nonadiabatic temperature.");
          prm.declare_entry ("Nonadiabatic temperature maximum refinement level", "4", Patterns::Integer(0),
                             "The level of refinement of the mesh, cannot be more than the maximum level of refinement prescribed");                             
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    NonadiabaticTemperatureThresholdwithLevel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        min_level = prm.get_integer("Minimum refinement level");
        max_level = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");        
        prm.enter_subsection("Nonadiabatic temperature threshold with level");
        {
          threshold = prm.get_double("Threshold");
          max_non_adiabatic_temperature_level = prm.get_integer("Nonadiabatic temperature maximum refinement level");

          if (prm.get ("Temperature anomaly type") == "negative only")
            temperature_anomaly_type = negative_only;
          else if (prm.get ("Temperature anomaly type") == "positive only")
            temperature_anomaly_type = positive_only;
          else if (prm.get ("Temperature anomaly type") == "absolute value")
            temperature_anomaly_type = absolute_value;
          else
            AssertThrow (false, ExcNotImplemented());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(NonadiabaticTemperatureThresholdwithLevel,
                                              "nonadiabatic temperature threshold with level",
                                              "A mesh refinement criterion that computes refinement "
                                              "indicators from the temperature difference between the "
                                              "actual temperature and the adiabatic conditions (the "
                                              "nonadiabatic temperature). If the temperature anomaly "
                                              "exceeds the threshold given in the input file, the cell "
                                              "is marked for refinement.")
  }
}
