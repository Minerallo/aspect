/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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

#include <aspect/mesh_deformation/fastscape_cpp_helpers.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <fastscapelib/flow/flow_router.hpp>
#include <fastscapelib/flow/sink_resolver.hpp>

#include <deal.II/base/data_out_base.h>
#include <deal.II/numerics/data_out.h>

#if __has_include(<xtensor/generators/xbuilder.hpp>)
#  include <xtensor/generators/xbuilder.hpp>
#  include <xtensor/core/xmath.hpp>
#else
#  include <xtensor/xbuilder.hpp>
#  include <xtensor/xmath.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int surface_dim>
    void
    SpatialErosionStrength<surface_dim>::initialize(const std::string &filename,
                                                    const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      values = xt::ones<double>({surface_coordinates.size()});
      if (filename.empty())
        return;

      lookup = std::make_unique<Utilities::StructuredDataLookup<surface_dim>>(1, 1.0);
      lookup->load_file(filename, MPI_COMM_SELF);
      sample(surface_coordinates);
    }

    template <int surface_dim>
    void
    SpatialErosionStrength<surface_dim>::sample(const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      if (!lookup)
        return;

      values.resize({surface_coordinates.size()});
      for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
        values[i] = std::max(0.0,
                             lookup->get_data(surface_coordinates[i], 0));
    }

    template <int surface_dim>
    const xt::xarray<double> &
    SpatialErosionStrength<surface_dim>::get_values() const
    {
      return values;
    }

    template <int surface_dim>
    void
    SpatialSurfaceRunoff<surface_dim>::initialize(const std::string &filename,
                                                  const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      values = xt::ones<double>({surface_coordinates.size()});
      if (filename.empty())
        return;

      lookup = std::make_unique<Utilities::StructuredDataLookup<surface_dim>>(1, 1.0);
      lookup->load_file(filename, MPI_COMM_SELF);
      sample(surface_coordinates);
    }

    template <int surface_dim>
    void
    SpatialSurfaceRunoff<surface_dim>::sample(const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      if (!lookup)
        return;

      values.resize({surface_coordinates.size()});
      for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
        values[i] = std::max(0.0,
                             lookup->get_data(surface_coordinates[i], 0));
    }

    template <int surface_dim>
    const xt::xarray<double> &
    SpatialSurfaceRunoff<surface_dim>::get_values() const
    {
      return values;
    }

    template <int surface_dim>
    void
    SpatialIceThickness<surface_dim>::initialize(const std::string &filename,
                                                 const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      values = xt::zeros<double>({surface_coordinates.size()});
      if (filename.empty())
        return;

      lookup = std::make_unique<Utilities::StructuredDataLookup<surface_dim>>(1, 1.0);
      lookup->load_file(filename, MPI_COMM_SELF);
      sample(surface_coordinates);
    }

    template <int surface_dim>
    void
    SpatialIceThickness<surface_dim>::sample(const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      if (!lookup)
        return;

      values.resize({surface_coordinates.size()});
      for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
        values[i] = std::max(0.0,
                             lookup->get_data(surface_coordinates[i], 0));
    }

    template <int surface_dim>
    const xt::xarray<double> &
    SpatialIceThickness<surface_dim>::get_values() const
    {
      return values;
    }

    template <int surface_dim>
    void
    SpatialBasalIceVelocity<surface_dim>::initialize(const std::string &filename,
                                                     const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      values = xt::zeros<double>({surface_coordinates.size()});
      if (filename.empty())
        return;

      lookup = std::make_unique<Utilities::StructuredDataLookup<surface_dim>>(1, 1.0);
      lookup->load_file(filename, MPI_COMM_SELF);
      sample(surface_coordinates);
    }

    template <int surface_dim>
    void
    SpatialBasalIceVelocity<surface_dim>::sample(const std::vector<Point<surface_dim>> &surface_coordinates)
    {
      if (!lookup)
        return;

      values.resize({surface_coordinates.size()});
      for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
        values[i] = std::max(0.0,
                             lookup->get_data(surface_coordinates[i], 0));
    }

    template <int surface_dim>
    const xt::xarray<double> &
    SpatialBasalIceVelocity<surface_dim>::get_values() const
    {
      return values;
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::initialize(SurfaceMesh &surface_mesh,
                                                          const bool closed_surface,
                                                          const xt::xarray<double> &initial_elevation,
                                                          const double river_incision_coefficient,
                                                          const double area_exponent,
                                                          const double surface_slope_exponent,
                                                          const std::string &routing_method,
                                                          const std::string &routing_connectivity,
                                                          const double multi_flow_slope_exponent,
                                                          const double solver_tolerance,
                                                          const double marine_transport_coefficient,
                                                          const double sediment_porosity,
                                                          const double transport_depth_scale,
                                                          const bool limit_deposition_to_available_accommodation,
                                                          const double maximum_deposition_above_sea_level,
                                                          const bool use_sea_level_base_level,
                                                          const bool restrict_ocean_connectivity,
                                                          const double submarine_incision_factor,
                                                          const std::string &marine_open_boundary,
                                                          const std::vector<std::string> &rock_names,
                                                          const std::vector<double> &rock_probabilities,
                                                          const std::vector<double> &rock_erodibility_factors,
                                                          const unsigned int rock_random_seed)
    {
      spherical_geometry = closed_surface;
      drainage_area_exponent = area_exponent;
      marine_sediment_transport_coefficient = marine_transport_coefficient;
      marine_sediment_porosity = sediment_porosity;
      marine_transport_depth_scale = transport_depth_scale;
      limit_marine_deposition_to_available_accommodation =
        limit_deposition_to_available_accommodation;
      maximum_marine_deposition_above_sea_level =
        maximum_deposition_above_sea_level;
      use_sea_level_as_drainage_base_level = use_sea_level_base_level;
      restrict_ocean_to_largest_connected_component =
        restrict_ocean_connectivity;
      submarine_river_incision_factor = submarine_incision_factor;
      open_marine_sediment_boundary = marine_open_boundary;
      grid = std::make_unique<Grid>(surface_mesh,
                                    closed_surface,
                                    routing_connectivity ==
                                    "face and vertex");
      if (routing_method == "multiple flow")
        flow_graph = std::make_unique<FlowGraph>(
                       *grid,
                       typename FlowGraph::operators_type
        {
          std::make_shared<fastscapelib::single_flow_router>(),
          std::make_shared<fastscapelib::mst_sink_resolver>(),
          std::make_shared<fastscapelib::multi_flow_router>(
            multi_flow_slope_exponent)
        });
      else
        flow_graph = std::make_unique<FlowGraph>(
                       *grid,
                       typename FlowGraph::operators_type
        {
          std::make_shared<fastscapelib::single_flow_router>(),
          std::make_shared<fastscapelib::mst_sink_resolver>()
        });
      eroder = std::make_unique<Eroder>(*flow_graph,
                                        river_incision_coefficient,
                                        area_exponent,
                                        surface_slope_exponent,
                                        solver_tolerance);

      elevation = initial_elevation;
      bedrock_elevation = initial_elevation;
      drainage_area = xt::zeros<double>(flow_graph->grid_shape());
      geometric_drainage_area = xt::zeros<double>(flow_graph->grid_shape());
      unit_surface_runoff = xt::ones<double>(flow_graph->grid_shape());
      dominant_drainage_basin = xt::zeros<double>(flow_graph->grid_shape());
      dominant_outlet_node = xt::zeros<double>(flow_graph->grid_shape());
      primary_receiver_node = xt::zeros<double>(flow_graph->grid_shape());
      primary_receiver_fraction = xt::zeros<double>(flow_graph->grid_shape());
      flow_receiver_count = xt::zeros<double>(flow_graph->grid_shape());
      drainage_outlet_mask = xt::zeros<double>(flow_graph->grid_shape());
      coastal_outlet_mask = xt::zeros<double>(flow_graph->grid_shape());
      erosion = xt::zeros<double>(flow_graph->grid_shape());
      fluvial_erosion = xt::zeros<double>(flow_graph->grid_shape());
      glacial_erosion = xt::zeros<double>(flow_graph->grid_shape());
      accumulated_fluvial_erosion =
        xt::zeros<double>(flow_graph->grid_shape());
      accumulated_glacial_erosion =
        xt::zeros<double>(flow_graph->grid_shape());
      modeled_ice_thickness = xt::zeros<double>(flow_graph->grid_shape());
      modeled_basal_ice_velocity = xt::zeros<double>(flow_graph->grid_shape());
      routed_ice_discharge = xt::zeros<double>(flow_graph->grid_shape());
      routed_glacier_width = xt::zeros<double>(flow_graph->grid_shape());
      routed_ice_mass_balance = xt::zeros<double>(flow_graph->grid_shape());
      sediment_flux = xt::zeros<double>(flow_graph->grid_shape());
      marine_sediment_flux = xt::zeros<double>(flow_graph->grid_shape());
      sediment_thickness = xt::zeros<double>(flow_graph->grid_shape());
      deposition_rate = xt::zeros<double>(flow_graph->grid_shape());
      ocean_mask = xt::zeros<double>(flow_graph->grid_shape());

      lithology_names = rock_names;
      lithology_erodibility_factors = rock_erodibility_factors;
      bedrock_lithology.resize(initial_elevation.size());
      sediment_flux_by_lithology.clear();
      sediment_thickness_by_lithology.clear();
      sediment_thickness_at_last_output_by_lithology.clear();
      for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
        {
          sediment_flux_by_lithology.push_back(
            xt::zeros<double>(flow_graph->grid_shape()));
          sediment_thickness_by_lithology.push_back(
            xt::zeros<double>(flow_graph->grid_shape()));
          sediment_thickness_at_last_output_by_lithology.push_back(
            xt::zeros<double>(flow_graph->grid_shape()));
        }

      std::vector<double> cumulative_probability(rock_probabilities.size());
      std::partial_sum(rock_probabilities.begin(), rock_probabilities.end(),
                       cumulative_probability.begin());
      const double probability_sum = cumulative_probability.back();
      for (unsigned int i = 0; i < bedrock_lithology.size(); ++i)
        {
          // SplitMix64 gives a reproducible cell-wise realization without
          // depending on a standard-library random-number implementation.
          std::uint64_t value =
            static_cast<std::uint64_t>(i) +
            (static_cast<std::uint64_t>(rock_random_seed) << 32);
          value += 0x9e3779b97f4a7c15ULL;
          value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
          value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
          value ^= value >> 31;
          const double sample =
            probability_sum * static_cast<double>(value >> 11) /
            static_cast<double>(std::uint64_t(1) << 53);
          bedrock_lithology[i] = static_cast<unsigned int>(
                                   std::lower_bound(cumulative_probability.begin(),
                                                    cumulative_probability.end(), sample) -
                                   cumulative_probability.begin());
          bedrock_lithology[i] =
            std::min<unsigned int>(bedrock_lithology[i],
                                   lithology_names.size()-1);
        }
    }

    template <int surface_dim, int space_dim>
    typename FastscapeLandscape<surface_dim,space_dim>::StepResult
    FastscapeLandscape<surface_dim,space_dim>::advance(const xt::xarray<double> &uplift_rate,
                                                       const std::vector<SurfaceVelocity> &tangential_velocity,
                                                       const double total_time_years,
                                                       unsigned int number_of_steps,
                                                       const double maximum_step_years,
                                                       const double sea_level,
                                                       const xt::xarray<double> &erosion_strength,
                                                       const xt::xarray<double> &surface_runoff,
                                                       const xt::xarray<double> &ice_thickness,
                                                       const xt::xarray<double> &basal_ice_velocity,
                                                       const std::string &glacial_erosion_mode,
                                                       const double routed_glacier_ela,
                                                       const double routed_glacier_accumulation_gradient,
                                                       const double routed_glacier_maximum_accumulation,
                                                       const double routed_glacier_ablation_gradient,
                                                       const double routed_glacier_maximum_ablation,
                                                       const double routed_glacier_minimum_discharge,
                                                       const double routed_glacier_reference_discharge,
                                                       const double routed_glacier_reference_width,
                                                       const double routed_glacier_width_exponent,
                                                       const double routed_glacier_minimum_width,
                                                       const double routed_glacier_maximum_width,
                                                       const double routed_glacier_reference_thickness,
                                                       const double routed_glacier_thickness_exponent,
                                                       const double routed_glacier_minimum_thickness,
                                                       const double routed_glacier_maximum_thickness,
                                                       const double routed_glacier_thickness_slope_exponent,
                                                       const bool routed_glacier_spread_erosion,
                                                       const bool routed_glacier_terminate_at_sea_level,
                                                       const double glacial_erosion_coefficient,
                                                       const double glacial_velocity_exponent,
                                                       const double glacial_ice_thickness_scale,
                                                       const double minimum_ice_thickness,
                                                       const bool restrict_glacial_erosion_to_grounded_ice,
                                                       const double minimum_glacial_erosion_elevation,
                                                       const double ice_density,
                                                       const double seawater_density,
                                                       const bool advect_surface_state,
                                                       const double maximum_advection_courant,
                                                       const double hillslope_diffusivity,
                                                       const double submarine_hillslope_diffusivity,
                                                       const double maximum_marine_transport_courant,
                                                       const double maximum_diffusion_courant)
    {
      Assert(grid && flow_graph && eroder, ExcInternalError());
      AssertDimension(uplift_rate.size(), elevation.size());
      AssertDimension(tangential_velocity.size(), elevation.size());
      AssertDimension(erosion_strength.size(), elevation.size());
      AssertDimension(surface_runoff.size(), elevation.size());
      AssertDimension(ice_thickness.size(), elevation.size());
      AssertDimension(basal_ice_velocity.size(), elevation.size());

      number_of_steps = std::max(1u, number_of_steps);
      while (total_time_years / number_of_steps > maximum_step_years)
        number_of_steps *= 2;
      const double step_years = total_time_years / number_of_steps;

      StepResult result;
      result.previous_elevation = elevation;
      accumulated_fluvial_erosion.fill(0.0);
      accumulated_glacial_erosion.fill(0.0);

      for (unsigned int step = 0; step < number_of_steps; ++step)
        {
          if (advect_surface_state)
            advect_surface_fields(tangential_velocity,
                                  step_years,
                                  maximum_advection_courant);
          bedrock_elevation += step_years * uplift_rate;
          const xt::xarray<double> uplifted =
            bedrock_elevation + sediment_thickness;
          set_base_levels(uplifted, sea_level);
          flow_graph->update_routes(uplifted);
          flow_graph->accumulate(drainage_area, surface_runoff);
          flow_graph->accumulate(geometric_drainage_area,
                                 unit_surface_runoff);
          update_drainage_diagnostics();

          if (glacial_erosion_mode == "routed glacier")
            update_routed_glacier_fields(
              uplifted,
              routed_glacier_ela,
              routed_glacier_accumulation_gradient,
              routed_glacier_maximum_accumulation,
              routed_glacier_ablation_gradient,
              routed_glacier_maximum_ablation,
              routed_glacier_minimum_discharge,
              routed_glacier_reference_discharge,
              routed_glacier_reference_width,
              routed_glacier_width_exponent,
              routed_glacier_minimum_width,
              routed_glacier_maximum_width,
              routed_glacier_reference_thickness,
              routed_glacier_thickness_exponent,
              routed_glacier_minimum_thickness,
              routed_glacier_maximum_thickness,
              routed_glacier_thickness_slope_exponent,
              sea_level,
              routed_glacier_terminate_at_sea_level);
          else
            {
              modeled_ice_thickness = ice_thickness;
              modeled_basal_ice_velocity = basal_ice_velocity;
              routed_ice_discharge.fill(0.0);
              routed_glacier_width.fill(0.0);
              routed_ice_mass_balance.fill(0.0);
            }

          // Erosion strength is a local multiplier. Surface runoff is
          // handled separately above so that water supplied upstream is
          // carried through the drainage network.
          xt::xarray<double> effective_drainage_area = drainage_area;
          xt::xarray<double> exposed_erodibility =
            xt::ones<double>(flow_graph->grid_shape());
          for (unsigned int i = 0; i < effective_drainage_area.size(); ++i)
            {
              exposed_erodibility[i] =
                lithology_erodibility_factors[bedrock_lithology[i]];
              if (sediment_thickness[i] > 0.0)
                {
                  exposed_erodibility[i] = 0.0;
                  for (unsigned int rock = 0;
                       rock < lithology_names.size(); ++rock)
                    exposed_erodibility[i] +=
                      lithology_erodibility_factors[rock] *
                      sediment_thickness_by_lithology[rock][i] /
                      sediment_thickness[i];
                }
              effective_drainage_area[i] *=
                std::pow(erosion_strength[i] * exposed_erodibility[i],
                         1.0 / drainage_area_exponent);
            }

          fluvial_erosion =
            eroder->erode(uplifted, effective_drainage_area, step_years);
          for (unsigned int i = 0; i < fluvial_erosion.size(); ++i)
            if (ocean_mask[i] >= 0.5)
              fluvial_erosion[i] *= submarine_river_incision_factor;
          for (unsigned int i = 0; i < glacial_erosion.size(); ++i)
            {
              const bool above_minimum_elevation =
                uplifted[i] >= minimum_glacial_erosion_elevation;
              const bool grounded =
                uplifted[i] >= sea_level ||
                ice_density * modeled_ice_thickness[i] >=
                seawater_density * (sea_level - uplifted[i]);
              glacial_erosion[i] =
                modeled_ice_thickness[i] >= minimum_ice_thickness &&
                above_minimum_elevation &&
                (!restrict_glacial_erosion_to_grounded_ice || grounded)
                ? step_years * glacial_erosion_coefficient *
                exposed_erodibility[i] *
                std::pow(modeled_basal_ice_velocity[i],
                         glacial_velocity_exponent) *
                (glacial_ice_thickness_scale > 0.0
                 ? 1.0 - std::exp(-modeled_ice_thickness[i] /
                                  glacial_ice_thickness_scale)
                 : 1.0)
                : 0.0;
            }
          if (glacial_erosion_mode == "routed glacier" &&
              routed_glacier_spread_erosion)
            glacial_erosion = spread_routed_glacial_erosion(
                                glacial_erosion,
                                routed_glacier_width);
          erosion = fluvial_erosion + glacial_erosion;
          accumulated_fluvial_erosion += fluvial_erosion;
          accumulated_glacial_erosion += glacial_erosion;

          const auto areas = grid->nodes_areas();
          std::vector<xt::xarray<double>> local_source_by_lithology;
          local_source_by_lithology.reserve(lithology_names.size());
          for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
            local_source_by_lithology.push_back(
              xt::zeros<double>(flow_graph->grid_shape()));

          for (unsigned int i = 0; i < erosion.size(); ++i)
            {
              const double sediment_erosion =
                std::min(sediment_thickness[i], erosion[i]);
              if (sediment_erosion > 0.0)
                for (unsigned int rock = 0;
                     rock < lithology_names.size(); ++rock)
                  {
                    const double removed = sediment_erosion *
                                           sediment_thickness_by_lithology[rock][i] /
                                           sediment_thickness[i];
                    sediment_thickness_by_lithology[rock][i] -= removed;
                    // Mobile sediment thickness includes pore space, whereas
                    // routed sediment flux is solid volume. Without this
                    // conversion, every erosion/deposition cycle expands
                    // recycled sediment by 1/(1-porosity).
                    local_source_by_lithology[rock][i] +=
                      removed * (1.0 - marine_sediment_porosity);
                  }
              const double bedrock_erosion = erosion[i] - sediment_erosion;
              local_source_by_lithology[bedrock_lithology[i]][i] +=
                bedrock_erosion;
              sediment_thickness[i] -= sediment_erosion;
              bedrock_elevation[i] -= bedrock_erosion;
              result.eroded_volume += erosion[i] * areas[i];
              result.fluvial_eroded_volume += fluvial_erosion[i] * areas[i];
              result.glacial_eroded_volume += glacial_erosion[i] * areas[i];
            }

          sediment_flux.fill(0.0);
          for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
            {
              sediment_flux_by_lithology[rock] =
                flow_graph->accumulate(
                  local_source_by_lithology[rock] / step_years);
              sediment_flux += sediment_flux_by_lithology[rock];
            }

          elevation = bedrock_elevation + sediment_thickness;
          if (hillslope_diffusivity > 0.0)
            diffuse_hillslopes(step_years,
                               hillslope_diffusivity,
                               submarine_hillslope_diffusivity,
                               maximum_diffusion_courant);
          if (marine_sediment_transport_coefficient > 0.0)
            {
              const xt::xarray<double> sediment_before_deposition =
                sediment_thickness;
              for (const auto index : flow_graph->base_levels())
                if (elevation[index] <= sea_level)
                  {
                    const double solid_volume = sediment_flux[index] * step_years;
                    double deposited_solid_volume = solid_volume;
                    if (limit_marine_deposition_to_available_accommodation)
                      {
                        const double maximum_surface_elevation =
                          sea_level + maximum_marine_deposition_above_sea_level;
                        const double available_bulk_volume =
                          std::max(0.0,
                                   maximum_surface_elevation - elevation[index]) *
                          areas[index];
                        const double available_solid_volume =
                          available_bulk_volume *
                          (1.0 - marine_sediment_porosity);
                        deposited_solid_volume =
                          std::min(solid_volume, available_solid_volume);
                      }
                    const double deposited_fraction = solid_volume > 0.0
                                                      ? deposited_solid_volume / solid_volume
                                                      : 0.0;
                    for (unsigned int rock = 0;
                         rock < lithology_names.size(); ++rock)
                      {
                        const double deposited_thickness =
                          sediment_flux_by_lithology[rock][index] *
                          step_years * deposited_fraction /
                          ((1.0 - marine_sediment_porosity) * areas[index]);
                        sediment_thickness_by_lithology[rock][index] +=
                          deposited_thickness;
                        sediment_thickness[index] += deposited_thickness;
                      }
                    const double accommodation_limited_export =
                      solid_volume - deposited_solid_volume;
                    result.exported_sediment_flux +=
                      accommodation_limited_export;
                    result.accommodation_limited_exported_sediment_flux +=
                      accommodation_limited_export;
                    result.coastal_sediment_flux += solid_volume;
                    result.deposited_sediment_volume +=
                      deposited_solid_volume /
                      (1.0 - marine_sediment_porosity);
                  }

              elevation = bedrock_elevation + sediment_thickness;
              transport_marine_sediment(step_years,
                                        sea_level,
                                        maximum_marine_transport_courant);
              export_marine_sediment(result);
              elevation = bedrock_elevation + sediment_thickness;
              for (unsigned int i = 0; i < sediment_thickness.size(); ++i)
                deposition_rate[i] =
                  (sediment_thickness[i] -
                   sediment_before_deposition[i]) / step_years;
            }
          else
            {
              deposition_rate.fill(0.0);
              marine_sediment_flux.fill(0.0);
            }
        }

      if (marine_sediment_transport_coefficient > 0.0 &&
          total_time_years > 0.0)
        {
          result.coastal_sediment_flux /= total_time_years;
          result.exported_sediment_flux /= total_time_years;
          result.accommodation_limited_exported_sediment_flux /=
            total_time_years;
        }

      if (marine_sediment_transport_coefficient == 0.0)
        for (const auto index : flow_graph->base_levels())
          result.exported_sediment_flux += sediment_flux[index];

      const auto areas = grid->nodes_areas();
      for (unsigned int i = 0; i < sediment_thickness.size(); ++i)
        result.stored_sediment_volume +=
          sediment_thickness[i] * areas[i] *
          (1.0 - marine_sediment_porosity);

      return result;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_elevation() const
    {
      return elevation;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_drainage_area() const
    {
      return drainage_area;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_geometric_drainage_area() const
    {
      return geometric_drainage_area;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_dominant_drainage_basin() const
    {
      return dominant_drainage_basin;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_dominant_outlet_node() const
    {
      return dominant_outlet_node;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_primary_receiver_node() const
    {
      return primary_receiver_node;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_primary_receiver_fraction() const
    {
      return primary_receiver_fraction;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_flow_receiver_count() const
    {
      return flow_receiver_count;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_drainage_outlet_mask() const
    {
      return drainage_outlet_mask;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_coastal_outlet_mask() const
    {
      return coastal_outlet_mask;
    }

    template <int surface_dim, int space_dim>
    xt::xarray<double> FastscapeLandscape<surface_dim,space_dim>::get_cell_areas() const
    {
      return grid->nodes_areas();
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_erosion() const
    {
      return erosion;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_fluvial_erosion() const
    {
      return fluvial_erosion;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_glacial_erosion() const
    {
      return glacial_erosion;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_accumulated_fluvial_erosion() const
    {
      return accumulated_fluvial_erosion;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_accumulated_glacial_erosion() const
    {
      return accumulated_glacial_erosion;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_modeled_ice_thickness() const
    {
      return modeled_ice_thickness;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_modeled_basal_ice_velocity() const
    {
      return modeled_basal_ice_velocity;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_routed_ice_discharge() const
    {
      return routed_ice_discharge;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_routed_glacier_width() const
    {
      return routed_glacier_width;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_routed_ice_mass_balance() const
    {
      return routed_ice_mass_balance;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_sediment_flux() const
    {
      return sediment_flux;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_marine_sediment_flux() const
    {
      return marine_sediment_flux;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_sediment_thickness() const
    {
      return sediment_thickness;
    }

    template <int surface_dim, int space_dim>
    const std::vector<std::string> &FastscapeLandscape<surface_dim,space_dim>::get_lithology_names() const
    {
      return lithology_names;
    }

    template <int surface_dim, int space_dim>
    const std::vector<unsigned int> &FastscapeLandscape<surface_dim,space_dim>::get_bedrock_lithology() const
    {
      return bedrock_lithology;
    }

    template <int surface_dim, int space_dim>
    const std::vector<xt::xarray<double>> &
    FastscapeLandscape<surface_dim,space_dim>::get_sediment_flux_by_lithology() const
    {
      return sediment_flux_by_lithology;
    }

    template <int surface_dim, int space_dim>
    const std::vector<xt::xarray<double>> &
    FastscapeLandscape<surface_dim,space_dim>::get_sediment_thickness_by_lithology() const
    {
      return sediment_thickness_by_lithology;
    }

    template <int surface_dim, int space_dim>
    std::vector<xt::xarray<double>>
    FastscapeLandscape<surface_dim,space_dim>::take_deposited_thickness_by_lithology()
    {
      std::vector<xt::xarray<double>> result;
      result.reserve(lithology_names.size());
      for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
        {
          result.push_back(xt::maximum(
                             sediment_thickness_by_lithology[rock] -
                             sediment_thickness_at_last_output_by_lithology[rock], 0.0));
          sediment_thickness_at_last_output_by_lithology[rock] =
            sediment_thickness_by_lithology[rock];
        }
      return result;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_deposition_rate() const
    {
      return deposition_rate;
    }

    template <int surface_dim, int space_dim>
    const xt::xarray<double> &FastscapeLandscape<surface_dim,space_dim>::get_ocean_mask() const
    {
      return ocean_mask;
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::set_elevation(const std::vector<double> &values)
    {
      AssertDimension(values.size(), elevation.size());
      std::copy(values.begin(), values.end(), elevation.begin());
      for (unsigned int i = 0; i < elevation.size(); ++i)
        bedrock_elevation[i] = elevation[i] - sediment_thickness[i];
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::set_sediment_thickness(const std::vector<double> &values)
    {
      AssertDimension(values.size(), sediment_thickness.size());
      std::copy(values.begin(), values.end(), sediment_thickness.begin());
      for (auto &field : sediment_thickness_by_lithology)
        field.fill(0.0);
      if (!sediment_thickness_by_lithology.empty())
        {
          std::copy(values.begin(), values.end(),
                    sediment_thickness_by_lithology[0].begin());
          sediment_thickness_at_last_output_by_lithology[0] =
            sediment_thickness_by_lithology[0];
        }
      for (unsigned int i = 0; i < elevation.size(); ++i)
        bedrock_elevation[i] = elevation[i] - sediment_thickness[i];
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::set_lithology_state(
      const std::vector<unsigned int> &stored_bedrock_lithology,
      const std::vector<std::vector<double>> &stored_sediment_thickness)
    {
      AssertDimension(stored_bedrock_lithology.size(),
                      bedrock_lithology.size());
      AssertDimension(stored_sediment_thickness.size(),
                      sediment_thickness_by_lithology.size());
      bedrock_lithology = stored_bedrock_lithology;
      sediment_thickness.fill(0.0);
      for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
        {
          AssertDimension(stored_sediment_thickness[rock].size(),
                          sediment_thickness.size());
          std::copy(stored_sediment_thickness[rock].begin(),
                    stored_sediment_thickness[rock].end(),
                    sediment_thickness_by_lithology[rock].begin());
          sediment_thickness += sediment_thickness_by_lithology[rock];
          sediment_thickness_at_last_output_by_lithology[rock] =
            sediment_thickness_by_lithology[rock];
        }
      for (unsigned int i = 0; i < elevation.size(); ++i)
        bedrock_elevation[i] = elevation[i] - sediment_thickness[i];
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::update_routed_glacier_fields(
      const xt::xarray<double> &surface_elevation,
      const double ela,
      const double accumulation_gradient,
      const double maximum_accumulation,
      const double ablation_gradient,
      const double maximum_ablation,
      const double minimum_discharge,
      const double reference_discharge,
      const double reference_width,
      const double width_exponent,
      const double minimum_width,
      const double maximum_width,
      const double reference_thickness,
      const double thickness_exponent,
      const double minimum_thickness,
      const double maximum_thickness,
      const double thickness_slope_exponent,
      const double sea_level,
      const bool terminate_at_sea_level)
    {
      const auto areas = grid->nodes_areas();
      routed_ice_discharge.fill(0.0);
      routed_glacier_width.fill(0.0);
      modeled_ice_thickness.fill(0.0);
      modeled_basal_ice_velocity.fill(0.0);

      for (unsigned int i = 0; i < surface_elevation.size(); ++i)
        routed_ice_mass_balance[i] =
          surface_elevation[i] >= ela
          ? std::min(maximum_accumulation,
                     accumulation_gradient *
                     (surface_elevation[i] - ela))
          : -std::min(maximum_ablation,
                      ablation_gradient *
                      (ela - surface_elevation[i]));

      const auto &graph = flow_graph->impl();
      const auto &receivers = graph.receivers();
      const auto &receiver_count = graph.receivers_count();
      const auto &receiver_weight = graph.receivers_weight();
      const auto nodes = graph.nodes_indices_bottomup();

      for (auto node = nodes.rbegin(); node != nodes.rend(); ++node)
        {
          const std::size_t i = *node;
          const double local_balance_volume =
            routed_ice_mass_balance[i] * areas[i];
          routed_ice_discharge[i] =
            std::max(0.0,
                     routed_ice_discharge[i] + local_balance_volume);

          if (terminate_at_sea_level && surface_elevation[i] <= sea_level)
            {
              routed_ice_discharge[i] = 0.0;
              continue;
            }

          const double discharge = routed_ice_discharge[i];
          if (discharge >= minimum_discharge)
            {
              const double discharge_ratio =
                discharge / reference_discharge;
              routed_glacier_width[i] =
                std::clamp(reference_width *
                           std::pow(discharge_ratio, width_exponent),
                           minimum_width,
                           maximum_width);

              double steepest_slope = 0.0;
              for (std::size_t n = 0;
                   n < grid->number_of_cell_neighbors(i); ++n)
                {
                  const std::size_t neighbor = grid->cell_neighbor(i, n);
                  const double distance = grid->cell_neighbor_distance(i, n);
                  if (distance > 0.0)
                    steepest_slope =
                      std::max(steepest_slope,
                               (surface_elevation[i] -
                                surface_elevation[neighbor]) / distance);
                }
              const double slope_factor =
                std::pow(0.05 / std::max(steepest_slope, 1e-3),
                         thickness_slope_exponent);
              modeled_ice_thickness[i] =
                std::clamp(reference_thickness *
                           std::pow(discharge_ratio,
                                    thickness_exponent) *
                           slope_factor,
                           minimum_thickness,
                           maximum_thickness);
              modeled_basal_ice_velocity[i] =
                discharge /
                (routed_glacier_width[i] *
                 modeled_ice_thickness[i]);
            }

          for (std::size_t receiver_number = 0;
               receiver_number < receiver_count[i]; ++receiver_number)
            {
              const std::size_t receiver =
                receivers(i, receiver_number);
              if (receiver != i)
                routed_ice_discharge[receiver] +=
                  discharge * receiver_weight(i, receiver_number);
            }
        }
    }

    template <int surface_dim, int space_dim>
    xt::xarray<double>
    FastscapeLandscape<surface_dim,space_dim>::spread_routed_glacial_erosion(
      const xt::xarray<double> &centerline_erosion,
      const xt::xarray<double> &glacier_width) const
    {
      const auto areas = grid->nodes_areas();
      xt::xarray<double> spread =
        xt::zeros<double>(flow_graph->grid_shape());
      const double infinity = std::numeric_limits<double>::infinity();
      std::vector<double> distance(centerline_erosion.size(), infinity);
      std::vector<std::size_t> touched;
      std::vector<std::pair<std::size_t,double>> recipients;
      using QueueEntry = std::pair<double,std::size_t>;

      for (std::size_t source = 0;
           source < centerline_erosion.size(); ++source)
        {
          if (centerline_erosion[source] <= 0.0)
            continue;

          const double radius = 0.5 * glacier_width[source];
          if (radius <= 0.0)
            {
              spread[source] += centerline_erosion[source];
              continue;
            }

          const double sigma =
            std::max(0.5 * radius,
                     0.25 * std::sqrt(areas[source]));
          std::priority_queue<QueueEntry,
              std::vector<QueueEntry>,
              std::greater<QueueEntry>> frontier;
          touched.clear();
          recipients.clear();
          distance[source] = 0.0;
          touched.push_back(source);
          frontier.push({0.0, source});

          while (!frontier.empty())
            {
              const auto [current_distance, i] = frontier.top();
              frontier.pop();
              if (current_distance != distance[i])
                continue;
              if (current_distance > radius)
                continue;

              const double weight =
                std::exp(-0.5 * current_distance * current_distance /
                         (sigma * sigma));
              recipients.push_back({i, weight});

              for (std::size_t n = 0;
                   n < grid->number_of_cell_neighbors(i); ++n)
                {
                  const std::size_t neighbor = grid->cell_neighbor(i, n);
                  const double candidate =
                    current_distance +
                    grid->cell_neighbor_distance(i, n);
                  if (candidate <= radius &&
                      candidate < distance[neighbor])
                    {
                      if (!std::isfinite(distance[neighbor]))
                        touched.push_back(neighbor);
                      distance[neighbor] = candidate;
                      frontier.push({candidate, neighbor});
                    }
                }
            }

          double normalization = 0.0;
          for (const auto &[i, weight] : recipients)
            normalization += weight * areas[i];
          const double eroded_volume =
            centerline_erosion[source] * areas[source];
          for (const auto &[i, weight] : recipients)
            spread[i] += eroded_volume * weight / normalization;

          for (const std::size_t i : touched)
            distance[i] = infinity;
        }
      return spread;
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::diffuse_hillslopes(const double step_years,
                                                                  const double diffusivity,
                                                                  const double submarine_diffusivity,
                                                                  const double maximum_courant)
    {
      AssertThrow(maximum_courant > 0.0,
                  ExcMessage("Maximum hillslope-diffusion Courant number "
                             "must be positive."));
      const auto areas = grid->nodes_areas();
      std::vector<double> conductance_sum(elevation.size(), 0.0);
      for (std::size_t i = 0; i < elevation.size(); ++i)
        for (std::size_t n = 0;
             n < grid->number_of_cell_neighbors(i); ++n)
          {
            const std::size_t j = grid->cell_neighbor(i, n);
            if (j <= i)
              continue;
            const double local_diffusivity =
              ocean_mask[i] >= 0.5 && ocean_mask[j] >= 0.5
              ? submarine_diffusivity
              : diffusivity;
            const double conductance =
              local_diffusivity *
              grid->cell_shared_face_measure(i, n) /
              grid->cell_neighbor_distance(i, n);
            conductance_sum[i] += conductance;
            conductance_sum[j] += conductance;
          }

      double largest_courant = 0.0;
      for (std::size_t i = 0; i < elevation.size(); ++i)
        largest_courant =
          std::max(largest_courant,
                   step_years * conductance_sum[i] / areas[i]);
      const unsigned int diffusion_steps =
        std::max(1u,
                 static_cast<unsigned int>(
                   std::ceil(largest_courant / maximum_courant)));
      const double diffusion_step_years =
        step_years / diffusion_steps;

      for (unsigned int step = 0; step < diffusion_steps; ++step)
        {
          std::vector<double> volume_change(elevation.size(), 0.0);
          for (std::size_t i = 0; i < elevation.size(); ++i)
            for (std::size_t n = 0;
                 n < grid->number_of_cell_neighbors(i); ++n)
              {
                const std::size_t j = grid->cell_neighbor(i, n);
                if (j <= i)
                  continue;
                const double local_diffusivity =
                  ocean_mask[i] >= 0.5 && ocean_mask[j] >= 0.5
                  ? submarine_diffusivity
                  : diffusivity;
                const double volume =
                  local_diffusivity *
                  grid->cell_shared_face_measure(i, n) /
                  grid->cell_neighbor_distance(i, n) *
                  (elevation[i] - elevation[j]) *
                  diffusion_step_years;
                volume_change[i] -= volume;
                volume_change[j] += volume;
              }

          for (std::size_t i = 0; i < elevation.size(); ++i)
            {
              const double previous_sediment_thickness = sediment_thickness[i];
              const double elevation_change =
                volume_change[i] / areas[i];
              if (elevation_change >= 0.0)
                {
                  sediment_thickness[i] += elevation_change;
                  // Hillslope diffusion does not retain an explicit donor
                  // list after face transfers are assembled. Attribute the
                  // local deposit to the exposed rock as a conservative
                  // fallback; river and marine routing preserve the full
                  // transported mixture.
                  sediment_thickness_by_lithology[bedrock_lithology[i]][i] +=
                    elevation_change;
                }
              else
                {
                  const double removal = -elevation_change;
                  const double sediment_removal =
                    std::min(sediment_thickness[i], removal);
                  if (sediment_removal > 0.0)
                    for (unsigned int rock = 0;
                         rock < lithology_names.size(); ++rock)
                      sediment_thickness_by_lithology[rock][i] *=
                        (previous_sediment_thickness - sediment_removal) /
                        previous_sediment_thickness;
                  sediment_thickness[i] -= sediment_removal;
                  bedrock_elevation[i] -= removal - sediment_removal;
                }
              elevation[i] =
                bedrock_elevation[i] + sediment_thickness[i];
            }
        }
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::export_marine_sediment(StepResult &result)
    {
      if (open_marine_sediment_boundary == "none")
        return;

      for (std::size_t i = 0; i < sediment_thickness.size(); ++i)
        {
          const bool on_open_boundary =
            (open_marine_sediment_boundary == "minimum x" &&
             grid->cell_touches_boundary(i, 0, false)) ||
            (open_marine_sediment_boundary == "maximum x" &&
             grid->cell_touches_boundary(i, 0, true)) ||
            (open_marine_sediment_boundary == "all nonperiodic" &&
             [&]()
          {
            for (unsigned int d = 0; d < surface_dim; ++d)
              if (grid->cell_touches_boundary(i, d, false) ||
                  grid->cell_touches_boundary(i, d, true))
                return true;
            return false;
          }());
          if (!on_open_boundary || ocean_mask[i] < 0.5 ||
              sediment_thickness[i] <= 0.0)
            continue;

          const double bulk_volume =
            sediment_thickness[i] * grid->nodes_areas(i);
          result.exported_sediment_flux +=
            bulk_volume * (1.0 - marine_sediment_porosity);
          sediment_thickness[i] = 0.0;
          for (auto &field : sediment_thickness_by_lithology)
            field[i] = 0.0;
        }
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::advect_surface_fields(const std::vector<SurfaceVelocity> &velocity,
                                                                     const double step_years,
                                                                     const double maximum_courant)
    {
      AssertThrow(maximum_courant > 0.0,
                  ExcMessage("Maximum surface-advection Courant number "
                             "must be positive."));
      const auto areas = grid->nodes_areas();
      std::vector<double> absolute_face_rate(elevation.size(), 0.0);

      for (std::size_t i = 0; i < elevation.size(); ++i)
        for (std::size_t n = 0;
             n < grid->number_of_cell_neighbors(i); ++n)
          {
            const std::size_t j = grid->cell_neighbor(i, n);
            if (j <= i)
              continue;
            const SurfaceVelocity face_velocity =
              0.5 * (velocity[i] + velocity[j]);
            const double rate =
              (face_velocity * grid->cell_neighbor_direction(i, n)) *
              grid->cell_shared_face_measure(i, n);
            absolute_face_rate[i] += std::abs(rate);
            absolute_face_rate[j] += std::abs(rate);
          }

      double maximum_cell_courant = 0.0;
      for (std::size_t i = 0; i < elevation.size(); ++i)
        maximum_cell_courant =
          std::max(maximum_cell_courant,
                   step_years * absolute_face_rate[i] / areas[i]);
      const unsigned int advection_steps =
        std::max(1u,
                 static_cast<unsigned int>(
                   std::ceil(maximum_cell_courant /
                             maximum_courant)));
      const double advection_step_years =
        step_years / advection_steps;

      for (unsigned int step = 0; step < advection_steps; ++step)
        {
          advect_field(bedrock_elevation, velocity,
                       advection_step_years, areas,
                       true);
          sediment_thickness.fill(0.0);
          for (unsigned int rock = 0;
               rock < lithology_names.size(); ++rock)
            {
              advect_field(sediment_thickness_by_lithology[rock], velocity,
                           advection_step_years, areas,
                           false);
              for (double &thickness :
                   sediment_thickness_by_lithology[rock])
                thickness = std::max(0.0, thickness);
              sediment_thickness += sediment_thickness_by_lithology[rock];
            }
        }
      elevation = bedrock_elevation + sediment_thickness;
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::advect_field(xt::xarray<double> &field,
                                                            const std::vector<SurfaceVelocity> &velocity,
                                                            const double step_years,
                                                            const typename Grid::container_type &areas,
                                                            const bool preserve_constant_field)
    {
      std::vector<double> extensive_change(field.size(), 0.0);
      std::vector<double> net_outward_face_rate(field.size(), 0.0);
      for (std::size_t i = 0; i < field.size(); ++i)
        for (std::size_t n = 0;
             n < grid->number_of_cell_neighbors(i); ++n)
          {
            const std::size_t j = grid->cell_neighbor(i, n);
            if (j <= i)
              continue;
            const SurfaceVelocity face_velocity =
              0.5 * (velocity[i] + velocity[j]);
            const double signed_rate =
              (face_velocity * grid->cell_neighbor_direction(i, n)) *
              grid->cell_shared_face_measure(i, n);
            net_outward_face_rate[i] += signed_rate;
            net_outward_face_rate[j] -= signed_rate;
            const std::size_t donor = signed_rate >= 0.0 ? i : j;
            const std::size_t receiver = signed_rate >= 0.0 ? j : i;
            const double transported =
              std::abs(signed_rate) * step_years * field[donor];
            extensive_change[donor] -= transported;
            extensive_change[receiver] += transported;
          }

      for (std::size_t i = 0; i < field.size(); ++i)
        {
          // Bedrock elevation is a tracer, not a conserved volume. This
          // correction changes the conservative flux form into the
          // advective form and exactly preserves a constant field even when
          // the face-interpolated velocity has small discrete divergence.
          // Sediment thickness remains in conservative flux form.
          if (preserve_constant_field)
            extensive_change[i] +=
              step_years * net_outward_face_rate[i] * field[i];
          field[i] += extensive_change[i] / areas[i];
        }
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::transport_marine_sediment(const double step_years,
                                                                         const double sea_level,
                                                                         const double maximum_courant)
    {
      unsigned int transport_steps = 1;
      if (maximum_courant > 0.0)
        {
          const auto areas = grid->nodes_areas();
          std::vector<double> conductance_sum(elevation.size(), 0.0);
          for (std::size_t i = 0; i < elevation.size(); ++i)
            if (ocean_mask[i] >= 0.5)
              for (std::size_t n = 0;
                   n < grid->number_of_cell_neighbors(i); ++n)
                {
                  const std::size_t j = grid->cell_neighbor(i, n);
                  if (j <= i || ocean_mask[j] < 0.5)
                    continue;
                  const double conductance =
                    marine_sediment_transport_coefficient *
                    grid->cell_shared_face_measure(i, n) /
                    grid->cell_neighbor_distance(i, n);
                  conductance_sum[i] += conductance;
                  conductance_sum[j] += conductance;
                }

          double largest_courant = 0.0;
          for (std::size_t i = 0; i < elevation.size(); ++i)
            largest_courant =
              std::max(largest_courant,
                       step_years * conductance_sum[i] / areas[i]);
          transport_steps = std::max(
                              1u,
                              static_cast<unsigned int>(
                                std::ceil(largest_courant / maximum_courant)));
        }

      const double transport_step_years = step_years / transport_steps;
      marine_sediment_flux.fill(0.0);
      for (unsigned int step = 0; step < transport_steps; ++step)
        transport_marine_sediment_substep(transport_step_years, sea_level);
      marine_sediment_flux /= transport_steps;
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::transport_marine_sediment_substep(const double step_years,
                                                                                 const double sea_level)
    {
      struct Transfer
      {
        std::size_t donor;
        std::size_t receiver;
        double volume;
      };

      const auto areas = grid->nodes_areas();
      std::vector<Transfer> transfers;
      std::vector<double> requested_outflow(elevation.size(), 0.0);
      for (std::size_t i = 0; i < elevation.size(); ++i)
        {
          if (ocean_mask[i] < 0.5)
            continue;

          for (std::size_t neighbor_number = 0;
               neighbor_number < grid->number_of_cell_neighbors(i);
               ++neighbor_number)
            {
              const std::size_t j =
                grid->cell_neighbor(i, neighbor_number);
              if (j <= i || ocean_mask[j] < 0.5)
                continue;

              const double elevation_difference =
                elevation[i] - elevation[j];
              const double distance =
                grid->cell_neighbor_distance(i, neighbor_number);
              if (elevation_difference == 0.0 || distance <= 0.0)
                continue;

              const std::size_t donor =
                elevation_difference > 0.0 ? i : j;
              const std::size_t receiver =
                elevation_difference > 0.0 ? j : i;
              const double mean_water_depth =
                std::max(0.0,
                         sea_level -
                         0.5 * (elevation[i] + elevation[j]));
              const double depth_factor =
                marine_transport_depth_scale > 0.0
                ? std::exp(-mean_water_depth /
                           marine_transport_depth_scale)
                : 1.0;
              const double volume =
                marine_sediment_transport_coefficient * depth_factor *
                grid->cell_shared_face_measure(i, neighbor_number) /
                distance * std::abs(elevation_difference) * step_years;

              if (volume > 0.0)
                {
                  transfers.push_back({donor, receiver, volume});
                  requested_outflow[donor] += volume;
                }
            }
        }

      std::vector<double> outflow_scale(elevation.size(), 1.0);
      for (unsigned int i = 0; i < elevation.size(); ++i)
        if (requested_outflow[i] > 0.0)
          outflow_scale[i] =
            std::min(1.0,
                     sediment_thickness[i] * areas[i] /
                     requested_outflow[i]);

      std::vector<double> volume_change(elevation.size(), 0.0);
      std::vector<std::vector<double>> lithology_volume_change(
        lithology_names.size(),
        std::vector<double>(elevation.size(), 0.0));
      for (const Transfer &transfer : transfers)
        {
          const double volume =
            transfer.volume * outflow_scale[transfer.donor];
          volume_change[transfer.donor] -= volume;
          volume_change[transfer.receiver] += volume;
          if (sediment_thickness[transfer.donor] > 0.0)
            for (unsigned int rock = 0;
                 rock < lithology_names.size(); ++rock)
              {
                const double rock_volume = volume *
                                           sediment_thickness_by_lithology[rock][transfer.donor] /
                                           sediment_thickness[transfer.donor];
                lithology_volume_change[rock][transfer.donor] -= rock_volume;
                lithology_volume_change[rock][transfer.receiver] += rock_volume;
              }
          marine_sediment_flux[transfer.donor] += volume / step_years;
          marine_sediment_flux[transfer.receiver] += volume / step_years;
        }

      for (unsigned int i = 0; i < elevation.size(); ++i)
        {
          sediment_thickness[i] += volume_change[i] / areas[i];
          sediment_thickness[i] =
            std::max(0.0, sediment_thickness[i]);
          for (unsigned int rock = 0;
               rock < lithology_names.size(); ++rock)
            sediment_thickness_by_lithology[rock][i] =
              std::max(0.0,
                       sediment_thickness_by_lithology[rock][i] +
                       lithology_volume_change[rock][i] / areas[i]);
        }
      elevation = bedrock_elevation + sediment_thickness;
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::update_drainage_diagnostics()
    {
      const auto &graph = flow_graph->impl();
      const auto &receivers = graph.receivers();
      const auto &receiver_count = graph.receivers_count();
      const auto &receiver_weight = graph.receivers_weight();
      const auto nodes = graph.nodes_indices_bottomup();

      std::vector<std::size_t> outlet_node(drainage_area.size(), 0);
      std::vector<unsigned int> basin_id(drainage_area.size(), 0);
      unsigned int next_basin_id = 1;

      // nodes_indices_bottomup() visits receivers before their donors. For
      // multiple-flow routing a unique catchment is not intrinsic, so the
      // diagnostic basin follows the receiver carrying the largest weight.
      for (const std::size_t i : nodes)
        {
          std::size_t primary_receiver = i;
          double largest_weight = -1.0;
          unsigned int nonself_receiver_count = 0;
          for (std::size_t receiver_number = 0;
               receiver_number < receiver_count[i]; ++receiver_number)
            {
              const std::size_t receiver = receivers(i, receiver_number);
              if (receiver == i)
                continue;
              ++nonself_receiver_count;
              const double weight = receiver_weight(i, receiver_number);
              if (weight > largest_weight)
                {
                  largest_weight = weight;
                  primary_receiver = receiver;
                }
            }

          if (primary_receiver == i)
            {
              outlet_node[i] = i;
              basin_id[i] = next_basin_id++;
            }
          else
            {
              Assert(basin_id[primary_receiver] != 0, ExcInternalError());
              outlet_node[i] = outlet_node[primary_receiver];
              basin_id[i] = basin_id[primary_receiver];
            }

          dominant_drainage_basin[i] = basin_id[i];
          // One-based node identifiers leave zero available for missing data.
          dominant_outlet_node[i] = outlet_node[i] + 1;
          primary_receiver_node[i] = primary_receiver + 1;
          primary_receiver_fraction[i] = primary_receiver == i
                                         ? 1.0 : largest_weight;
          flow_receiver_count[i] = nonself_receiver_count;
          drainage_outlet_mask[i] = primary_receiver == i ? 1.0 : 0.0;
          coastal_outlet_mask[i] =
            primary_receiver == i && ocean_mask[i] >= 0.5 ? 1.0 : 0.0;
        }
    }

    template <int surface_dim, int space_dim>
    void
    FastscapeLandscape<surface_dim,space_dim>::set_base_levels(const xt::xarray<double> &surface_elevation,
                                                               const double sea_level)
    {
      ocean_mask.fill(0.0);
      std::vector<std::size_t> wet_nodes;
      for (std::size_t i = 0; i < surface_elevation.size(); ++i)
        if (surface_elevation[i] <= sea_level)
          wet_nodes.push_back(i);

      std::vector<std::size_t> base_levels;
      if (!restrict_ocean_to_largest_connected_component)
        base_levels = wet_nodes;
      else if (!wet_nodes.empty())
        {
          const auto areas = grid->nodes_areas();
          std::vector<bool> wet(surface_elevation.size(), false);
          std::vector<bool> visited(surface_elevation.size(), false);
          for (const std::size_t index : wet_nodes)
            wet[index] = true;

          double largest_area = -1.0;
          for (const std::size_t seed : wet_nodes)
            if (!visited[seed])
              {
                std::vector<std::size_t> component;
                std::vector<std::size_t> frontier(1, seed);
                visited[seed] = true;
                double component_area = 0.0;
                while (!frontier.empty())
                  {
                    const std::size_t index = frontier.back();
                    frontier.pop_back();
                    component.push_back(index);
                    component_area += areas[index];
                    for (std::size_t n = 0;
                         n < grid->number_of_cell_neighbors(index);
                         ++n)
                      {
                        const std::size_t neighbor =
                          grid->cell_neighbor(index, n);
                        if (wet[neighbor] && !visited[neighbor])
                          {
                            visited[neighbor] = true;
                            frontier.push_back(neighbor);
                          }
                      }
                  }

                if (component_area > largest_area)
                  {
                    largest_area = component_area;
                    base_levels = std::move(component);
                  }
              }
        }

      if (base_levels.empty())
        base_levels.push_back(static_cast<std::size_t>(
                                std::distance(
                                  surface_elevation.begin(),
                                  std::min_element(surface_elevation.begin(),
                                                   surface_elevation.end()))));
      else
        for (const std::size_t index : base_levels)
          ocean_mask[index] = 1.0;
      if (spherical_geometry || use_sea_level_as_drainage_base_level)
        flow_graph->set_base_levels(base_levels);
    }

    template <int surface_dim, int space_dim>
    void
    SurfaceResults<surface_dim,space_dim>::initialize(const xt::xarray<double> &initial_elevation)
    {
      reference_elevation = initial_elevation;
    }

    template <int surface_dim, int space_dim>
    std::vector<double>
    SurfaceResults<surface_dim,space_dim>::get_reference_elevation() const
    {
      return std::vector<double>(reference_elevation.begin(),
                                 reference_elevation.end());
    }

    template <int surface_dim, int space_dim>
    const std::vector<std::pair<double,std::string>> &
    SurfaceResults<surface_dim,space_dim>::get_output_history() const
    {
      return output_history;
    }

    template <int surface_dim, int space_dim>
    void
    SurfaceResults<surface_dim,space_dim>::restore_output_state(
      const std::vector<double> &stored_reference_elevation,
      const std::vector<std::pair<double,std::string>> &stored_output_history)
    {
      AssertDimension(stored_reference_elevation.size(),
                      reference_elevation.size());
      std::copy(stored_reference_elevation.begin(),
                stored_reference_elevation.end(),
                reference_elevation.begin());
      output_history = stored_output_history;
      budget_output_initialized = !output_history.empty();
      basin_budget_output_initialized = !output_history.empty();
    }

    template <int surface_dim, int space_dim>
    void
    SurfaceResults<surface_dim,space_dim>::write_budget(const std::string &output_directory,
                                                        const unsigned int timestep_number,
                                                        const double time_years,
                                                        const double sea_level,
                                                        const double eroded_volume,
                                                        const double fluvial_eroded_volume,
                                                        const double glacial_eroded_volume,
                                                        const double exported_sediment_flux,
                                                        const double accommodation_limited_exported_sediment_flux,
                                                        const double coastal_sediment_flux,
                                                        const double deposited_sediment_volume,
                                                        const double stored_sediment_volume,
                                                        const FastscapeLandscape<surface_dim,space_dim> &landscape)
    {
      const std::string directory =
        output_directory + "fastscape_surface_evolution/";
      std::filesystem::create_directories(directory);

      const auto &elevation = landscape.get_elevation();
      const auto &drainage_area = landscape.get_drainage_area();
      const std::string budget_file = directory + "sediment_budget.csv";
      std::ofstream output(budget_file,
                           budget_output_initialized
                           ? std::ios::app
                           : std::ios::trunc);
      if (!budget_output_initialized)
        output << "timestep,time_years,sea_level_m,eroded_volume_m3,"
               << "fluvial_eroded_volume_m3,glacial_eroded_volume_m3,"
               << "sediment_outflux_m3_per_year,"
               << "accommodation_limited_export_m3_per_year,"
               << "coastal_sediment_flux_m3_per_year,"
               << "deposited_sediment_volume_m3,"
               << "stored_sediment_solid_volume_m3,"
               << "max_drainage_area_m2,"
               << "min_elevation_m,max_elevation_m\n";
      output << timestep_number << ','
             << std::setprecision(16) << time_years << ','
             << sea_level << ','
             << eroded_volume << ','
             << fluvial_eroded_volume << ','
             << glacial_eroded_volume << ','
             << exported_sediment_flux << ','
             << accommodation_limited_exported_sediment_flux << ','
             << coastal_sediment_flux << ','
             << deposited_sediment_volume << ','
             << stored_sediment_volume << ','
             << *std::max_element(drainage_area.begin(),
                                  drainage_area.end()) << ','
             << *std::min_element(elevation.begin(), elevation.end()) << ','
             << *std::max_element(elevation.begin(), elevation.end()) << '\n';
      budget_output_initialized = true;
    }

    template <int surface_dim, int space_dim>
    void
    SurfaceResults<surface_dim,space_dim>::write(const std::string &output_directory,
                                                 const unsigned int timestep_number,
                                                 const double time_years,
                                                 const bool write_visualization,
                                                 const bool output_drainage_diagnostics,
                                                 const SurfaceMesh &surface_mesh,
                                                 const std::vector<Point<space_dim>> &surface_points,
                                                 FastscapeLandscape<surface_dim,space_dim> &landscape,
                                                 const xt::xarray<double> &erosion_strength,
                                                 const xt::xarray<double> &surface_runoff,
                                                 const xt::xarray<double> &,
                                                 const xt::xarray<double> &,
                                                 const std::vector<double> &regional_ice_load_displacement,
                                                 const std::vector<double> &regional_ice_load_velocity)
    {
      const std::string directory =
        output_directory + "fastscape_surface_evolution/";
      std::filesystem::create_directories(directory);

      const auto &elevation = landscape.get_elevation();
      const auto &drainage_area = landscape.get_drainage_area();
      const auto &geometric_drainage_area =
        landscape.get_geometric_drainage_area();
      const auto &dominant_drainage_basin =
        landscape.get_dominant_drainage_basin();
      const auto &dominant_outlet_node =
        landscape.get_dominant_outlet_node();
      const auto &primary_receiver_node =
        landscape.get_primary_receiver_node();
      const auto &primary_receiver_fraction =
        landscape.get_primary_receiver_fraction();
      const auto &flow_receiver_count =
        landscape.get_flow_receiver_count();
      const auto &drainage_outlet_mask =
        landscape.get_drainage_outlet_mask();
      const auto &coastal_outlet_mask =
        landscape.get_coastal_outlet_mask();
      const auto &erosion = landscape.get_erosion();
      const auto &fluvial_erosion = landscape.get_fluvial_erosion();
      const auto &glacial_erosion = landscape.get_glacial_erosion();
      const auto &accumulated_fluvial_erosion =
        landscape.get_accumulated_fluvial_erosion();
      const auto &accumulated_glacial_erosion =
        landscape.get_accumulated_glacial_erosion();
      const auto cell_areas = landscape.get_cell_areas();
      const auto &modeled_ice_thickness =
        landscape.get_modeled_ice_thickness();
      const auto &modeled_basal_ice_velocity =
        landscape.get_modeled_basal_ice_velocity();
      const auto &routed_ice_discharge =
        landscape.get_routed_ice_discharge();
      const auto &routed_glacier_width =
        landscape.get_routed_glacier_width();
      const auto &routed_ice_mass_balance =
        landscape.get_routed_ice_mass_balance();
      const auto &sediment_flux = landscape.get_sediment_flux();
      const auto &marine_sediment_flux =
        landscape.get_marine_sediment_flux();
      const auto &sediment_thickness =
        landscape.get_sediment_thickness();
      const auto &deposition_rate =
        landscape.get_deposition_rate();
      const auto &ocean_mask = landscape.get_ocean_mask();
      const auto &lithology_names = landscape.get_lithology_names();
      const auto &bedrock_lithology = landscape.get_bedrock_lithology();
      const auto &sediment_flux_by_lithology =
        landscape.get_sediment_flux_by_lithology();
      const auto &sediment_thickness_by_lithology =
        landscape.get_sediment_thickness_by_lithology();
      const auto deposited_thickness_by_lithology =
        landscape.take_deposited_thickness_by_lithology();
      AssertDimension(regional_ice_load_displacement.size(), elevation.size());
      AssertDimension(regional_ice_load_velocity.size(), elevation.size());

      const std::string surface_file =
        directory + "surface-" +
        Utilities::int_to_string(timestep_number, 5) + ".csv";
      std::ofstream surface(surface_file);
      surface << "longitude_deg,latitude_deg,surface_x_m,surface_y_m,"
              << "elevation_m,erosion_m,"
              << "fluvial_erosion_m,glacial_erosion_m,"
              << "drainage_area_m2";
      if (output_drainage_diagnostics)
        surface << ",runoff_weighted_drainage_area_m2,"
                << "geometric_drainage_area_m2,"
                << "dominant_drainage_basin_id,"
                << "dominant_outlet_node_one_based,"
                << "primary_receiver_node_one_based,"
                << "primary_receiver_fraction,flow_receiver_count,"
                << "is_drainage_outlet,is_coastal_outlet";
      surface << ",sediment_flux_m3_per_year,"
              << "marine_sediment_flux_m3_per_year,"
              << "sediment_thickness_m,deposition_rate_m_per_year,"
              << "is_connected_ocean,bedrock_lithology,"
              << "erosion_strength,surface_runoff_factor,"
              << "ice_thickness_m,basal_ice_velocity_m_per_year,"
              << "routed_ice_discharge_m3_per_year,"
              << "routed_glacier_width_m,"
              << "routed_ice_mass_balance_m_per_year,"
              << "regional_ice_load_displacement_m,"
              << "regional_ice_load_velocity_m_per_year,"
              << "elevation_change_m";
      for (const std::string &name : lithology_names)
        surface << ",sediment_flux_" << name << "_m3_per_year"
                << ",sediment_thickness_" << name << "_m";
      surface << '\n';
      surface << std::setprecision(16);
      for (unsigned int i = 0; i < elevation.size(); ++i)
        {
          const Point<space_dim> &point = surface_points[i];
          const double longitude =
            std::atan2(point[1], point[0]) * 180.0 / numbers::PI;
          double latitude = 0.0;
          if constexpr (space_dim == 3)
            latitude =
              std::asin(point[2] / point.norm()) * 180.0 / numbers::PI;
          surface << longitude << ',' << latitude << ','
                  << point[0] << ',' << point[1] << ',' << elevation[i] << ','
                  << erosion[i] << ',' << fluvial_erosion[i] << ','
                  << glacial_erosion[i] << ',' << drainage_area[i];
          if (output_drainage_diagnostics)
            surface << ',' << drainage_area[i]
                    << ',' << geometric_drainage_area[i]
                    << ',' << dominant_drainage_basin[i]
                    << ',' << dominant_outlet_node[i]
                    << ',' << primary_receiver_node[i]
                    << ',' << primary_receiver_fraction[i]
                    << ',' << flow_receiver_count[i]
                    << ',' << drainage_outlet_mask[i]
                    << ',' << coastal_outlet_mask[i];
          surface << ',' << sediment_flux[i] << ',' << marine_sediment_flux[i] << ','
                  << sediment_thickness[i] << ',' << deposition_rate[i] << ','
                  << ocean_mask[i] << ','
                  << lithology_names[bedrock_lithology[i]] << ','
                  << erosion_strength[i] << ','
                  << surface_runoff[i] << ','
                  << modeled_ice_thickness[i] << ','
                  << modeled_basal_ice_velocity[i] << ','
                  << routed_ice_discharge[i] << ','
                  << routed_glacier_width[i] << ','
                  << routed_ice_mass_balance[i] << ','
                  << regional_ice_load_displacement[i] << ','
                  << regional_ice_load_velocity[i] << ','
                  << elevation[i] - reference_elevation[i];
          for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
            surface << ',' << sediment_flux_by_lithology[rock][i]
                    << ',' << sediment_thickness_by_lithology[rock][i];
          surface << '\n';
        }

      if (output_drainage_diagnostics)
        {
          struct BasinBudget
          {
            std::size_t outlet_node = 0;
            double surface_area = 0.0;
            double receiver_fraction_area_sum = 0.0;
            double fluvial_eroded_volume = 0.0;
            double glacial_eroded_volume = 0.0;
            double land_net_deposition_rate = 0.0;
            double marine_net_deposition_rate = 0.0;
            double stored_sediment_bulk_volume = 0.0;
          };
          std::map<unsigned int,BasinBudget> basin_budgets;
          for (unsigned int i = 0; i < elevation.size(); ++i)
            {
              const unsigned int basin = static_cast<unsigned int>(
                                           std::lround(dominant_drainage_basin[i]));
              BasinBudget &budget = basin_budgets[basin];
              budget.outlet_node = static_cast<std::size_t>(
                                     std::lround(dominant_outlet_node[i])) - 1;
              budget.surface_area += cell_areas[i];
              budget.receiver_fraction_area_sum +=
                primary_receiver_fraction[i] * cell_areas[i];
              budget.fluvial_eroded_volume +=
                accumulated_fluvial_erosion[i] * cell_areas[i];
              budget.glacial_eroded_volume +=
                accumulated_glacial_erosion[i] * cell_areas[i];
              const double net_deposition_rate =
                deposition_rate[i] * cell_areas[i];
              if (ocean_mask[i] >= 0.5)
                budget.marine_net_deposition_rate += net_deposition_rate;
              else
                budget.land_net_deposition_rate += net_deposition_rate;
              budget.stored_sediment_bulk_volume +=
                sediment_thickness[i] * cell_areas[i];
            }

          const std::string basin_budget_file =
            directory + "drainage_basin_budget.csv";
          const bool basin_budget_file_exists =
            std::filesystem::exists(basin_budget_file) &&
            std::filesystem::file_size(basin_budget_file) > 0;
          std::ofstream basin_output(
            basin_budget_file,
            basin_budget_output_initialized && basin_budget_file_exists
            ? std::ios::app : std::ios::trunc);
          if (!basin_budget_output_initialized || !basin_budget_file_exists)
            basin_output
                << "timestep,time_years,dominant_drainage_basin_id,"
                << "outlet_node_one_based,outlet_x_m,outlet_y_m,"
                << "outlet_elevation_m,is_coastal_outlet,"
                << "dominant_basin_surface_area_m2,"
                << "outlet_geometric_drainage_area_m2,"
                << "outlet_runoff_weighted_drainage_area_m2,"
                << "area_weighted_mean_primary_receiver_fraction,"
                << "fluvial_eroded_volume_m3_geodynamic_step,"
                << "glacial_eroded_volume_m3_geodynamic_step,"
                << "land_net_deposition_rate_m3_per_year,"
                << "marine_net_deposition_rate_m3_per_year,"
                << "sediment_flux_at_outlet_m3_per_year,"
                << "marine_sediment_flux_at_outlet_m3_per_year,"
                << "stored_sediment_bulk_volume_m3\n";
          basin_output << std::setprecision(16);
          for (const auto &[basin,budget] : basin_budgets)
            {
              const std::size_t outlet = budget.outlet_node;
              basin_output
                  << timestep_number << ',' << time_years << ',' << basin << ','
                  << outlet + 1 << ',' << surface_points[outlet][0] << ','
                  << surface_points[outlet][1] << ',' << elevation[outlet] << ','
                  << coastal_outlet_mask[outlet] << ',' << budget.surface_area << ','
                  << geometric_drainage_area[outlet] << ','
                  << drainage_area[outlet] << ','
                  << budget.receiver_fraction_area_sum / budget.surface_area << ','
                  << budget.fluvial_eroded_volume << ','
                  << budget.glacial_eroded_volume << ','
                  << budget.land_net_deposition_rate << ','
                  << budget.marine_net_deposition_rate << ','
                  << sediment_flux[outlet] << ','
                  << marine_sediment_flux[outlet] << ','
                  << budget.stored_sediment_bulk_volume << '\n';
            }
          basin_budget_output_initialized = true;
        }

      // Each output interval is one dated depositional layer. Rows contain
      // sediment preserved since the preceding result, split by source-rock
      // class; together the files form a basin stratigraphy.
      const std::string stratigraphy_file =
        directory + "stratigraphy-" +
        Utilities::int_to_string(timestep_number, 5) + ".csv";
      std::ofstream stratigraphy(stratigraphy_file);
      stratigraphy << "longitude_deg,latitude_deg,time_years,"
                   << "bedrock_lithology,deposited_thickness_m";
      for (const std::string &name : lithology_names)
        stratigraphy << ",fraction_" << name;
      stratigraphy << '\n' << std::setprecision(16);
      for (unsigned int i = 0; i < elevation.size(); ++i)
        {
          double deposited_thickness = 0.0;
          for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
            deposited_thickness += deposited_thickness_by_lithology[rock][i];
          if (deposited_thickness <= 0.0)
            continue;
          const Point<space_dim> &point = surface_points[i];
          const double longitude =
            std::atan2(point[1], point[0]) * 180.0 / numbers::PI;
          double latitude = 0.0;
          if constexpr (space_dim == 3)
            latitude =
              std::asin(point[2] / point.norm()) * 180.0 / numbers::PI;
          stratigraphy << longitude << ',' << latitude << ',' << time_years
                       << ',' << lithology_names[bedrock_lithology[i]]
                       << ',' << deposited_thickness;
          for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
            stratigraphy << ','
                         << deposited_thickness_by_lithology[rock][i] /
                         deposited_thickness;
          stratigraphy << '\n';
        }

      if (!write_visualization)
        return;

      Vector<double> elevation_output(elevation.size());
      Vector<double> erosion_output(erosion.size());
      Vector<double> fluvial_erosion_output(fluvial_erosion.size());
      Vector<double> glacial_erosion_output(glacial_erosion.size());
      Vector<double> drainage_output(drainage_area.size());
      Vector<double> geometric_drainage_output(drainage_area.size());
      Vector<double> basin_output(drainage_area.size());
      Vector<double> outlet_output(drainage_area.size());
      Vector<double> receiver_output(drainage_area.size());
      Vector<double> receiver_fraction_output(drainage_area.size());
      Vector<double> receiver_count_output(drainage_area.size());
      Vector<double> drainage_outlet_output(drainage_area.size());
      Vector<double> coastal_outlet_output(drainage_area.size());
      Vector<double> flux_output(sediment_flux.size());
      Vector<double> marine_flux_output(marine_sediment_flux.size());
      Vector<double> sediment_thickness_output(sediment_thickness.size());
      Vector<double> deposition_rate_output(deposition_rate.size());
      Vector<double> ocean_mask_output(ocean_mask.size());
      Vector<double> ice_thickness_output(modeled_ice_thickness.size());
      Vector<double> basal_ice_velocity_output(modeled_basal_ice_velocity.size());
      Vector<double> routed_ice_discharge_output(routed_ice_discharge.size());
      Vector<double> routed_glacier_width_output(routed_glacier_width.size());
      Vector<double> routed_ice_mass_balance_output(routed_ice_mass_balance.size());
      Vector<double> regional_displacement_output(elevation.size());
      Vector<double> regional_velocity_output(elevation.size());
      Vector<double> bedrock_lithology_output(elevation.size());
      for (unsigned int i = 0; i < elevation.size(); ++i)
        {
          elevation_output[i] = elevation[i];
          erosion_output[i] = erosion[i];
          fluvial_erosion_output[i] = fluvial_erosion[i];
          glacial_erosion_output[i] = glacial_erosion[i];
          drainage_output[i] = drainage_area[i];
          geometric_drainage_output[i] = geometric_drainage_area[i];
          basin_output[i] = dominant_drainage_basin[i];
          outlet_output[i] = dominant_outlet_node[i];
          receiver_output[i] = primary_receiver_node[i];
          receiver_fraction_output[i] = primary_receiver_fraction[i];
          receiver_count_output[i] = flow_receiver_count[i];
          drainage_outlet_output[i] = drainage_outlet_mask[i];
          coastal_outlet_output[i] = coastal_outlet_mask[i];
          flux_output[i] = sediment_flux[i];
          marine_flux_output[i] = marine_sediment_flux[i];
          sediment_thickness_output[i] = sediment_thickness[i];
          deposition_rate_output[i] = deposition_rate[i];
          ocean_mask_output[i] = ocean_mask[i];
          ice_thickness_output[i] = modeled_ice_thickness[i];
          basal_ice_velocity_output[i] = modeled_basal_ice_velocity[i];
          routed_ice_discharge_output[i] = routed_ice_discharge[i];
          routed_glacier_width_output[i] = routed_glacier_width[i];
          routed_ice_mass_balance_output[i] = routed_ice_mass_balance[i];
          regional_displacement_output[i] =
            regional_ice_load_displacement[i];
          regional_velocity_output[i] = regional_ice_load_velocity[i];
          bedrock_lithology_output[i] = bedrock_lithology[i];
        }

      DataOut<surface_dim,space_dim> data_out;
      data_out.attach_triangulation(surface_mesh);
      data_out.add_data_vector(elevation_output, "elevation",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(erosion_output, "erosion",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(fluvial_erosion_output, "fluvial_erosion",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(glacial_erosion_output, "glacial_erosion",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(drainage_output, "drainage_area",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      if (output_drainage_diagnostics)
        {
          data_out.add_data_vector(drainage_output,
                                   "runoff_weighted_drainage_area",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(geometric_drainage_output,
                                   "geometric_drainage_area",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(basin_output,
                                   "dominant_drainage_basin_id",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(outlet_output,
                                   "dominant_outlet_node_one_based",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(receiver_output,
                                   "primary_receiver_node_one_based",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(receiver_fraction_output,
                                   "primary_receiver_fraction",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(receiver_count_output,
                                   "flow_receiver_count",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(drainage_outlet_output,
                                   "is_drainage_outlet",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
          data_out.add_data_vector(coastal_outlet_output,
                                   "is_coastal_outlet",
                                   DataOut<surface_dim,space_dim>::type_cell_data);
        }
      data_out.add_data_vector(flux_output, "sediment_flux",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(marine_flux_output, "marine_sediment_flux",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(sediment_thickness_output, "sediment_thickness",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(deposition_rate_output, "deposition_rate",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(ocean_mask_output, "is_connected_ocean",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(ice_thickness_output, "ice_thickness",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(basal_ice_velocity_output,
                               "basal_ice_velocity",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(routed_ice_discharge_output,
                               "routed_ice_discharge",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(routed_glacier_width_output,
                               "routed_glacier_width",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(routed_ice_mass_balance_output,
                               "routed_ice_mass_balance",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(regional_displacement_output,
                               "regional_ice_load_displacement",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(regional_velocity_output,
                               "regional_ice_load_velocity",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      data_out.add_data_vector(bedrock_lithology_output,
                               "bedrock_lithology",
                               DataOut<surface_dim,space_dim>::type_cell_data);
      for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
        {
          Vector<double> rock_fraction(elevation.size());
          for (unsigned int i = 0; i < elevation.size(); ++i)
            rock_fraction[i] = sediment_thickness[i] > 0.0
                               ? sediment_thickness_by_lithology[rock][i] /
                               sediment_thickness[i]
                               : 0.0;
          data_out.add_data_vector(rock_fraction,
                                   "sediment_fraction_" +
                                   lithology_names[rock],
                                   DataOut<surface_dim,space_dim>::type_cell_data);
        }
      data_out.build_patches();

      const std::string basename =
        "fastscape-" + Utilities::int_to_string(timestep_number, 5) + ".vtu";
      std::ofstream visualization_file(directory + basename);
      data_out.write_vtu(visualization_file);

      output_history.emplace_back(time_years, basename);
      std::ofstream collection_file(directory + "fastscape.pvd");
      DataOutBase::write_pvd_record(collection_file, output_history);
    }


    template class SpatialErosionStrength<1>;
    template class SpatialErosionStrength<2>;
    template class SpatialSurfaceRunoff<1>;
    template class SpatialSurfaceRunoff<2>;
    template class SpatialIceThickness<1>;
    template class SpatialIceThickness<2>;
    template class SpatialBasalIceVelocity<1>;
    template class SpatialBasalIceVelocity<2>;
    template class FastscapeLandscape<1,2>;
    template class FastscapeLandscape<2,3>;
    template class SurfaceResults<1,2>;
    template class SurfaceResults<2,3>;

  }
}

#endif
