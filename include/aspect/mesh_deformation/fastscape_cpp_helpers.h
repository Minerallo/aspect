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

#ifndef _aspect_mesh_deformation_fastscape_cpp_helpers_h
#define _aspect_mesh_deformation_fastscape_cpp_helpers_h

#include <aspect/config.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <aspect/mesh_deformation/fastscape_cpp_grid.h>
#include <aspect/structured_data.h>

#include <fastscapelib/eroders/spl.hpp>
#include <fastscapelib/flow/flow_graph.hpp>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

/**
 * @file
 * Declare helper classes used internally by the FastscapeCpp mesh deformation
 * plugin. Implementations are in
 * source/mesh_deformation/fastscape_cpp_helpers.cc.
 */

namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * Read a generic erosion-strength field and interpolate it onto the
     * independent landscape mesh. This class deliberately has no knowledge of
     * the program or physical formula that produced the field.
     */
    template <int surface_dim>
    class SpatialErosionStrength
    {
      public:
        void
        initialize(const std::string &filename,
                   const std::vector<Point<surface_dim>> &surface_coordinates);

        void
        sample(const std::vector<Point<surface_dim>> &surface_coordinates);

        const xt::xarray<double> &
        get_values() const;

      private:
        xt::xarray<double> values;
        std::unique_ptr<Utilities::StructuredDataLookup<surface_dim>> lookup;
    };

    /**
     * Read a dimensionless local surface-runoff field and interpolate it onto
     * the independent landscape mesh. FastScape multiplies each surface-cell
     * area by this value before accumulating water supply downstream. A
     * uniform value of one therefore recovers ordinary drainage area.
     */
    template <int surface_dim>
    class SpatialSurfaceRunoff
    {
      public:
        void
        initialize(const std::string &filename,
                   const std::vector<Point<surface_dim>> &surface_coordinates);

        void
        sample(const std::vector<Point<surface_dim>> &surface_coordinates);

        const xt::xarray<double> &
        get_values() const;

      private:
        xt::xarray<double> values;
        std::unique_ptr<Utilities::StructuredDataLookup<surface_dim>> lookup;
    };


    /**
     * Read ice thickness in meters and interpolate it onto the independent
     * landscape mesh. An omitted file represents an ice-free surface.
     */
    template <int surface_dim>
    class SpatialIceThickness
    {
      public:
        void
        initialize(const std::string &filename,
                   const std::vector<Point<surface_dim>> &surface_coordinates);

        void
        sample(const std::vector<Point<surface_dim>> &surface_coordinates);

        const xt::xarray<double> &
        get_values() const;

      private:
        xt::xarray<double> values;
        std::unique_ptr<Utilities::StructuredDataLookup<surface_dim>> lookup;
    };


    /**
     * Read basal ice velocity in meters per year and interpolate it onto the
     * independent landscape mesh. An omitted file represents no sliding.
     */
    template <int surface_dim>
    class SpatialBasalIceVelocity
    {
      public:
        void
        initialize(const std::string &filename,
                   const std::vector<Point<surface_dim>> &surface_coordinates);

        void
        sample(const std::vector<Point<surface_dim>> &surface_coordinates);

        const xt::xarray<double> &
        get_values() const;

      private:
        xt::xarray<double> values;
        std::unique_ptr<Utilities::StructuredDataLookup<surface_dim>> lookup;
    };


    /**
     * Own the FastScape grid, routing graph, erosion solver, and evolving
     * landscape fields. ASPECT-specific field transfer and file output remain
     * outside this class.
     */
    template <int surface_dim, int space_dim>
    class FastscapeLandscape
    {
      public:
        using SurfaceMesh = Triangulation<surface_dim,space_dim>;
        using Grid = fastscapelib::dealii_surface_grid<SurfaceMesh>;
        using FlowGraph = fastscapelib::flow_graph<Grid>;
        using Eroder = fastscapelib::spl_eroder<FlowGraph>;
        using SurfaceVelocity = Tensor<1,space_dim>;

        struct StepResult
        {
          xt::xarray<double> previous_elevation;
          double eroded_volume = 0.0;
          double fluvial_eroded_volume = 0.0;
          double glacial_eroded_volume = 0.0;
          double exported_sediment_flux = 0.0;
          double accommodation_limited_exported_sediment_flux = 0.0;
          double coastal_sediment_flux = 0.0;
          double deposited_sediment_volume = 0.0;
          double stored_sediment_volume = 0.0;
        };

        void
        initialize(SurfaceMesh &surface_mesh,
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
                   const unsigned int rock_random_seed);

        StepResult
        advance(const xt::xarray<double> &uplift_rate,
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
                const double maximum_diffusion_courant);

        const xt::xarray<double> &get_elevation() const;

        const xt::xarray<double> &get_drainage_area() const;

        const xt::xarray<double> &get_geometric_drainage_area() const;

        const xt::xarray<double> &get_dominant_drainage_basin() const;

        const xt::xarray<double> &get_dominant_outlet_node() const;

        const xt::xarray<double> &get_primary_receiver_node() const;

        const xt::xarray<double> &get_primary_receiver_fraction() const;

        const xt::xarray<double> &get_flow_receiver_count() const;

        const xt::xarray<double> &get_drainage_outlet_mask() const;

        const xt::xarray<double> &get_coastal_outlet_mask() const;

        xt::xarray<double> get_cell_areas() const;

        const xt::xarray<double> &get_erosion() const;

        const xt::xarray<double> &get_fluvial_erosion() const;

        const xt::xarray<double> &get_glacial_erosion() const;

        const xt::xarray<double> &get_accumulated_fluvial_erosion() const;

        const xt::xarray<double> &get_accumulated_glacial_erosion() const;

        const xt::xarray<double> &get_modeled_ice_thickness() const;

        const xt::xarray<double> &get_modeled_basal_ice_velocity() const;

        const xt::xarray<double> &get_routed_ice_discharge() const;

        const xt::xarray<double> &get_routed_glacier_width() const;

        const xt::xarray<double> &get_routed_ice_mass_balance() const;

        const xt::xarray<double> &get_sediment_flux() const;

        const xt::xarray<double> &get_marine_sediment_flux() const;

        const xt::xarray<double> &get_sediment_thickness() const;

        const std::vector<std::string> &get_lithology_names() const;

        const std::vector<unsigned int> &get_bedrock_lithology() const;

        const std::vector<xt::xarray<double>> &
        get_sediment_flux_by_lithology() const;

        const std::vector<xt::xarray<double>> &
        get_sediment_thickness_by_lithology() const;

        std::vector<xt::xarray<double>>
        take_deposited_thickness_by_lithology();

        const xt::xarray<double> &get_deposition_rate() const;

        const xt::xarray<double> &get_ocean_mask() const;

        void
        set_elevation(const std::vector<double> &values);

        void
        set_sediment_thickness(const std::vector<double> &values);


        void
        set_lithology_state(
          const std::vector<unsigned int> &stored_bedrock_lithology,
          const std::vector<std::vector<double>> &stored_sediment_thickness);

      private:
        /**
         * Build a deliberately inexpensive glacier proxy on the current drainage
         * graph. Positive degree-day-style mass balance above the ELA is routed
         * downstream and is progressively removed by ablation below the ELA.
         * Width and thickness follow configurable discharge power laws; basal
         * speed then follows directly from Q = u H W. This is not an ice-dynamics
         * solver, but it produces evolving glacier corridors that respond to
         * topography and conserve the routed ice flux.
         */
        void
        update_routed_glacier_fields(
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
          const bool terminate_at_sea_level);

        /**
         * Spread centerline glacial erosion over the empirical glacier width.
         * A compact Gaussian graph-distance kernel conserves eroded volume and
         * provides valley widening without solving continuum ice stresses.
         */
        xt::xarray<double>
        spread_routed_glacial_erosion(
          const xt::xarray<double> &centerline_erosion,
          const xt::xarray<double> &glacier_width) const;

        /**
         * Apply conservative linear hillslope transport on the unstructured
         * surface grid. Material removed from a high cell is first taken from
         * mobile sediment and then bedrock; deposition becomes mobile sediment.
         */
        void
        diffuse_hillslopes(const double step_years,
                           const double diffusivity,
                           const double submarine_diffusivity,
                           const double maximum_courant);

        void
        export_marine_sediment(StepResult &result);

        /**
         * Advect bedrock elevation as a tracer and mobile-sediment thickness as a
         * conserved volume over the fixed landscape grid using a first-order
         * upwind finite-volume scheme. The velocity is tangential to the ASPECT
         * surface and expressed in meters per year. Courant substeps keep the
         * mobile thickness positive.
         */
        void
        advect_surface_fields(const std::vector<SurfaceVelocity> &velocity,
                              const double step_years,
                              const double maximum_courant);

        void
        advect_field(xt::xarray<double> &field,
                     const std::vector<SurfaceVelocity> &velocity,
                     const double step_years,
                     const typename Grid::container_type &areas,
                     const bool preserve_constant_field);

        /**
         * Move deposited sediment down the seafloor gradient. Each shared face
         * is visited once, and equal volumes are removed from the donor and
         * added to the receiver. A donor-wide limiter prevents transport from
         * removing more sediment than is locally available.
         */
        void
        transport_marine_sediment(const double step_years,
                                  const double sea_level,
                                  const double maximum_courant);

        void
        transport_marine_sediment_substep(const double step_years,
                                          const double sea_level);

        void
        update_drainage_diagnostics();

        void
        set_base_levels(const xt::xarray<double> &surface_elevation,
                        const double sea_level);

        bool spherical_geometry = false;
        double drainage_area_exponent = 0.4;
        double marine_sediment_transport_coefficient = 0.0;
        double marine_sediment_porosity = 0.4;
        double marine_transport_depth_scale = 0.0;
        bool limit_marine_deposition_to_available_accommodation = false;
        double maximum_marine_deposition_above_sea_level = 0.0;
        double submarine_river_incision_factor = 1.0;
        bool use_sea_level_as_drainage_base_level = false;
        std::string open_marine_sediment_boundary = "none";
        bool restrict_ocean_to_largest_connected_component = true;
        std::unique_ptr<Grid> grid;
        std::unique_ptr<FlowGraph> flow_graph;
        std::unique_ptr<Eroder> eroder;
        xt::xarray<double> elevation;
        xt::xarray<double> bedrock_elevation;
        xt::xarray<double> drainage_area;
        xt::xarray<double> geometric_drainage_area;
        xt::xarray<double> unit_surface_runoff;
        xt::xarray<double> dominant_drainage_basin;
        xt::xarray<double> dominant_outlet_node;
        xt::xarray<double> primary_receiver_node;
        xt::xarray<double> primary_receiver_fraction;
        xt::xarray<double> flow_receiver_count;
        xt::xarray<double> drainage_outlet_mask;
        xt::xarray<double> coastal_outlet_mask;
        xt::xarray<double> erosion;
        xt::xarray<double> fluvial_erosion;
        xt::xarray<double> glacial_erosion;
        xt::xarray<double> accumulated_fluvial_erosion;
        xt::xarray<double> accumulated_glacial_erosion;
        xt::xarray<double> modeled_ice_thickness;
        xt::xarray<double> modeled_basal_ice_velocity;
        xt::xarray<double> routed_ice_discharge;
        xt::xarray<double> routed_glacier_width;
        xt::xarray<double> routed_ice_mass_balance;
        xt::xarray<double> sediment_flux;
        xt::xarray<double> marine_sediment_flux;
        xt::xarray<double> sediment_thickness;
        xt::xarray<double> deposition_rate;
        xt::xarray<double> ocean_mask;
        std::vector<std::string> lithology_names;
        std::vector<double> lithology_erodibility_factors;
        std::vector<unsigned int> bedrock_lithology;
        std::vector<xt::xarray<double>> sediment_flux_by_lithology;
        std::vector<xt::xarray<double>> sediment_thickness_by_lithology;
        std::vector<xt::xarray<double>>
        sediment_thickness_at_last_output_by_lithology;
    };


    /**
     * Write landscape budgets and spatial result files. Keeping this code out
     * of the registered model makes the coupling calculation easier to follow.
     */
    template <int surface_dim, int space_dim>
    class SurfaceResults
    {
      public:
        using SurfaceMesh = Triangulation<surface_dim,space_dim>;

        void
        initialize(const xt::xarray<double> &initial_elevation);

        std::vector<double>
        get_reference_elevation() const;

        const std::vector<std::pair<double,std::string>> &
        get_output_history() const;

        void
        restore_output_state(
          const std::vector<double> &stored_reference_elevation,
          const std::vector<std::pair<double,std::string>> &stored_output_history);

        void
        write_budget(const std::string &output_directory,
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
                     const FastscapeLandscape<surface_dim,space_dim> &landscape);

        void
        write(const std::string &output_directory,
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
              const std::vector<double> &regional_ice_load_velocity);

      private:
        xt::xarray<double> reference_elevation;
        std::vector<std::pair<double,std::string>> output_history;
        bool budget_output_initialized = false;
        bool basin_budget_output_initialized = false;
    };


  }
}

#endif
#endif
