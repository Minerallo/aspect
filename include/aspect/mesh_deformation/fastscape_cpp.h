/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#ifndef _aspect_mesh_deformation_fastscape_cpp_h
#define _aspect_mesh_deformation_fastscape_cpp_h

#include <aspect/config.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <aspect/mesh_deformation/fastscape_cpp_grid.h>
#include <aspect/mesh_deformation/parallel_unstructured_interface.h>

#include <fastscapelib/eroders/spl.hpp>
#include <fastscapelib/flow/flow_graph.hpp>

#include <deal.II/base/parsed_function.h>
#include <deal.II/grid/tria.h>

#include <map>
#include <memory>
#include <string>


namespace aspect
{
namespace MeshDeformation
{
template <int dim>
class FastscapeLandscape;

template <int dim>
class SpatialErosionStrength;

template <int dim>
class SpatialSurfaceRunoff;

template <int dim>
class SpatialIceThickness;

template <int dim>
class SpatialBasalIceVelocity;

template <int dim>
class SurfaceResults;

/**
 * Use the FastScape library to evolve an independent surface mesh.
 *
 * The plugin uses the ParallelUnstructuredInterface introduced in
 * geodynamics/aspect#7083 for MPI-safe transfer between ASPECT's volume
 * mesh and an independent FastScape surface grid. The FastScape grid is
 * a planar surface mesh for Box geometries and a closed surface mesh for
 * SphericalShell geometries. FastScape runs on the first process; the
 * transfer framework distributes ASPECT values to its evaluation points
 * and transfers the resulting normal velocities back to the surface.
 */
template <int dim>
class FastscapeCpp : public ParallelUnstructuredInterface<dim>
{
public:
    FastscapeCpp();
    ~FastscapeCpp() override;

    void initialize() override;
    void update() override;

    std::vector<Tensor<1,dim>>
                            compute_updated_velocities_at_points(
                                const std::vector<std::vector<double>> &current_solution_at_points) const override;

    bool needs_surface_stabilization() const override;

    void save(std::map<std::string, std::string> &status_strings) const override;
    void load(const std::map<std::string, std::string> &status_strings) override;

    static void declare_parameters(ParameterHandler &prm);
    void parse_parameters(ParameterHandler &prm) override;

private:
    using SurfaceMesh = Triangulation<dim-1,dim>;

    void build_surface_mesh();
    Point<dim> reference_surface_point(const Point<dim> &point) const;
    Point<dim-1> natural_surface_coordinates(const Point<dim> &point) const;
    Point<dim-1> climate_surface_coordinates(const Point<dim> &point) const;
    Tensor<1,dim> outward_direction(const Point<dim> &point) const;
    void update_true_polar_wander();
    void resample_climate_fields();
    SymmetricTensor<2,dim> ice_load_moment_of_inertia() const;
    SymmetricTensor<2,dim> apply_degree_two_self_gravity(
      const SymmetricTensor<2,dim> &rigid_ice_load);
    std::vector<double> update_regional_ice_load_response(
      const double time_step_years) const;
    void write_true_polar_wander_state() const;

    mutable SurfaceMesh surface_mesh;
    std::vector<Point<dim>> fastscape_points;
    std::vector<double> fastscape_point_areas;
    mutable std::unique_ptr<FastscapeLandscape<dim>> landscape;
    std::unique_ptr<SpatialErosionStrength<dim>> spatial_erosion_strength;
    std::unique_ptr<SpatialSurfaceRunoff<dim>> spatial_surface_runoff;
    std::unique_ptr<SpatialIceThickness<dim>> spatial_ice_thickness;
    std::unique_ptr<SpatialBasalIceVelocity<dim>> spatial_basal_ice_velocity;
    mutable std::unique_ptr<SurfaceResults<dim>> surface_results;

    unsigned int box_repetitions = 8;
    unsigned int surface_refinement = 2;
    std::string surface_transfer_scheme = "conservative";
    unsigned int surface_transfer_neighbors = 8;
    unsigned int landscape_steps_per_geodynamic_step = 4;
    double maximum_landscape_step_years = 10000.0;
    double incision_rate = 5e-5;
    double drainage_area_exponent = 0.4;
    double slope_exponent = 1.0;
    double nonlinear_tolerance = 1e-5;
    double glacial_erosion_coefficient = 0.0;
    double glacial_velocity_exponent = 1.0;
    double minimum_ice_thickness = 1.0;
    double initial_relief = 0.0;
    double sea_level = 0.0;
    Functions::ParsedFunction<1> sea_level_function;
    double marine_sediment_transport_coefficient = 0.0;
    double marine_sediment_porosity = 0.4;
    double marine_transport_depth_scale = 0.0;
    double maximum_surface_advection_courant = 0.5;
    double hillslope_diffusion_coefficient = 0.0;
    double maximum_hillslope_diffusion_courant = 0.25;
    std::vector<std::string> lithology_names = {"upper_crust"};
    std::vector<double> lithology_probabilities = {1.0};
    std::vector<double> lithology_erodibility_factors = {1.0};
    unsigned int lithology_random_seed = 1;
    std::string spatial_erosion_strength_file;
    std::string spatial_surface_runoff_file;
    std::string spatial_ice_thickness_file;
    std::string spatial_basal_ice_velocity_file;
    unsigned int result_interval = 1;
    bool write_visualization_results = true;
    bool advect_surface_state = false;
    bool use_sea_level_function = false;
    bool restrict_ocean_to_largest_connected_component = true;
    bool spherical_geometry = false;

    bool true_polar_wander_enabled = false;
    bool include_ice_load_in_true_polar_wander = true;
    bool degree_two_self_gravity_enabled = false;
    bool initialize_self_gravity_in_equilibrium = true;
    double ice_density = 917.0;
    double elastic_degree_two_load_love_number = -0.3;
    double fluid_degree_two_load_love_number = -0.9;
    double self_gravity_relaxation_time = 1e4;
    double rotational_bulge_inertia_difference = 2.6e35;
    double polar_wander_relaxation_time = 1e6;
    double maximum_polar_wander_rate = 10.0;
    Tensor<1,dim> spin_axis;
    Tensor<1,dim> equilibrium_spin_axis;
    SymmetricTensor<2,dim> reference_moment_of_inertia;
    SymmetricTensor<2,dim> delayed_self_gravity_ice_load;
    bool reference_moment_of_inertia_is_initialized = false;
    bool self_gravity_state_is_initialized = false;
    double rigid_ice_load_norm = 0.0;
    double effective_ice_load_norm = 0.0;
    mutable double last_polar_wander_output_time = -1.0;

    bool regional_ice_load_response_enabled = false;
    bool initialize_regional_ice_load_in_equilibrium = true;
    double regional_compensation_density = 3300.0;
    double regional_immediate_response_fraction = 0.0;
    double regional_ice_load_relaxation_time = 1e4;
    mutable std::vector<double> regional_delayed_ice_load_displacement;
    mutable std::vector<double> regional_total_ice_load_displacement;
    mutable std::vector<double> regional_ice_load_velocity;
    mutable bool regional_ice_load_state_is_initialized = false;
};
}
}

#endif
#endif
