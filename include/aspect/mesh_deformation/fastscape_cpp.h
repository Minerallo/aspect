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
    Tensor<1,dim> outward_direction(const Point<dim> &point) const;

    mutable SurfaceMesh surface_mesh;
    std::vector<Point<dim>> fastscape_points;
    mutable std::unique_ptr<FastscapeLandscape<dim>> landscape;
    std::unique_ptr<SpatialErosionStrength<dim>> spatial_erosion_strength;
    std::unique_ptr<SpatialSurfaceRunoff<dim>> spatial_surface_runoff;
    mutable std::unique_ptr<SurfaceResults<dim>> surface_results;

    unsigned int box_repetitions = 8;
    unsigned int surface_refinement = 2;
    unsigned int landscape_steps_per_geodynamic_step = 4;
    double maximum_landscape_step_years = 10000.0;
    double incision_rate = 5e-5;
    double drainage_area_exponent = 0.4;
    double slope_exponent = 1.0;
    double nonlinear_tolerance = 1e-5;
    double initial_relief = 0.0;
    double sea_level = 0.0;
    Functions::ParsedFunction<1> sea_level_function;
    double marine_sediment_transport_coefficient = 0.0;
    double marine_sediment_porosity = 0.4;
    double marine_transport_depth_scale = 0.0;
    double maximum_surface_advection_courant = 0.5;
    double hillslope_diffusion_coefficient = 0.0;
    double maximum_hillslope_diffusion_courant = 0.25;
    std::string spatial_erosion_strength_file;
    std::string spatial_surface_runoff_file;
    unsigned int result_interval = 1;
    bool write_visualization_results = true;
    bool advect_surface_state = false;
    bool use_sea_level_function = false;
    bool restrict_ocean_to_largest_connected_component = true;
    bool spherical_geometry = false;
};
}
}

#endif
#endif
