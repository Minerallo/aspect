/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#include <aspect/mesh_deformation/fastscape_cpp.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator.h>
#include <aspect/structured_data.h>

#include <fastscapelib/flow/flow_router.hpp>
#include <fastscapelib/flow/sink_resolver.hpp>

#include <deal.II/base/data_out_base.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/data_out.h>

#if __has_include(<xtensor/generators/xbuilder.hpp>)
#  include <xtensor/generators/xbuilder.hpp>
#  include <xtensor/core/xmath.hpp>
#else
#  include <xtensor/xbuilder.hpp>
#  include <xtensor/xmath.hpp>
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>


namespace aspect
{
namespace MeshDeformation
{
/**
 * Read a generic erosion-strength field and interpolate it onto the
 * independent landscape mesh. This class deliberately has no knowledge of
 * the program or physical formula that produced the field.
 */
template <int dim>
class SpatialErosionStrength
{
public:
    void
    initialize(const std::string &filename,
               const std::vector<Point<dim-1>> &surface_coordinates)
    {
        values = xt::ones<double>({surface_coordinates.size()});
        if (filename.empty())
            return;

        lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(1, 1.0);
        lookup->load_file(filename, MPI_COMM_SELF);
        sample(surface_coordinates);
    }

    void
    sample(const std::vector<Point<dim-1>> &surface_coordinates)
    {
        if (!lookup)
            return;

        values.resize({surface_coordinates.size()});
        for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
            values[i] = std::max(0.0,
                                 lookup->get_data(surface_coordinates[i], 0));
    }

    const xt::xarray<double> &
    get_values() const
    {
        return values;
    }

private:
    xt::xarray<double> values;
    std::unique_ptr<Utilities::StructuredDataLookup<dim-1>> lookup;
};

/**
 * Read a dimensionless local surface-runoff field and interpolate it onto
 * the independent landscape mesh. FastScape multiplies each surface-cell
 * area by this value before accumulating water supply downstream. A
 * uniform value of one therefore recovers ordinary drainage area.
 */
template <int dim>
class SpatialSurfaceRunoff
{
public:
    void
    initialize(const std::string &filename,
               const std::vector<Point<dim-1>> &surface_coordinates)
    {
        values = xt::ones<double>({surface_coordinates.size()});
        if (filename.empty())
            return;

        lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(1, 1.0);
        lookup->load_file(filename, MPI_COMM_SELF);
        sample(surface_coordinates);
    }

    void
    sample(const std::vector<Point<dim-1>> &surface_coordinates)
    {
        if (!lookup)
            return;

        values.resize({surface_coordinates.size()});
        for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
            values[i] = std::max(0.0,
                                 lookup->get_data(surface_coordinates[i], 0));
    }

    const xt::xarray<double> &
    get_values() const
    {
        return values;
    }

private:
    xt::xarray<double> values;
    std::unique_ptr<Utilities::StructuredDataLookup<dim-1>> lookup;
};


/**
 * Read ice thickness in meters and interpolate it onto the independent
 * landscape mesh. An omitted file represents an ice-free surface.
 */
template <int dim>
class SpatialIceThickness
{
public:
    void
    initialize(const std::string &filename,
               const std::vector<Point<dim-1>> &surface_coordinates)
    {
        values = xt::zeros<double>({surface_coordinates.size()});
        if (filename.empty())
            return;

        lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(1, 1.0);
        lookup->load_file(filename, MPI_COMM_SELF);
        sample(surface_coordinates);
    }

    void
    sample(const std::vector<Point<dim-1>> &surface_coordinates)
    {
        if (!lookup)
            return;

        values.resize({surface_coordinates.size()});
        for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
            values[i] = std::max(0.0,
                                 lookup->get_data(surface_coordinates[i], 0));
    }

    const xt::xarray<double> &
    get_values() const
    {
        return values;
    }

private:
    xt::xarray<double> values;
    std::unique_ptr<Utilities::StructuredDataLookup<dim-1>> lookup;
};


/**
 * Read basal ice velocity in meters per year and interpolate it onto the
 * independent landscape mesh. An omitted file represents no sliding.
 */
template <int dim>
class SpatialBasalIceVelocity
{
public:
    void
    initialize(const std::string &filename,
               const std::vector<Point<dim-1>> &surface_coordinates)
    {
        values = xt::zeros<double>({surface_coordinates.size()});
        if (filename.empty())
            return;

        lookup = std::make_unique<Utilities::StructuredDataLookup<dim-1>>(1, 1.0);
        lookup->load_file(filename, MPI_COMM_SELF);
        sample(surface_coordinates);
    }

    void
    sample(const std::vector<Point<dim-1>> &surface_coordinates)
    {
        if (!lookup)
            return;

        values.resize({surface_coordinates.size()});
        for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
            values[i] = std::max(0.0,
                                 lookup->get_data(surface_coordinates[i], 0));
    }

    const xt::xarray<double> &
    get_values() const
    {
        return values;
    }

private:
    xt::xarray<double> values;
    std::unique_ptr<Utilities::StructuredDataLookup<dim-1>> lookup;
};


/**
 * Own the FastScape grid, routing graph, erosion solver, and evolving
 * landscape fields. ASPECT-specific field transfer and file output remain
 * outside this class.
 */
template <int dim>
class FastscapeLandscape
{
public:
    using SurfaceMesh = Triangulation<dim-1,dim>;
    using Grid = fastscapelib::dealii_surface_grid<SurfaceMesh>;
    using FlowGraph = fastscapelib::flow_graph<Grid>;
    using Eroder = fastscapelib::spl_eroder<FlowGraph>;
    using SurfaceVelocity = Tensor<1,dim>;

    struct StepResult
    {
        xt::xarray<double> previous_elevation;
        double eroded_volume = 0.0;
        double fluvial_eroded_volume = 0.0;
        double glacial_eroded_volume = 0.0;
        double exported_sediment_flux = 0.0;
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
               const double solver_tolerance,
               const double marine_transport_coefficient,
               const double sediment_porosity,
               const double transport_depth_scale,
               const bool restrict_ocean_connectivity,
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
        restrict_ocean_to_largest_connected_component =
            restrict_ocean_connectivity;
        grid = std::make_unique<Grid>(surface_mesh, closed_surface);
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
        erosion = xt::zeros<double>(flow_graph->grid_shape());
        fluvial_erosion = xt::zeros<double>(flow_graph->grid_shape());
        glacial_erosion = xt::zeros<double>(flow_graph->grid_shape());
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
            const double glacial_erosion_coefficient,
            const double glacial_velocity_exponent,
            const double minimum_ice_thickness,
            const bool advect_surface_state,
            const double maximum_advection_courant,
            const double hillslope_diffusivity,
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
            for (unsigned int i = 0; i < glacial_erosion.size(); ++i)
                glacial_erosion[i] =
                    ice_thickness[i] >= minimum_ice_thickness
                    ? step_years * glacial_erosion_coefficient *
                    exposed_erodibility[i] *
                    std::pow(basal_ice_velocity[i],
                             glacial_velocity_exponent)
                    : 0.0;
            erosion = fluvial_erosion + glacial_erosion;

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
                        local_source_by_lithology[rock][i] += removed;
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
                                   maximum_diffusion_courant);
            if (marine_sediment_transport_coefficient > 0.0)
            {
                const xt::xarray<double> sediment_before_deposition =
                    sediment_thickness;
                for (const auto index : flow_graph->base_levels())
                    if (elevation[index] <= sea_level)
                    {
                        const double solid_volume = sediment_flux[index] * step_years;
                        for (unsigned int rock = 0;
                                rock < lithology_names.size(); ++rock)
                        {
                            const double deposited_thickness =
                                sediment_flux_by_lithology[rock][index] *
                                step_years /
                                ((1.0 - marine_sediment_porosity) * areas[index]);
                            sediment_thickness_by_lithology[rock][index] +=
                                deposited_thickness;
                            sediment_thickness[index] += deposited_thickness;
                        }
                        result.coastal_sediment_flux += solid_volume;
                        result.deposited_sediment_volume +=
                            solid_volume / (1.0 - marine_sediment_porosity);
                    }

                elevation = bedrock_elevation + sediment_thickness;
                transport_marine_sediment(step_years, sea_level);
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
            result.coastal_sediment_flux /= total_time_years;

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

    const xt::xarray<double> &get_elevation() const
    {
        return elevation;
    }

    const xt::xarray<double> &get_drainage_area() const
    {
        return drainage_area;
    }

    const xt::xarray<double> &get_erosion() const
    {
        return erosion;
    }

    const xt::xarray<double> &get_fluvial_erosion() const
    {
        return fluvial_erosion;
    }

    const xt::xarray<double> &get_glacial_erosion() const
    {
        return glacial_erosion;
    }

    const xt::xarray<double> &get_sediment_flux() const
    {
        return sediment_flux;
    }

    const xt::xarray<double> &get_marine_sediment_flux() const
    {
        return marine_sediment_flux;
    }

    const xt::xarray<double> &get_sediment_thickness() const
    {
        return sediment_thickness;
    }

    const std::vector<std::string> &get_lithology_names() const
    {
        return lithology_names;
    }

    const std::vector<unsigned int> &get_bedrock_lithology() const
    {
        return bedrock_lithology;
    }

    const std::vector<xt::xarray<double>> &
    get_sediment_flux_by_lithology() const
    {
        return sediment_flux_by_lithology;
    }

    const std::vector<xt::xarray<double>> &
    get_sediment_thickness_by_lithology() const
    {
        return sediment_thickness_by_lithology;
    }

    std::vector<xt::xarray<double>>
    take_deposited_thickness_by_lithology()
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

    const xt::xarray<double> &get_deposition_rate() const
    {
        return deposition_rate;
    }

    const xt::xarray<double> &get_ocean_mask() const
    {
        return ocean_mask;
    }

    void
    set_elevation(const std::vector<double> &values)
    {
        AssertDimension(values.size(), elevation.size());
        std::copy(values.begin(), values.end(), elevation.begin());
        for (unsigned int i = 0; i < elevation.size(); ++i)
            bedrock_elevation[i] = elevation[i] - sediment_thickness[i];
    }

    void
    set_sediment_thickness(const std::vector<double> &values)
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


    void
    set_lithology_state(
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

private:
    /**
     * Apply conservative linear hillslope transport on the unstructured
     * surface grid. Material removed from a high cell is first taken from
     * mobile sediment and then bedrock; deposition becomes mobile sediment.
     */
    void
    diffuse_hillslopes(const double step_years,
                       const double diffusivity,
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
                const double conductance =
                    diffusivity *
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
                    const double volume =
                        diffusivity *
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

    /**
     * Conservatively advect bedrock elevation and mobile-sediment thickness
     * over the fixed landscape grid using a first-order upwind finite-volume
     * scheme. The velocity is tangential to the ASPECT surface and expressed
     * in meters per year. Courant substeps keep the mobile thickness positive.
     */
    void
    advect_surface_fields(const std::vector<SurfaceVelocity> &velocity,
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
                         advection_step_years, areas);
            sediment_thickness.fill(0.0);
            for (unsigned int rock = 0;
                    rock < lithology_names.size(); ++rock)
            {
                advect_field(sediment_thickness_by_lithology[rock], velocity,
                             advection_step_years, areas);
                for (double &thickness :
                        sediment_thickness_by_lithology[rock])
                    thickness = std::max(0.0, thickness);
                sediment_thickness += sediment_thickness_by_lithology[rock];
            }
        }
        elevation = bedrock_elevation + sediment_thickness;
    }

    void
    advect_field(xt::xarray<double> &field,
                 const std::vector<SurfaceVelocity> &velocity,
                 const double step_years,
                 const typename Grid::container_type &areas)
    {
        std::vector<double> extensive_change(field.size(), 0.0);
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
                const std::size_t donor = signed_rate >= 0.0 ? i : j;
                const std::size_t receiver = signed_rate >= 0.0 ? j : i;
                const double transported =
                    std::abs(signed_rate) * step_years * field[donor];
                extensive_change[donor] -= transported;
                extensive_change[receiver] += transported;
            }

        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] += extensive_change[i] / areas[i];
    }

    /**
     * Move deposited sediment down the seafloor gradient. Each shared face
     * is visited once, and equal volumes are removed from the donor and
     * added to the receiver. A donor-wide limiter prevents transport from
     * removing more sediment than is locally available.
     */
    void
    transport_marine_sediment(const double step_years,
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
        marine_sediment_flux.fill(0.0);

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
    }

    void
    set_base_levels(const xt::xarray<double> &surface_elevation,
                    const double sea_level)
    {
        ocean_mask.fill(0.0);
        std::vector<std::size_t> wet_nodes;
        for (std::size_t i = 0; i < surface_elevation.size(); ++i)
            if (surface_elevation[i] <= sea_level)
                wet_nodes.push_back(i);

        if (!spherical_geometry)
        {
            for (const std::size_t index : wet_nodes)
                ocean_mask[index] = 1.0;
            return;
        }

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
        flow_graph->set_base_levels(base_levels);
    }

    bool spherical_geometry = false;
    double drainage_area_exponent = 0.4;
    double marine_sediment_transport_coefficient = 0.0;
    double marine_sediment_porosity = 0.4;
    double marine_transport_depth_scale = 0.0;
    bool restrict_ocean_to_largest_connected_component = true;
    std::unique_ptr<Grid> grid;
    std::unique_ptr<FlowGraph> flow_graph;
    std::unique_ptr<Eroder> eroder;
    xt::xarray<double> elevation;
    xt::xarray<double> bedrock_elevation;
    xt::xarray<double> drainage_area;
    xt::xarray<double> erosion;
    xt::xarray<double> fluvial_erosion;
    xt::xarray<double> glacial_erosion;
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
template <int dim>
class SurfaceResults
{
public:
    using SurfaceMesh = Triangulation<dim-1,dim>;

    void
    initialize(const xt::xarray<double> &initial_elevation)
    {
        reference_elevation = initial_elevation;
    }

    std::vector<double>
    get_reference_elevation() const
    {
        return std::vector<double>(reference_elevation.begin(),
                                   reference_elevation.end());
    }

    const std::vector<std::pair<double,std::string>> &
    get_output_history() const
    {
        return output_history;
    }

    void
    restore_output_state(
        const std::vector<double> &stored_reference_elevation,
        const std::vector<std::pair<double,std::string>> &stored_output_history)
    {
        AssertDimension(stored_reference_elevation.size(),
                        reference_elevation.size());
        std::copy(stored_reference_elevation.begin(),
                  stored_reference_elevation.end(),
                  reference_elevation.begin());
        output_history = stored_output_history;
    }

    void
    write(const std::string &output_directory,
          const unsigned int timestep_number,
          const double time_years,
          const bool write_visualization,
          const SurfaceMesh &surface_mesh,
          const std::vector<Point<dim>> &surface_points,
          FastscapeLandscape<dim> &landscape,
          const xt::xarray<double> &erosion_strength,
          const xt::xarray<double> &surface_runoff,
          const xt::xarray<double> &ice_thickness,
          const xt::xarray<double> &basal_ice_velocity,
          const double sea_level,
          const double eroded_volume,
          const double fluvial_eroded_volume,
          const double glacial_eroded_volume,
          const double exported_sediment_flux,
          const double coastal_sediment_flux,
          const double deposited_sediment_volume,
          const double stored_sediment_volume)
    {
        const std::string directory =
            output_directory + "fastscape_surface_evolution/";
        std::filesystem::create_directories(directory);

        const auto &elevation = landscape.get_elevation();
        const auto &drainage_area = landscape.get_drainage_area();
        const auto &erosion = landscape.get_erosion();
        const auto &fluvial_erosion = landscape.get_fluvial_erosion();
        const auto &glacial_erosion = landscape.get_glacial_erosion();
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

        const std::string budget_file =
            directory + "sediment_budget.csv";
        const bool write_header = !std::filesystem::exists(budget_file);
        {
            std::ofstream output(budget_file, std::ios::app);
            if (write_header)
                output << "timestep,time_years,sea_level_m,eroded_volume_m3,"
                       << "fluvial_eroded_volume_m3,glacial_eroded_volume_m3,"
                       << "sediment_outflux_m3_per_year,"
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
                   << coastal_sediment_flux << ','
                   << deposited_sediment_volume << ','
                   << stored_sediment_volume << ','
                   << *std::max_element(drainage_area.begin(),
                                        drainage_area.end()) << ','
                   << *std::min_element(elevation.begin(), elevation.end()) << ','
                   << *std::max_element(elevation.begin(), elevation.end()) << '\n';
        }

        const std::string surface_file =
            directory + "surface-" +
            Utilities::int_to_string(timestep_number, 5) + ".csv";
        std::ofstream surface(surface_file);
        surface << "longitude_deg,latitude_deg,elevation_m,erosion_m,"
                << "fluvial_erosion_m,glacial_erosion_m,"
                << "drainage_area_m2,sediment_flux_m3_per_year,"
                << "marine_sediment_flux_m3_per_year,"
                << "sediment_thickness_m,deposition_rate_m_per_year,"
                << "is_connected_ocean,bedrock_lithology,"
                << "erosion_strength,surface_runoff_factor,"
                << "ice_thickness_m,basal_ice_velocity_m_per_year,"
                << "elevation_change_m";
        for (const std::string &name : lithology_names)
            surface << ",sediment_flux_" << name << "_m3_per_year"
                    << ",sediment_thickness_" << name << "_m";
        surface << '\n';
        surface << std::setprecision(16);
        for (unsigned int i = 0; i < elevation.size(); ++i)
        {
            const Point<dim> &point = surface_points[i];
            const double longitude =
                std::atan2(point[1], point[0]) * 180.0 / numbers::PI;
            double latitude = 0.0;
            if constexpr (dim == 3)
                latitude =
                    std::asin(point[2] / point.norm()) * 180.0 / numbers::PI;
            surface << longitude << ',' << latitude << ',' << elevation[i] << ','
                    << erosion[i] << ',' << fluvial_erosion[i] << ','
                    << glacial_erosion[i] << ',' << drainage_area[i] << ','
                    << sediment_flux[i] << ',' << marine_sediment_flux[i] << ','
                    << sediment_thickness[i] << ',' << deposition_rate[i] << ','
                    << ocean_mask[i] << ','
                    << lithology_names[bedrock_lithology[i]] << ','
                    << erosion_strength[i] << ','
                    << surface_runoff[i] << ','
                    << ice_thickness[i] << ','
                    << basal_ice_velocity[i] << ','
                    << elevation[i] - reference_elevation[i];
            for (unsigned int rock = 0; rock < lithology_names.size(); ++rock)
                surface << ',' << sediment_flux_by_lithology[rock][i]
                        << ',' << sediment_thickness_by_lithology[rock][i];
            surface << '\n';
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
            const Point<dim> &point = surface_points[i];
            const double longitude =
                std::atan2(point[1], point[0]) * 180.0 / numbers::PI;
            double latitude = 0.0;
            if constexpr (dim == 3)
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
        Vector<double> flux_output(sediment_flux.size());
        Vector<double> marine_flux_output(marine_sediment_flux.size());
        Vector<double> sediment_thickness_output(sediment_thickness.size());
        Vector<double> deposition_rate_output(deposition_rate.size());
        Vector<double> ocean_mask_output(ocean_mask.size());
        Vector<double> ice_thickness_output(ice_thickness.size());
        Vector<double> basal_ice_velocity_output(basal_ice_velocity.size());
        Vector<double> bedrock_lithology_output(elevation.size());
        for (unsigned int i = 0; i < elevation.size(); ++i)
        {
            elevation_output[i] = elevation[i];
            erosion_output[i] = erosion[i];
            fluvial_erosion_output[i] = fluvial_erosion[i];
            glacial_erosion_output[i] = glacial_erosion[i];
            drainage_output[i] = drainage_area[i];
            flux_output[i] = sediment_flux[i];
            marine_flux_output[i] = marine_sediment_flux[i];
            sediment_thickness_output[i] = sediment_thickness[i];
            deposition_rate_output[i] = deposition_rate[i];
            ocean_mask_output[i] = ocean_mask[i];
            ice_thickness_output[i] = ice_thickness[i];
            basal_ice_velocity_output[i] = basal_ice_velocity[i];
            bedrock_lithology_output[i] = bedrock_lithology[i];
        }

        DataOut<dim-1,dim> data_out;
        data_out.attach_triangulation(surface_mesh);
        data_out.add_data_vector(elevation_output, "elevation",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(erosion_output, "erosion",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(fluvial_erosion_output, "fluvial_erosion",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(glacial_erosion_output, "glacial_erosion",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(drainage_output, "drainage_area",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(flux_output, "sediment_flux",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(marine_flux_output, "marine_sediment_flux",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(sediment_thickness_output, "sediment_thickness",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(deposition_rate_output, "deposition_rate",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(ocean_mask_output, "is_connected_ocean",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(ice_thickness_output, "ice_thickness",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(basal_ice_velocity_output,
                                 "basal_ice_velocity",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(bedrock_lithology_output,
                                 "bedrock_lithology",
                                 DataOut<dim-1,dim>::type_cell_data);
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
                                     DataOut<dim-1,dim>::type_cell_data);
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

private:
    xt::xarray<double> reference_elevation;
    std::vector<std::pair<double,std::string>> output_history;
};


template <int dim>
FastscapeCpp<dim>::FastscapeCpp()
    :
    landscape(std::make_unique<FastscapeLandscape<dim>>()),
    spatial_erosion_strength(
        std::make_unique<SpatialErosionStrength<dim>>()),
    spatial_surface_runoff(
        std::make_unique<SpatialSurfaceRunoff<dim>>()),
    spatial_ice_thickness(
        std::make_unique<SpatialIceThickness<dim>>()),
    spatial_basal_ice_velocity(
        std::make_unique<SpatialBasalIceVelocity<dim>>()),
    surface_results(std::make_unique<SurfaceResults<dim>>())
{
    spin_axis[dim-1] = 1.0;
    equilibrium_spin_axis = spin_axis;
}


template <int dim>
FastscapeCpp<dim>::~FastscapeCpp() = default;


template <int dim>
Point<dim>
FastscapeCpp<dim>::reference_surface_point(const Point<dim> &point) const
{
    if (!spherical_geometry)
        return point;

    const auto *shell =
        dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model());
    Assert(shell != nullptr, ExcInternalError());
    Assert(point.norm() > 0.0, ExcInternalError());
    return point / point.norm() * shell->outer_radius();
}


template <int dim>
Point<dim-1>
FastscapeCpp<dim>::natural_surface_coordinates(const Point<dim> &point) const
{
    const std::array<double,dim> natural =
        this->get_geometry_model().cartesian_to_natural_coordinates(point);
    Point<dim-1> surface_point;
    const unsigned int offset = spherical_geometry ? 1 : 0;
    for (unsigned int d = 0; d < dim-1; ++d)
        surface_point[d] = natural[d+offset];
    return surface_point;
}


template <int dim>
Point<dim-1>
FastscapeCpp<dim>::climate_surface_coordinates(const Point<dim> &point) const
{
    if constexpr (dim != 3)
        return natural_surface_coordinates(point);
    else
    {
        const Tensor<1,3> radial_direction = point / point.norm();
        const Tensor<1,3> climate_north = spin_axis / spin_axis.norm();

        Tensor<1,3> climate_zero_longitude;
        climate_zero_longitude[0] = 1.0;
        climate_zero_longitude -=
            (climate_zero_longitude * climate_north) * climate_north;
        if (climate_zero_longitude.norm() < 1e-12)
        {
            climate_zero_longitude = Tensor<1,3>();
            climate_zero_longitude[1] = 1.0;
            climate_zero_longitude -=
                (climate_zero_longitude * climate_north) * climate_north;
        }
        climate_zero_longitude /= climate_zero_longitude.norm();
        const Tensor<1,3> climate_east =
            cross_product_3d(climate_north, climate_zero_longitude);

        double longitude = std::atan2(radial_direction * climate_east,
                                      radial_direction * climate_zero_longitude);
        if (longitude < 0.0)
            longitude += 2.0 * numbers::PI;
        const double colatitude =
            std::acos(std::clamp(radial_direction * climate_north, -1.0, 1.0));
        return Point<2>(longitude, colatitude);
    }
}


template <int dim>
void
FastscapeCpp<dim>::resample_climate_fields()
{
    if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

    std::vector<Point<dim-1>> climate_coordinates(fastscape_points.size());
    for (unsigned int i = 0; i < fastscape_points.size(); ++i)
        climate_coordinates[i] = climate_surface_coordinates(fastscape_points[i]);

    spatial_erosion_strength->sample(climate_coordinates);
    spatial_surface_runoff->sample(climate_coordinates);
    spatial_ice_thickness->sample(climate_coordinates);
    spatial_basal_ice_velocity->sample(climate_coordinates);
}


template <int dim>
SymmetricTensor<2,dim>
FastscapeCpp<dim>::ice_load_moment_of_inertia() const
{
    SymmetricTensor<2,dim> local_moment;
    if constexpr (dim == 3)
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0 &&
            include_ice_load_in_true_polar_wander)
        {
            const xt::xarray<double> &ice_thickness =
                spatial_ice_thickness->get_values();
            AssertDimension(ice_thickness.size(), fastscape_points.size());
            for (unsigned int i = 0; i < fastscape_points.size(); ++i)
            {
                const Tensor<1,dim> radial_direction =
                    fastscape_points[i] / fastscape_points[i].norm();
                const double surface_mass =
                    ice_density * ice_thickness[i] * fastscape_point_areas[i];
                const double radius_squared = fastscape_points[i].norm_square();
                local_moment += surface_mass * radius_squared *
                    (unit_symmetric_tensor<dim>() -
                     symmetrize(outer_product(radial_direction,
                                              radial_direction)));
            }
        }

    return Utilities::MPI::sum(local_moment, this->get_mpi_communicator());
}


template <int dim>
SymmetricTensor<2,dim>
FastscapeCpp<dim>::apply_degree_two_self_gravity(
    const SymmetricTensor<2,dim> &rigid_ice_load)
{
    rigid_ice_load_norm = rigid_ice_load.norm();
    if (!degree_two_self_gravity_enabled)
    {
        effective_ice_load_norm = rigid_ice_load_norm;
        return rigid_ice_load;
    }

    const SymmetricTensor<2,dim> delayed_equilibrium =
        (fluid_degree_two_load_love_number
         - elastic_degree_two_load_love_number) * rigid_ice_load;
    if (!self_gravity_state_is_initialized)
    {
        delayed_self_gravity_ice_load =
            initialize_self_gravity_in_equilibrium
            ? delayed_equilibrium
            : SymmetricTensor<2,dim>();
        self_gravity_state_is_initialized = true;
    }
    else if (self_gravity_relaxation_time == 0.0)
        delayed_self_gravity_ice_load = delayed_equilibrium;
    else if (this->get_timestep() > 0.0)
    {
        const double time_step_years = this->get_timestep() / year_in_seconds;
        const double relaxed_fraction =
            1.0 - std::exp(-time_step_years / self_gravity_relaxation_time);
        delayed_self_gravity_ice_load +=
            relaxed_fraction
            * (delayed_equilibrium - delayed_self_gravity_ice_load);
    }

    const SymmetricTensor<2,dim> effective_ice_load =
        (1.0 + elastic_degree_two_load_love_number) * rigid_ice_load
        + delayed_self_gravity_ice_load;
    effective_ice_load_norm = effective_ice_load.norm();
    return effective_ice_load;
}


template <int dim>
void
FastscapeCpp<dim>::write_true_polar_wander_state() const
{
    if constexpr (dim == 3)
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0 &&
            this->get_time() != last_polar_wander_output_time)
        {
            const std::string filename =
                this->get_output_directory() + "true_polar_wander.csv";
            const bool write_header = !std::filesystem::exists(filename);
            std::ofstream output(filename, std::ios::app);
            if (write_header)
                output << "time_years,pole_longitude_degrees,pole_latitude_degrees,"
                       << "equilibrium_longitude_degrees,equilibrium_latitude_degrees,"
                       << "rigid_ice_load_kg_m2,effective_ice_load_kg_m2\n";

            const auto longitude = [](const Tensor<1,3> &axis)
            {
                double value = std::atan2(axis[1], axis[0]) * 180.0 / numbers::PI;
                if (value < 0.0)
                    value += 360.0;
                return value;
            };
            const auto latitude = [](const Tensor<1,3> &axis)
            {
                return std::asin(std::clamp(axis[2] / axis.norm(), -1.0, 1.0))
                       * 180.0 / numbers::PI;
            };
            output << std::setprecision(16)
                   << this->get_time() / year_in_seconds << ','
                   << longitude(spin_axis) << ',' << latitude(spin_axis) << ','
                   << longitude(equilibrium_spin_axis) << ','
                   << latitude(equilibrium_spin_axis) << ','
                   << rigid_ice_load_norm << ',' << effective_ice_load_norm << '\n';
            last_polar_wander_output_time = this->get_time();
        }
}


template <int dim>
void
FastscapeCpp<dim>::update_true_polar_wander()
{
    if (!true_polar_wander_enabled)
        return;

    AssertThrow(dim == 3 && spherical_geometry,
                ExcMessage("True polar wander requires a three-dimensional "
                           "spherical-shell model."));

    if constexpr (dim == 3)
    {
        std::vector<double> saved_state(23, 0.0);
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
            for (unsigned int d = 0; d < 3; ++d)
            {
                saved_state[d] = spin_axis[d];
                saved_state[3+d] = equilibrium_spin_axis[d];
            }
            unsigned int entry = 6;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = i; j < 3; ++j)
                    saved_state[entry++] = reference_moment_of_inertia[i][j];
            saved_state[12] = reference_moment_of_inertia_is_initialized ? 1.0 : 0.0;
            saved_state[13] = last_polar_wander_output_time;
            entry = 14;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = i; j < 3; ++j)
                    saved_state[entry++] = delayed_self_gravity_ice_load[i][j];
            saved_state[20] = self_gravity_state_is_initialized ? 1.0 : 0.0;
            saved_state[21] = rigid_ice_load_norm;
            saved_state[22] = effective_ice_load_norm;
        }
        saved_state = Utilities::MPI::broadcast(this->get_mpi_communicator(),
                                                saved_state,
                                                0);
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        {
            for (unsigned int d = 0; d < 3; ++d)
            {
                spin_axis[d] = saved_state[d];
                equilibrium_spin_axis[d] = saved_state[3+d];
            }
            unsigned int entry = 6;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = i; j < 3; ++j)
                    reference_moment_of_inertia[i][j] = saved_state[entry++];
            reference_moment_of_inertia_is_initialized = saved_state[12] > 0.5;
            last_polar_wander_output_time = saved_state[13];
            entry = 14;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = i; j < 3; ++j)
                    delayed_self_gravity_ice_load[i][j] = saved_state[entry++];
            self_gravity_state_is_initialized = saved_state[20] > 0.5;
            rigid_ice_load_norm = saved_state[21];
            effective_ice_load_norm = saved_state[22];
        }

        // Climate fields are inexpensive to reconstruct from their source
        // files. Reconstruct them before evaluating the ice load so a restart
        // uses the restored pole rather than the geographic pole used during
        // surface-mesh initialization.
        resample_climate_fields();

        const RotationProperties<dim> rotation =
            this->compute_net_angular_momentum(false,
                                               this->get_solution(),
                                               false);
        if (!reference_moment_of_inertia_is_initialized)
        {
            reference_moment_of_inertia = rotation.tensor_moment_of_inertia;
            reference_moment_of_inertia_is_initialized = true;
        }

        SymmetricTensor<2,dim> reorientation_moment =
            rotation.tensor_moment_of_inertia - reference_moment_of_inertia;
        reorientation_moment +=
            apply_degree_two_self_gravity(ice_load_moment_of_inertia());
        reorientation_moment += rotational_bulge_inertia_difference *
            symmetrize(outer_product(spin_axis, spin_axis));

        const auto principal_axes =
            eigenvectors(reorientation_moment,
                         SymmetricTensorEigenvectorMethod::jacobi);
        unsigned int maximum_axis = 0;
        for (unsigned int d = 1; d < dim; ++d)
            if (principal_axes[d].first > principal_axes[maximum_axis].first)
                maximum_axis = d;
        equilibrium_spin_axis = principal_axes[maximum_axis].second;
        if (equilibrium_spin_axis * spin_axis < 0.0)
            equilibrium_spin_axis *= -1.0;

        const double cosine =
            std::clamp(equilibrium_spin_axis * spin_axis, -1.0, 1.0);
        const double separation = std::acos(cosine);
        const double time_step_years = this->get_timestep() / year_in_seconds;
        if (separation > 0.0 && time_step_years > 0.0)
        {
            double angular_step = separation;
            if (polar_wander_relaxation_time > 0.0)
                angular_step *=
                    1.0 - std::exp(-time_step_years / polar_wander_relaxation_time);
            if (maximum_polar_wander_rate > 0.0)
                angular_step = std::min(
                    angular_step,
                    maximum_polar_wander_rate * time_step_years / 1e6
                    * numbers::PI / 180.0);

            const double fraction = angular_step / separation;
            if (separation < 1e-10)
                spin_axis = (1.0-fraction) * spin_axis
                            + fraction * equilibrium_spin_axis;
            else
                spin_axis =
                    std::sin((1.0-fraction)*separation) / std::sin(separation)
                    * spin_axis
                    + std::sin(fraction*separation) / std::sin(separation)
                    * equilibrium_spin_axis;
            spin_axis /= spin_axis.norm();
        }

        resample_climate_fields();
        write_true_polar_wander_state();
    }
}


template <int dim>
Tensor<1,dim>
FastscapeCpp<dim>::outward_direction(const Point<dim> &point) const
{
    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(point);
    if (gravity.norm() > 0.0)
        return -gravity / gravity.norm();

    Tensor<1,dim> outward;
    if (spherical_geometry)
        outward = point / point.norm();
    else
        outward[dim-1] = 1.0;
    return outward;
}


template <int dim>
void
FastscapeCpp<dim>::build_surface_mesh()
{
    fastscape_points.clear();
    fastscape_point_areas.clear();

    spherical_geometry =
        (dynamic_cast<const GeometryModel::SphericalShell<dim> *>(
             &this->get_geometry_model()) != nullptr);
    AssertThrow(spherical_geometry ||
                dynamic_cast<const GeometryModel::Box<dim> *>(
                    &this->get_geometry_model()) != nullptr,
                ExcMessage("The FastScape C++ coupling supports only Box and "
                           "SphericalShell geometry models."));
    this->set_surface_transfer_options(surface_transfer_scheme,
                                       surface_transfer_neighbors,
                                       spherical_geometry);

    if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

    if (const auto *box =
                dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model()))
    {
        spherical_geometry = false;
        const Point<dim> origin = box->get_origin();
        const Point<dim> extents = box->get_extents();
        Point<dim-1> lower;
        Point<dim-1> upper;
        for (unsigned int d = 0; d < dim-1; ++d)
        {
            lower[d] = origin[d];
            upper[d] = origin[d] + extents[d];
        }

        const unsigned int refinement_factor = Utilities::pow(2, surface_refinement);
        std::vector<unsigned int> repetitions(dim-1,
                                              box_repetitions * refinement_factor);
        GridGenerator::subdivided_hyper_rectangle(surface_mesh,
                repetitions,
                lower,
                upper,
                true);
        GridTools::transform(
            [&origin, &extents](const Point<dim> &point)
        {
            Point<dim> embedded_point = point;
            embedded_point[dim-1] = origin[dim-1] + extents[dim-1];
            return embedded_point;
        },
        surface_mesh);
    }
    else if (const auto *shell =
                 dynamic_cast<const GeometryModel::SphericalShell<dim> *>(
                     &this->get_geometry_model()))
    {
        spherical_geometry = true;
        GridGenerator::hyper_sphere(surface_mesh, Point<dim>(), shell->outer_radius());
        surface_mesh.refine_global(surface_refinement);
    }
    else
        Assert(false, ExcInternalError());

    for (const auto &cell : surface_mesh.active_cell_iterators())
    {
        fastscape_points.push_back(reference_surface_point(cell->center()));
        fastscape_point_areas.push_back(cell->measure());
    }

    AssertThrow(drainage_area_exponent > 0.0,
                ExcMessage("The drainage area exponent must be positive."));

    xt::xarray<double> initial_elevation =
        xt::zeros<double>({fastscape_points.size()});
    std::vector<Point<dim-1>> surface_coordinates(fastscape_points.size());
    for (unsigned int i = 0; i < fastscape_points.size(); ++i)
    {
        surface_coordinates[i] =
            natural_surface_coordinates(fastscape_points[i]);
        initial_elevation[i] =
            this->get_initial_topography_model().value(surface_coordinates[i]);

        if (initial_relief != 0.0)
        {
            if (spherical_geometry)
            {
                const Point<dim-1> surface =
                    surface_coordinates[i];
                if constexpr (dim == 3)
                    initial_elevation[i] +=
                        initial_relief * std::sin(surface[0]) * std::cos(surface[1]);
                else
                    initial_elevation[i] += initial_relief * std::cos(surface[0]);
            }
            else
            {
                const auto *box =
                    dynamic_cast<const GeometryModel::Box<dim> *>(
                        &this->get_geometry_model());
                Assert(box != nullptr, ExcInternalError());
                double shape_value = 1.0;
                for (unsigned int d = 0; d < dim-1; ++d)
                    shape_value *= std::cos(2.0 * numbers::PI
                                            * (fastscape_points[i][d] - box->get_origin()[d])
                                            / box->get_extents()[d]);
                initial_elevation[i] += initial_relief * shape_value;
            }
        }
    }

    landscape->initialize(surface_mesh,
                          spherical_geometry,
                          initial_elevation,
                          incision_rate,
                          drainage_area_exponent,
                          slope_exponent,
                          nonlinear_tolerance,
                          marine_sediment_transport_coefficient,
                          marine_sediment_porosity,
                          marine_transport_depth_scale,
                          restrict_ocean_to_largest_connected_component,
                          lithology_names,
                          lithology_probabilities,
                          lithology_erodibility_factors,
                          lithology_random_seed);
    spatial_erosion_strength->initialize(spatial_erosion_strength_file,
                                         surface_coordinates);
    spatial_surface_runoff->initialize(spatial_surface_runoff_file,
                                       surface_coordinates);
    spatial_ice_thickness->initialize(spatial_ice_thickness_file,
                                      surface_coordinates);
    spatial_basal_ice_velocity->initialize(spatial_basal_ice_velocity_file,
                                           surface_coordinates);
    surface_results->initialize(initial_elevation);

    // ASPECT applies the initial topography directly to the triangulation,
    // while the transfer framework deliberately uses the undeformed MappingQ
    // rather than the Eulerian mesh-deformation mapping. Evaluation exactly
    // on the reference surface can consequently lie outside valleys in the
    // initial mesh. Sample on a shallow, common subsurface instead. This
    // guarantees point ownership and is also more robust than evaluating
    // exactly on a partition boundary.
    const double inward_offset =
        1e-6 * this->get_geometry_model().length_scale();
    const double maximum_landscape_height =
        std::max(std::abs(*std::min_element(initial_elevation.begin(),
                                            initial_elevation.end())),
                 std::abs(*std::max_element(initial_elevation.begin(),
                                            initial_elevation.end())));
    const double sampling_depth =
        2.0 * maximum_landscape_height + inward_offset;
    for (unsigned int i = 0; i < fastscape_points.size(); ++i)
        fastscape_points[i] -=
            sampling_depth * outward_direction(fastscape_points[i]);
    this->get_pcout() << "   FastScape surface grid: "
                      << fastscape_points.size() << " nodes ("
                      << (spherical_geometry ? "global spherical" : "box")
                      << ")." << std::endl;
}


template <int dim>
void
FastscapeCpp<dim>::initialize()
{
    build_surface_mesh();
}


template <int dim>
void
FastscapeCpp<dim>::update()
{
    if (!this->remote_point_evaluator)
    {
        this->set_evaluation_point_areas(fastscape_point_areas);
        this->set_evaluation_points(fastscape_points);
    }

    if (use_sea_level_function)
        sea_level_function.set_time(this->get_time() / year_in_seconds);

    update_true_polar_wander();
}


template <int dim>
std::vector<Tensor<1,dim>>
                        FastscapeCpp<dim>::compute_updated_velocities_at_points(
                            const std::vector<std::vector<double>> &solution_at_points) const
{
    AssertDimension(solution_at_points.size(), this->evaluation_points.size());
    std::vector<Tensor<1,dim>> result(solution_at_points.size());

    if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0 ||
            solution_at_points.empty() ||
            this->get_timestep() <= 0.0)
        return result;

    Assert(landscape && spatial_erosion_strength &&
           spatial_surface_runoff && spatial_ice_thickness &&
           spatial_basal_ice_velocity && surface_results,
           ExcInternalError());
    AssertDimension(solution_at_points.size(),
                    landscape->get_elevation().size());

    const double aspect_dt_years = this->get_timestep() / year_in_seconds;
    const double current_sea_level =
        use_sea_level_function
        ? sea_level_function.value(Point<1>())
        : sea_level;
    xt::xarray<double> uplift_rate =
        xt::zeros<double>(landscape->get_elevation().shape());
    std::vector<Tensor<1,dim>> tangential_velocity(
        solution_at_points.size());
    for (unsigned int i = 0; i < solution_at_points.size(); ++i)
    {
        Tensor<1,dim> material_velocity;
        for (unsigned int d = 0; d < dim; ++d)
            material_velocity[d] =
                solution_at_points[i][this->introspection().component_indices.velocities[d]];
        const Tensor<1,dim> surface_normal =
            outward_direction(this->evaluation_points[i]);
        const double normal_material_velocity =
            material_velocity * surface_normal;
        uplift_rate[i] = normal_material_velocity * year_in_seconds;
        tangential_velocity[i] =
            (material_velocity -
             normal_material_velocity * surface_normal) * year_in_seconds;
    }

    const typename FastscapeLandscape<dim>::StepResult landscape_step =
        landscape->advance(
            uplift_rate,
            tangential_velocity,
            aspect_dt_years,
            landscape_steps_per_geodynamic_step,
            maximum_landscape_step_years,
            current_sea_level,
            spatial_erosion_strength->get_values(),
            spatial_surface_runoff->get_values(),
            spatial_ice_thickness->get_values(),
            spatial_basal_ice_velocity->get_values(),
            glacial_erosion_coefficient,
            glacial_velocity_exponent,
            minimum_ice_thickness,
            advect_surface_state,
            maximum_surface_advection_courant,
            hillslope_diffusion_coefficient,
            maximum_hillslope_diffusion_courant);

    for (unsigned int i = 0; i < result.size(); ++i)
    {
        const double normal_velocity =
            (landscape->get_elevation()[i]
             - landscape_step.previous_elevation[i])
            / this->get_timestep();
        result[i] = normal_velocity * outward_direction(this->evaluation_points[i]);
    }

    if (result_interval > 0 &&
            this->get_timestep_number() % result_interval == 0)
        surface_results->write(
            this->get_output_directory(),
            this->get_timestep_number(),
            this->get_time() / year_in_seconds,
            write_visualization_results,
            surface_mesh,
            fastscape_points,
            *landscape,
            spatial_erosion_strength->get_values(),
            spatial_surface_runoff->get_values(),
            spatial_ice_thickness->get_values(),
            spatial_basal_ice_velocity->get_values(),
            current_sea_level,
            landscape_step.eroded_volume,
            landscape_step.fluvial_eroded_volume,
            landscape_step.glacial_eroded_volume,
            landscape_step.exported_sediment_flux,
            landscape_step.coastal_sediment_flux,
            landscape_step.deposited_sediment_volume,
            landscape_step.stored_sediment_volume);

    return result;
}


template <int dim>
bool
FastscapeCpp<dim>::needs_surface_stabilization() const
{
    return true;
}


template <int dim>
void
FastscapeCpp<dim>::save(
    std::map<std::string, std::string> &status_strings) const
{
    if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

    std::vector<double> elevation_values(landscape->get_elevation().begin(),
                                         landscape->get_elevation().end());
    std::vector<double> sediment_thickness_values(
        landscape->get_sediment_thickness().begin(),
        landscape->get_sediment_thickness().end());
    std::ostringstream stream;
    {
        aspect::oarchive archive(stream);
        archive << elevation_values;
    }
    status_strings["FastscapeSurfaceEvolution"] = stream.str();

    std::ostringstream sediment_stream;
    {
        aspect::oarchive archive(sediment_stream);
        archive << sediment_thickness_values;
    }
    status_strings["FastscapeMarineSediment"] = sediment_stream.str();

    std::vector<std::vector<double>> sediment_by_lithology;
    for (const auto &field : landscape->get_sediment_thickness_by_lithology())
        sediment_by_lithology.emplace_back(field.begin(), field.end());
    std::ostringstream lithology_stream;
    {
        aspect::oarchive archive(lithology_stream);
        archive << landscape->get_lithology_names();
        archive << landscape->get_bedrock_lithology();
        archive << sediment_by_lithology;
    }
    status_strings["FastscapeLithologyProvenance"] =
        lithology_stream.str();

    std::ostringstream output_state_stream;
    {
        aspect::oarchive archive(output_state_stream);
        archive << surface_results->get_reference_elevation();
        archive << surface_results->get_output_history();
    }
    status_strings["FastscapeSurfaceOutputState"] =
        output_state_stream.str();

    std::ostringstream polar_wander_stream;
    {
        aspect::oarchive archive(polar_wander_stream);
        archive << spin_axis;
        archive << equilibrium_spin_axis;
        archive << reference_moment_of_inertia;
        archive << delayed_self_gravity_ice_load;
        archive << reference_moment_of_inertia_is_initialized;
        archive << self_gravity_state_is_initialized;
        archive << rigid_ice_load_norm;
        archive << effective_ice_load_norm;
        archive << last_polar_wander_output_time;
    }
    status_strings["FastscapeTruePolarWanderDegreeTwoSelfGravity"] =
        polar_wander_stream.str();
}


template <int dim>
void
FastscapeCpp<dim>::load(
    const std::map<std::string, std::string> &status_strings)
{
    if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

    const auto state =
        status_strings.find("FastscapeSurfaceEvolution");
    AssertThrow(state != status_strings.end(),
                ExcMessage("No FastScape surface-evolution state was found "
                           "in the checkpoint."));

    std::vector<double> elevation_values;
    std::istringstream stream(state->second);
    aspect::iarchive archive(stream);
    archive >> elevation_values;
    landscape->set_elevation(elevation_values);

    const auto sediment_state =
        status_strings.find("FastscapeMarineSediment");
    if (sediment_state != status_strings.end())
    {
        std::vector<double> sediment_thickness_values;
        std::istringstream sediment_stream(sediment_state->second);
        aspect::iarchive sediment_archive(sediment_stream);
        sediment_archive >> sediment_thickness_values;
        landscape->set_sediment_thickness(sediment_thickness_values);
    }

    const auto lithology_state =
        status_strings.find("FastscapeLithologyProvenance");
    if (lithology_state != status_strings.end())
    {
        std::vector<std::string> stored_names;
        std::vector<unsigned int> stored_bedrock_lithology;
        std::vector<std::vector<double>> stored_sediment_thickness;
        std::istringstream lithology_stream(lithology_state->second);
        aspect::iarchive lithology_archive(lithology_stream);
        lithology_archive >> stored_names;
        lithology_archive >> stored_bedrock_lithology;
        lithology_archive >> stored_sediment_thickness;
        AssertThrow(stored_names == landscape->get_lithology_names(),
                    ExcMessage("The lithology names in the checkpoint differ "
                               "from the current parameter file."));
        landscape->set_lithology_state(stored_bedrock_lithology,
                                       stored_sediment_thickness);
    }

    const auto output_state =
        status_strings.find("FastscapeSurfaceOutputState");
    if (output_state != status_strings.end())
    {
        std::vector<double> reference_elevation;
        std::vector<std::pair<double,std::string>> output_history;
        std::istringstream output_state_stream(output_state->second);
        aspect::iarchive output_state_archive(output_state_stream);
        output_state_archive >> reference_elevation;
        output_state_archive >> output_history;
        surface_results->restore_output_state(reference_elevation,
                                              output_history);
    }

    const auto polar_wander_state = status_strings.find(
        "FastscapeTruePolarWanderDegreeTwoSelfGravity");
    if (polar_wander_state != status_strings.end())
    {
        std::istringstream polar_wander_stream(polar_wander_state->second);
        aspect::iarchive archive(polar_wander_stream);
        archive >> spin_axis;
        archive >> equilibrium_spin_axis;
        archive >> reference_moment_of_inertia;
        archive >> delayed_self_gravity_ice_load;
        archive >> reference_moment_of_inertia_is_initialized;
        archive >> self_gravity_state_is_initialized;
        archive >> rigid_ice_load_norm;
        archive >> effective_ice_load_norm;
        archive >> last_polar_wander_output_time;
    }
    else
    {
        // Checkpoints written before the degree-two load response was added
        // contain only the pole and reference inertia state.
        const auto legacy_polar_wander_state =
            status_strings.find("FastscapeTruePolarWander");
        if (legacy_polar_wander_state != status_strings.end())
        {
            std::istringstream polar_wander_stream(
                legacy_polar_wander_state->second);
            aspect::iarchive archive(polar_wander_stream);
            archive >> spin_axis;
            archive >> equilibrium_spin_axis;
            archive >> reference_moment_of_inertia;
            archive >> reference_moment_of_inertia_is_initialized;
            archive >> last_polar_wander_output_time;
        }
    }
}


template <int dim>
void
FastscapeCpp<dim>::declare_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Mesh deformation");
    prm.enter_subsection("FastScape surface evolution");
    {
        prm.declare_entry("Box repetitions", "8", Patterns::Integer(1),
                          "Number of FastScape cells in each horizontal box "
                          "direction before surface refinement.");
        prm.declare_entry("Surface refinement level", "2", Patterns::Integer(0),
                          "Global refinement of the independent surface grid.");
        prm.declare_entry("Surface transfer scheme", "conservative",
                          Patterns::Selection("nearest|weighted|conservative"),
                          "Method used to transfer FastScape surface motion "
                          "back to ASPECT. 'nearest' copies the closest cell, "
                          "'weighted' uses an inverse-distance stencil, and "
                          "'conservative' additionally preserves the separate "
                          "area-integrated positive and negative surface motions.");
        prm.declare_entry("Surface transfer neighbors", "8",
                          Patterns::Integer(1),
                          "Number of nearby FastScape cells used by the "
                          "weighted and conservative transfer methods.");
        prm.declare_entry("Landscape steps per geodynamic step", "4",
                          Patterns::Integer(1),
                          "Minimum number of landscape-evolution steps during "
                          "one geodynamic step.");
        prm.declare_entry("Maximum landscape step", "10000",
                          Patterns::Double(0),
                          "Maximum landscape-evolution step in years.");
        prm.declare_entry("River incision coefficient", "5e-5",
                          Patterns::Double(0),
                          "Coefficient controlling river incision.");
        prm.declare_entry("Drainage area exponent", "0.4", Patterns::Double(0),
                          "Exponent controlling the influence of upstream "
                          "drainage area on river incision.");
        prm.declare_entry("Slope exponent", "1.0", Patterns::Double(0),
                          "Exponent controlling the influence of surface slope "
                          "on river incision.");
        prm.declare_entry("Nonlinear tolerance", "1e-5", Patterns::Double(0),
                          "Convergence tolerance used by the nonlinear river "
                          "incision solver.");
        prm.declare_entry("Glacial erosion coefficient", "0",
                          Patterns::Double(0),
                          "Coefficient in the glacial erosion law E = K u^m, "
                          "where E is erosion rate in meters per year, u is "
                          "basal ice velocity in meters per year, and m is the "
                          "glacial velocity exponent. Zero disables glacial "
                          "erosion.");
        prm.declare_entry("Glacial velocity exponent", "1",
                          Patterns::Double(0),
                          "Exponent m applied to basal ice velocity in the "
                          "glacial erosion law.");
        prm.declare_entry("Minimum ice thickness", "1",
                          Patterns::Double(0),
                          "Minimum ice thickness in meters required for "
                          "glacial erosion.");
        prm.declare_entry("Initial relief", "0", Patterns::Double(0),
                          "Amplitude in meters of a deterministic initial "
                          "FastScape relief field. This is independent of "
                          "ASPECT's geometry-level initial topography and is "
                          "recommended for global spherical models.");
        prm.declare_entry("Sea level", "0", Patterns::Double(),
                          "For a global closed surface, nodes at or below this "
                          "elevation are drainage base levels.");
        prm.declare_entry("Use sea level function", "false",
                          Patterns::Bool(),
                          "Use a time-dependent sea-level function instead of "
                          "the constant sea-level value.");
        prm.enter_subsection("Sea level function");
        {
            Functions::ParsedFunction<1>::declare_parameters(prm, 1);
        }
        prm.leave_subsection();
        prm.declare_entry("Marine sediment transport coefficient", "0",
                          Patterns::Double(0),
                          "Diffusive transport coefficient in square meters "
                          "per year for deposited ocean sediment. Zero keeps "
                          "the earlier behavior in which river sediment leaves "
                          "the landscape at drainage base levels.");
        prm.declare_entry("Marine sediment porosity", "0.4",
                          Patterns::Double(0, 0.999999),
                          "Pore-space fraction used to convert solid sediment "
                          "delivered by rivers into deposited bulk thickness.");
        prm.declare_entry("Marine transport depth scale", "0",
                          Patterns::Double(0),
                          "Water-depth scale in meters over which marine "
                          "transport decreases exponentially. Zero uses a "
                          "depth-independent transport coefficient.");
        prm.declare_entry("Restrict ocean to largest connected water body",
                          "true",
                          Patterns::Bool(),
                          "Treat only the largest face-connected group of "
                          "below-sea-level cells as the global ocean. This "
                          "prevents disconnected inland depressions from "
                          "receiving marine sediment.");
        prm.declare_entry("Advect surface state", "false",
                          Patterns::Bool(),
                          "Conservatively advect bedrock elevation and mobile "
                          "sediment over the fixed FastScape grid using the "
                          "tangential ASPECT material velocity.");
        prm.declare_entry("Maximum surface advection Courant number", "0.5",
                          Patterns::Double(0),
                          "Maximum finite-volume Courant number used while "
                          "advecting FastScape surface state.");
        prm.declare_entry("Hillslope diffusion coefficient", "0",
                          Patterns::Double(0),
                          "Conservative linear hillslope diffusivity in "
                          "square meters per year on the unstructured "
                          "FastScape surface grid.");
        prm.declare_entry("Maximum hillslope diffusion Courant number", "0.25",
                          Patterns::Double(0),
                          "Maximum explicit Courant number used for "
                          "unstructured hillslope diffusion.");
        prm.declare_entry("Lithology names", "upper_crust",
                          Patterns::List(Patterns::Anything(), 1),
                          "Comma-separated names of source-rock classes. The "
                          "names are used in provenance and stratigraphy "
                          "output, for example 'granite, limestone'.");
        prm.declare_entry("Lithology probabilities", "1",
                          Patterns::List(Patterns::Double(0), 1),
                          "Relative probabilities used once to make a fixed, "
                          "reproducible bedrock-lithology map. Supply one "
                          "nonnegative value for every lithology name.");
        prm.declare_entry("Lithology erodibility factors", "1",
                          Patterns::List(Patterns::Double(0), 1),
                          "Relative multiplier of river-incision strength for "
                          "each lithology. Values are relative, so one can be "
                          "used as the reference rock.");
        prm.declare_entry("Lithology random seed", "1",
                          Patterns::Integer(0),
                          "Seed for the reproducible initial lithology map. "
                          "Rock identity is not redrawn during evolution.");
        prm.declare_entry("Spatial erosion strength file", "",
                          Patterns::Anything(),
                          "Optional ASPECT structured text-data file containing "
                          "one nonnegative erosion-strength field. An empty "
                          "filename uses a uniform value of one.");
        prm.declare_entry("Spatial surface runoff file", "",
                          Patterns::Anything(),
                          "Optional ASPECT structured text-data file containing "
                          "dimensionless local surface runoff. FastScape "
                          "multiplies cell area by this value before accumulating "
                          "water supply downstream. An empty filename uses a "
                          "uniform value of one.");
        prm.declare_entry("Spatial ice thickness file", "",
                          Patterns::Anything(),
                          "Optional ASPECT structured text-data file containing "
                          "ice thickness in meters. An empty filename represents "
                          "an ice-free surface.");
        prm.declare_entry("Spatial basal ice velocity file", "",
                          Patterns::Anything(),
                          "Optional ASPECT structured text-data file containing "
                          "basal ice velocity in meters per year. An empty "
                          "filename represents no basal sliding.");
        prm.declare_entry("Enable true polar wander", "false",
                          Patterns::Bool(),
                          "Let changes in ASPECT's directly integrated moment "
                          "of inertia and the prescribed ice load move the spin "
                          "axis. Climate fields are then sampled in coordinates "
                          "defined by that axis, so rain and ice belts move over "
                          "the body-fixed surface. This option is available only "
                          "for three-dimensional spherical shells.");
        prm.declare_entry("Include ice load in true polar wander", "true",
                          Patterns::Bool(),
                          "Include the moment of inertia of the prescribed ice "
                          "thickness. Topography evolved by FastScape already "
                          "changes ASPECT's volume integral and is not added a "
                          "second time as a separate surface load.");
        prm.declare_entry("Enable degree two self gravity", "false",
                          Patterns::Bool(),
                          "Correct the rigid ice-load inertia tensor for the "
                          "self-gravitating deformation of the solid Earth. "
                          "The correction uses an immediate elastic degree-two "
                          "load Love number and one delayed viscous mode. This "
                          "is a reduced load-response model, not a self-gravity "
                          "body force in the ASPECT Stokes equations.");
        prm.declare_entry("Elastic degree two load Love number", "-0.3",
                          Patterns::Double(-1.0, 0.0),
                          "Immediate elastic degree-two gravity-potential "
                          "response to a surface load. The effective initial "
                          "load tensor is multiplied by one plus this value.");
        prm.declare_entry("Fluid degree two load Love number", "-0.9",
                          Patterns::Double(-1.0, 0.0),
                          "Long-time degree-two response after viscous "
                          "relaxation. A value approaching minus one represents "
                          "nearly complete compensation of the surface load.");
        prm.declare_entry("Self gravity relaxation time", "10000",
                          Patterns::Double(0),
                          "Relaxation time in years of the single delayed "
                          "degree-two load-response mode. Zero applies the "
                          "fluid response immediately.");
        prm.declare_entry("Initialize self gravity in equilibrium", "true",
                          Patterns::Bool(),
                          "Assume the initial ice load has already reached its "
                          "long-time compensated state. Disable this for an "
                          "ice load emplaced at the beginning of the model.");
        prm.declare_entry("Ice density", "917",
                          Patterns::Double(0),
                          "Density in kilograms per cubic meter used to convert "
                          "ice thickness into surface mass.");
        prm.declare_entry("Rotational bulge inertia difference", "2.6e35",
                          Patterns::Double(0),
                          "Difference between polar and equatorial moments of "
                          "inertia in kilograms square meters. This stabilizing "
                          "term represents the hydrostatic or fossil rotational "
                          "bulge. It is essential because ASPECT does not solve "
                          "self-gravitating deformation of that bulge.");
        prm.declare_entry("Polar wander relaxation time", "1e6",
                          Patterns::Double(0),
                          "Time in years over which the spin axis approaches "
                          "the maximum principal inertia axis. Zero applies the "
                          "change immediately, subject to the rate limit.");
        prm.declare_entry("Maximum polar wander rate", "10",
                          Patterns::Double(0),
                          "Maximum spin-axis motion in degrees per million "
                          "years. Zero disables this limit.");
        prm.declare_entry("Result interval", "1", Patterns::Integer(0),
                          "Number of geodynamic steps between landscape result "
                          "files. Zero disables result files.");
        prm.declare_entry("Write visualization results", "true",
                          Patterns::Bool(),
                          "Write cell-based elevation, erosion, drainage-area, "
                          "and sediment-flux visualization files.");
    }
    prm.leave_subsection();
    prm.leave_subsection();
}


template <int dim>
void
FastscapeCpp<dim>::parse_parameters(ParameterHandler &prm)
{
    prm.enter_subsection("Mesh deformation");
    prm.enter_subsection("FastScape surface evolution");
    {
        box_repetitions = prm.get_integer("Box repetitions");
        surface_refinement = prm.get_integer("Surface refinement level");
        surface_transfer_scheme = prm.get("Surface transfer scheme");
        surface_transfer_neighbors =
            prm.get_integer("Surface transfer neighbors");
        landscape_steps_per_geodynamic_step =
            prm.get_integer("Landscape steps per geodynamic step");
        maximum_landscape_step_years =
            prm.get_double("Maximum landscape step");
        incision_rate = prm.get_double("River incision coefficient");
        drainage_area_exponent = prm.get_double("Drainage area exponent");
        slope_exponent = prm.get_double("Slope exponent");
        nonlinear_tolerance = prm.get_double("Nonlinear tolerance");
        glacial_erosion_coefficient =
            prm.get_double("Glacial erosion coefficient");
        glacial_velocity_exponent =
            prm.get_double("Glacial velocity exponent");
        AssertThrow(glacial_velocity_exponent > 0.0,
                    ExcMessage("Glacial velocity exponent must be positive."));
        minimum_ice_thickness =
            prm.get_double("Minimum ice thickness");
        initial_relief = prm.get_double("Initial relief");
        sea_level = prm.get_double("Sea level");
        use_sea_level_function = prm.get_bool("Use sea level function");
        if (use_sea_level_function)
        {
            prm.enter_subsection("Sea level function");
            {
                sea_level_function.parse_parameters(prm);
            }
            prm.leave_subsection();
        }
        marine_sediment_transport_coefficient =
            prm.get_double("Marine sediment transport coefficient");
        marine_sediment_porosity =
            prm.get_double("Marine sediment porosity");
        marine_transport_depth_scale =
            prm.get_double("Marine transport depth scale");
        restrict_ocean_to_largest_connected_component =
            prm.get_bool("Restrict ocean to largest connected water body");
        advect_surface_state =
            prm.get_bool("Advect surface state");
        maximum_surface_advection_courant =
            prm.get_double("Maximum surface advection Courant number");
        hillslope_diffusion_coefficient =
            prm.get_double("Hillslope diffusion coefficient");
        maximum_hillslope_diffusion_courant =
            prm.get_double("Maximum hillslope diffusion Courant number");
        lithology_names =
            Utilities::split_string_list(prm.get("Lithology names"));
        lithology_probabilities = Utilities::string_to_double(
            Utilities::split_string_list(prm.get("Lithology probabilities")));
        lithology_erodibility_factors = Utilities::string_to_double(
            Utilities::split_string_list(
                prm.get("Lithology erodibility factors")));
        lithology_random_seed = prm.get_integer("Lithology random seed");
        AssertThrow(lithology_names.size() == lithology_probabilities.size() &&
                    lithology_names.size() ==
                    lithology_erodibility_factors.size(),
                    ExcMessage("Lithology names, probabilities, and "
                               "erodibility factors must contain the same "
                               "number of entries."));
        AssertThrow(std::accumulate(lithology_probabilities.begin(),
                                    lithology_probabilities.end(), 0.0) > 0.0,
                    ExcMessage("At least one lithology probability must be "
                               "positive."));
        for (const std::string &name : lithology_names)
            AssertThrow(!name.empty() &&
                        std::all_of(name.begin(), name.end(),
            [](const unsigned char character)
        {
            return std::isalnum(character) || character == '_';
        }),
        ExcMessage("Lithology names may contain only letters, numbers, and "
                   "underscores because they are used as output field names."));
        spatial_erosion_strength_file =
            prm.get("Spatial erosion strength file");
        spatial_surface_runoff_file =
            prm.get("Spatial surface runoff file");
        spatial_ice_thickness_file =
            prm.get("Spatial ice thickness file");
        spatial_basal_ice_velocity_file =
            prm.get("Spatial basal ice velocity file");
        true_polar_wander_enabled =
            prm.get_bool("Enable true polar wander");
        include_ice_load_in_true_polar_wander =
            prm.get_bool("Include ice load in true polar wander");
        degree_two_self_gravity_enabled =
            prm.get_bool("Enable degree two self gravity");
        elastic_degree_two_load_love_number =
            prm.get_double("Elastic degree two load Love number");
        fluid_degree_two_load_love_number =
            prm.get_double("Fluid degree two load Love number");
        AssertThrow(fluid_degree_two_load_love_number <=
                    elastic_degree_two_load_love_number,
                    ExcMessage("The fluid degree-two load Love number must be "
                               "less than or equal to the elastic value."));
        self_gravity_relaxation_time =
            prm.get_double("Self gravity relaxation time");
        initialize_self_gravity_in_equilibrium =
            prm.get_bool("Initialize self gravity in equilibrium");
        ice_density = prm.get_double("Ice density");
        rotational_bulge_inertia_difference =
            prm.get_double("Rotational bulge inertia difference");
        polar_wander_relaxation_time =
            prm.get_double("Polar wander relaxation time");
        maximum_polar_wander_rate =
            prm.get_double("Maximum polar wander rate");
        result_interval = prm.get_integer("Result interval");
        write_visualization_results =
            prm.get_bool("Write visualization results");
    }
    prm.leave_subsection();
    prm.leave_subsection();
}
}
}


namespace aspect
{
namespace MeshDeformation
{
ASPECT_REGISTER_MESH_DEFORMATION_MODEL(
    FastscapeCpp,
    "fastscape surface evolution",
    "Uses the FastScape library to evolve an independent surface mesh. "
    "The model transfers ASPECT material velocity to the landscape, computes "
    "river incision, sediment routing, and optional conservative marine "
    "sediment deposition, and returns surface-normal velocity to ASPECT. "
    "It supports box and global spherical-shell geometries.")
}
}

#endif
