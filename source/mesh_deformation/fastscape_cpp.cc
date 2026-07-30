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
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
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
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

        Utilities::StructuredDataLookup<dim-1> lookup(1, 1.0);
        lookup.load_file(filename, MPI_COMM_SELF);
        for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
            values[i] = std::max(0.0,
                                 lookup.get_data(surface_coordinates[i], 0));
    }

    const xt::xarray<double> &
    get_values() const
    {
        return values;
    }

private:
    xt::xarray<double> values;
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

        Utilities::StructuredDataLookup<dim-1> lookup(1, 1.0);
        lookup.load_file(filename, MPI_COMM_SELF);
        for (unsigned int i = 0; i < surface_coordinates.size(); ++i)
            values[i] = std::max(0.0,
                                 lookup.get_data(surface_coordinates[i], 0));
    }

    const xt::xarray<double> &
    get_values() const
    {
        return values;
    }

private:
    xt::xarray<double> values;
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

    struct StepResult
    {
        xt::xarray<double> previous_elevation;
        double eroded_volume = 0.0;
        double exported_sediment_flux = 0.0;
    };

    void
    initialize(SurfaceMesh &surface_mesh,
               const bool closed_surface,
               const xt::xarray<double> &initial_elevation,
               const double river_incision_coefficient,
               const double area_exponent,
               const double surface_slope_exponent,
               const double solver_tolerance)
    {
        spherical_geometry = closed_surface;
        drainage_area_exponent = area_exponent;
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
        drainage_area = xt::zeros<double>(flow_graph->grid_shape());
        erosion = xt::zeros<double>(flow_graph->grid_shape());
        sediment_flux = xt::zeros<double>(flow_graph->grid_shape());
    }

    StepResult
    advance(const xt::xarray<double> &uplift_rate,
            const double total_time_years,
            unsigned int number_of_steps,
            const double maximum_step_years,
            const double sea_level,
            const xt::xarray<double> &erosion_strength,
            const xt::xarray<double> &surface_runoff)
    {
        Assert(grid && flow_graph && eroder, ExcInternalError());
        AssertDimension(uplift_rate.size(), elevation.size());
        AssertDimension(erosion_strength.size(), elevation.size());
        AssertDimension(surface_runoff.size(), elevation.size());

        number_of_steps = std::max(1u, number_of_steps);
        while (total_time_years / number_of_steps > maximum_step_years)
            number_of_steps *= 2;
        const double step_years = total_time_years / number_of_steps;

        StepResult result;
        result.previous_elevation = elevation;

        for (unsigned int step = 0; step < number_of_steps; ++step)
        {
            const xt::xarray<double> uplifted =
                elevation + step_years * uplift_rate;
            set_base_levels(uplifted, sea_level);
            flow_graph->update_routes(uplifted);
            flow_graph->accumulate(drainage_area, surface_runoff);

            // Erosion strength is a local multiplier. Surface runoff is
            // handled separately above so that water supplied upstream is
            // carried through the drainage network.
            xt::xarray<double> effective_drainage_area = drainage_area;
            for (unsigned int i = 0; i < effective_drainage_area.size(); ++i)
                effective_drainage_area[i] *=
                    std::pow(erosion_strength[i],
                             1.0 / drainage_area_exponent);

            erosion =
                eroder->erode(uplifted, effective_drainage_area, step_years);
            sediment_flux = flow_graph->accumulate(erosion / step_years);
            elevation = uplifted - erosion;

            const auto areas = grid->nodes_areas();
            for (unsigned int i = 0; i < erosion.size(); ++i)
                result.eroded_volume += erosion[i] * areas[i];
        }

        for (const auto index : flow_graph->base_levels())
            result.exported_sediment_flux += sediment_flux[index];

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

    const xt::xarray<double> &get_sediment_flux() const
    {
        return sediment_flux;
    }

    void
    set_elevation(const std::vector<double> &values)
    {
        AssertDimension(values.size(), elevation.size());
        std::copy(values.begin(), values.end(), elevation.begin());
    }

private:
    void
    set_base_levels(const xt::xarray<double> &surface_elevation,
                    const double sea_level)
    {
        if (!spherical_geometry)
            return;

        std::vector<std::size_t> base_levels;
        for (std::size_t i = 0; i < surface_elevation.size(); ++i)
            if (surface_elevation[i] <= sea_level)
                base_levels.push_back(i);

        if (base_levels.empty())
            base_levels.push_back(static_cast<std::size_t>(
                                      std::distance(
                                          surface_elevation.begin(),
                                          std::min_element(surface_elevation.begin(),
                                                  surface_elevation.end()))));
        flow_graph->set_base_levels(base_levels);
    }

    bool spherical_geometry = false;
    double drainage_area_exponent = 0.4;
    std::unique_ptr<Grid> grid;
    std::unique_ptr<FlowGraph> flow_graph;
    std::unique_ptr<Eroder> eroder;
    xt::xarray<double> elevation;
    xt::xarray<double> drainage_area;
    xt::xarray<double> erosion;
    xt::xarray<double> sediment_flux;
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

    void
    reset_reference_elevation(const xt::xarray<double> &elevation)
    {
        reference_elevation = elevation;
    }

    void
    write(const std::string &output_directory,
          const unsigned int timestep_number,
          const double time_years,
          const bool write_visualization,
          const SurfaceMesh &surface_mesh,
          const std::vector<Point<dim>> &surface_points,
          const FastscapeLandscape<dim> &landscape,
          const xt::xarray<double> &erosion_strength,
          const xt::xarray<double> &surface_runoff,
          const double eroded_volume,
          const double exported_sediment_flux)
    {
        const std::string directory =
            output_directory + "fastscape_surface_evolution/";
        std::filesystem::create_directories(directory);

        const auto &elevation = landscape.get_elevation();
        const auto &drainage_area = landscape.get_drainage_area();
        const auto &erosion = landscape.get_erosion();
        const auto &sediment_flux = landscape.get_sediment_flux();

        const std::string budget_file =
            directory + "sediment_budget.csv";
        const bool write_header = !std::filesystem::exists(budget_file);
        {
            std::ofstream output(budget_file, std::ios::app);
            if (write_header)
                output << "timestep,time_years,eroded_volume_m3,"
                       << "sediment_outflux_m3_per_year,"
                       << "max_drainage_area_m2,"
                       << "min_elevation_m,max_elevation_m\n";
            output << timestep_number << ','
                   << std::setprecision(16) << time_years << ','
                   << eroded_volume << ','
                   << exported_sediment_flux << ','
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
                << "drainage_area_m2,sediment_flux_m3_per_year,"
                << "erosion_strength,surface_runoff_factor,"
                << "elevation_change_m\n";
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
                    << erosion[i] << ',' << drainage_area[i] << ','
                    << sediment_flux[i] << ','
                    << erosion_strength[i] << ','
                    << surface_runoff[i] << ','
                    << elevation[i] - reference_elevation[i] << '\n';
        }

        if (!write_visualization)
            return;

        Vector<double> elevation_output(elevation.size());
        Vector<double> erosion_output(erosion.size());
        Vector<double> drainage_output(drainage_area.size());
        Vector<double> flux_output(sediment_flux.size());
        for (unsigned int i = 0; i < elevation.size(); ++i)
        {
            elevation_output[i] = elevation[i];
            erosion_output[i] = erosion[i];
            drainage_output[i] = drainage_area[i];
            flux_output[i] = sediment_flux[i];
        }

        DataOut<dim-1,dim> data_out;
        data_out.attach_triangulation(surface_mesh);
        data_out.add_data_vector(elevation_output, "elevation",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(erosion_output, "erosion",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(drainage_output, "drainage_area",
                                 DataOut<dim-1,dim>::type_cell_data);
        data_out.add_data_vector(flux_output, "sediment_flux",
                                 DataOut<dim-1,dim>::type_cell_data);
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
    surface_results(std::make_unique<SurfaceResults<dim>>())
{}


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

    spherical_geometry =
        (dynamic_cast<const GeometryModel::SphericalShell<dim> *>(
             &this->get_geometry_model()) != nullptr);
    AssertThrow(spherical_geometry ||
                dynamic_cast<const GeometryModel::Box<dim> *>(
                    &this->get_geometry_model()) != nullptr,
                ExcMessage("The FastScape C++ coupling supports only Box and "
                           "SphericalShell geometry models."));
    this->set_normalized_surface_transfer(spherical_geometry);

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
        fastscape_points.push_back(reference_surface_point(cell->center()));

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
                          nonlinear_tolerance);
    spatial_erosion_strength->initialize(spatial_erosion_strength_file,
                                         surface_coordinates);
    spatial_surface_runoff->initialize(spatial_surface_runoff_file,
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
        this->set_evaluation_points(fastscape_points);
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
           spatial_surface_runoff && surface_results,
           ExcInternalError());
    AssertDimension(solution_at_points.size(),
                    landscape->get_elevation().size());

    const double aspect_dt_years = this->get_timestep() / year_in_seconds;
    xt::xarray<double> uplift_rate =
        xt::zeros<double>(landscape->get_elevation().shape());
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
    }

    const typename FastscapeLandscape<dim>::StepResult landscape_step =
        landscape->advance(
            uplift_rate,
            aspect_dt_years,
            landscape_steps_per_geodynamic_step,
            maximum_landscape_step_years,
            sea_level,
            spatial_erosion_strength->get_values(),
            spatial_surface_runoff->get_values());

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
            landscape_step.eroded_volume,
            landscape_step.exported_sediment_flux);

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
    std::ostringstream stream;
    {
        aspect::oarchive archive(stream);
        archive << elevation_values;
    }
    status_strings["FastscapeSurfaceEvolution"] = stream.str();
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
    surface_results->reset_reference_elevation(landscape->get_elevation());
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
        prm.declare_entry("Initial relief", "0", Patterns::Double(0),
                          "Amplitude in meters of a deterministic initial "
                          "FastScape relief field. This is independent of "
                          "ASPECT's geometry-level initial topography and is "
                          "recommended for global spherical models.");
        prm.declare_entry("Sea level", "0", Patterns::Double(),
                          "For a global closed surface, nodes at or below this "
                          "elevation are drainage base levels.");
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
        landscape_steps_per_geodynamic_step =
            prm.get_integer("Landscape steps per geodynamic step");
        maximum_landscape_step_years =
            prm.get_double("Maximum landscape step");
        incision_rate = prm.get_double("River incision coefficient");
        drainage_area_exponent = prm.get_double("Drainage area exponent");
        slope_exponent = prm.get_double("Slope exponent");
        nonlinear_tolerance = prm.get_double("Nonlinear tolerance");
        initial_relief = prm.get_double("Initial relief");
        sea_level = prm.get_double("Sea level");
        spatial_erosion_strength_file =
            prm.get("Spatial erosion strength file");
        spatial_surface_runoff_file =
            prm.get("Spatial surface runoff file");
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
    "river incision and sediment routing, and returns surface-normal velocity "
    "to ASPECT. "
    "It supports box and global spherical-shell geometries.")
}
}

#endif
