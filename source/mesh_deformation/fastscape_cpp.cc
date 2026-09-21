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

#include <aspect/mesh_deformation/fastscape_cpp.h>
#include <aspect/mesh_deformation/fastscape_cpp_helpers.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>


namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    FastscapeCpp<dim>::FastscapeCpp()
      :
      landscape(std::make_unique<FastscapeLandscape<dim-1,dim>>()),
      spatial_erosion_strength(
        std::make_unique<SpatialErosionStrength<dim-1>>()),
      spatial_surface_runoff(
        std::make_unique<SpatialSurfaceRunoff<dim-1>>()),
      spatial_ice_thickness(
        std::make_unique<SpatialIceThickness<dim-1>>()),
      spatial_basal_ice_velocity(
        std::make_unique<SpatialBasalIceVelocity<dim-1>>()),
      surface_results(std::make_unique<SurfaceResults<dim-1,dim>>()),
      t_landscape(std::make_unique<FastscapeLandscape<2,3>>()),
      t_spatial_erosion_strength(std::make_unique<SpatialErosionStrength<2>>()),
      t_spatial_surface_runoff(std::make_unique<SpatialSurfaceRunoff<2>>()),
      t_spatial_ice_thickness(std::make_unique<SpatialIceThickness<2>>()),
      t_spatial_basal_ice_velocity(std::make_unique<SpatialBasalIceVelocity<2>>()),
      t_surface_results(std::make_unique<SurfaceResults<2,3>>())
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
    std::vector<double>
    FastscapeCpp<dim>::update_regional_ice_load_response(
      const double time_step_years,
      const xt::xarray<double> &ice_thickness) const
    {
      std::vector<double> velocity(ice_thickness.size(), 0.0);
      if (!regional_ice_load_response_enabled)
        {
          regional_ice_load_velocity = velocity;
          return velocity;
        }
      if (time_step_years <= 0.0)
        return velocity;

      AssertThrow(!spherical_geometry,
                  ExcMessage("The regional ice-load response is intended for "
                             "box geometries. Use the degree-two self-gravity "
                             "option for a global spherical model."));
      AssertThrow(regional_compensation_density > 0.0,
                  ExcMessage("Regional compensation density must be positive."));

      if (!regional_ice_load_state_is_initialized)
        {
          regional_delayed_ice_load_displacement.assign(ice_thickness.size(), 0.0);
          regional_total_ice_load_displacement.assign(ice_thickness.size(), 0.0);
          regional_ice_load_velocity.assign(ice_thickness.size(), 0.0);
          if (initialize_regional_ice_load_in_equilibrium)
            for (unsigned int i = 0; i < ice_thickness.size(); ++i)
              {
                const double equilibrium_displacement =
                  -ice_density * ice_thickness[i] /
                  regional_compensation_density;
                regional_delayed_ice_load_displacement[i] =
                  (1.0 - regional_immediate_response_fraction) *
                  equilibrium_displacement;
                regional_total_ice_load_displacement[i] =
                  equilibrium_displacement;
              }
          regional_ice_load_state_is_initialized = true;
        }

      AssertDimension(regional_delayed_ice_load_displacement.size(),
                      ice_thickness.size());
      AssertDimension(regional_total_ice_load_displacement.size(),
                      ice_thickness.size());
      const double relaxed_fraction =
        regional_ice_load_relaxation_time == 0.0
        ? 1.0
        : 1.0 - std::exp(-time_step_years /
                         regional_ice_load_relaxation_time);
      for (unsigned int i = 0; i < ice_thickness.size(); ++i)
        {
          const double equilibrium_displacement =
            -ice_density * ice_thickness[i] /
            regional_compensation_density;
          const double delayed_equilibrium =
            (1.0 - regional_immediate_response_fraction) *
            equilibrium_displacement;
          regional_delayed_ice_load_displacement[i] +=
            relaxed_fraction *
            (delayed_equilibrium -
             regional_delayed_ice_load_displacement[i]);
          const double new_total_displacement =
            regional_immediate_response_fraction * equilibrium_displacement +
            regional_delayed_ice_load_displacement[i];
          velocity[i] =
            (new_total_displacement -
             regional_total_ice_load_displacement[i]) /
            time_step_years;
          regional_total_ice_load_displacement[i] = new_total_displacement;
        }
      regional_ice_load_velocity = velocity;
      return velocity;
    }


    template <int dim>
    std::pair<xt::xarray<double>,xt::xarray<double>>
    FastscapeCpp<dim>::effective_ice_fields(
      const xt::xarray<double> &surface_elevation,
      const xt::xarray<double> &spatial_ice_thickness_values,
      const xt::xarray<double> &spatial_basal_ice_velocity_values) const
    {
      AssertDimension(surface_elevation.size(),
                      spatial_ice_thickness_values.size());
      AssertDimension(surface_elevation.size(),
                      spatial_basal_ice_velocity_values.size());
      xt::xarray<double> thickness = spatial_ice_thickness_values;
      xt::xarray<double> velocity = spatial_basal_ice_velocity_values;

      if (ice_distribution_mode == "spatial field")
        return {thickness, velocity};

      for (unsigned int i = 0; i < surface_elevation.size(); ++i)
        {
          const double excess =
            surface_elevation[i] - glacial_elevation_threshold;
          const double factor = glacial_elevation_transition_width > 0.0
                                ? std::max(0.0, std::min(1.0,
                                                         excess / glacial_elevation_transition_width))
                                : (excess >= 0.0 ? 1.0 : 0.0);
          if (ice_distribution_mode == "elevation threshold")
            {
              thickness[i] = factor * elevation_threshold_ice_thickness;
              velocity[i] = factor * elevation_threshold_basal_ice_velocity;
            }
          else
            {
              thickness[i] *= factor;
              velocity[i] *= factor;
            }
        }
      return {thickness, velocity};
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
      const auto *box_geometry =
        dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());
      const auto *two_merged_boxes_geometry =
        dynamic_cast<const GeometryModel::TwoMergedBoxes<dim> *>(
          &this->get_geometry_model());
      const auto *chunk_geometry =
        dynamic_cast<const GeometryModel::Chunk<dim> *>(
          &this->get_geometry_model());
      AssertThrow(spherical_geometry ||
                  box_geometry != nullptr ||
                  two_merged_boxes_geometry != nullptr ||
                  chunk_geometry != nullptr,
                  ExcMessage("The FastScape C++ coupling supports only Box, "
                             "Box with lithosphere boundary indicators, and "
                             "Chunk and SphericalShell geometry models."));

      if (use_t_coupling)
        {
          build_t_coupling_surface_mesh();
          return;
        }
      this->set_surface_transfer_options(surface_transfer_scheme,
                                         surface_transfer_neighbors,
                                         spherical_geometry);

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

      if (box_geometry != nullptr || two_merged_boxes_geometry != nullptr)
        {
          spherical_geometry = false;
          const Point<dim> origin = box_geometry != nullptr
                                    ? box_geometry->get_origin()
                                    : two_merged_boxes_geometry->get_origin();
          const Point<dim> extents = box_geometry != nullptr
                                     ? box_geometry->get_extents()
                                     : two_merged_boxes_geometry->get_extents();
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

          std::vector<GridTools::PeriodicFacePair<typename SurfaceMesh::cell_iterator>>
          periodicity;
          for (const unsigned int direction : periodic_surface_dimensions)
            GridTools::collect_periodic_faces(surface_mesh,
                                              2*direction,
                                              2*direction+1,
                                              direction,
                                              periodicity);
          if (!periodicity.empty())
            surface_mesh.add_periodicity(periodicity);

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
                  Assert(box_geometry != nullptr ||
                         two_merged_boxes_geometry != nullptr,
                         ExcInternalError());
                  const Point<dim> origin = box_geometry != nullptr
                                            ? box_geometry->get_origin()
                                            : two_merged_boxes_geometry->get_origin();
                  const Point<dim> extents = box_geometry != nullptr
                                             ? box_geometry->get_extents()
                                             : two_merged_boxes_geometry->get_extents();
                  double shape_value = 1.0;
                  for (unsigned int d = 0; d < dim-1; ++d)
                    shape_value *= std::cos(2.0 * numbers::PI
                                            * (fastscape_points[i][d] - origin[d])
                                            / extents[d]);
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
                            flow_routing_method,
                            flow_routing_connectivity,
                            multiple_flow_slope_exponent,
                            nonlinear_tolerance,
                            marine_sediment_transport_coefficient,
                            marine_sediment_porosity,
                            marine_transport_depth_scale,
                            limit_marine_deposition_to_available_accommodation,
                            maximum_marine_deposition_above_sea_level,
                            use_sea_level_as_drainage_base_level,
                            restrict_ocean_to_largest_connected_component,
                            submarine_river_incision_factor,
                            open_marine_sediment_boundary,
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
      regional_delayed_ice_load_displacement.assign(initial_elevation.size(), 0.0);
      regional_total_ice_load_displacement.assign(initial_elevation.size(), 0.0);
      regional_ice_load_velocity.assign(initial_elevation.size(), 0.0);
      regional_ice_load_state_is_initialized = false;

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
    FastscapeCpp<dim>::build_t_coupling_surface_mesh()
    {
      AssertThrow(dim == 2,
                  ExcMessage("FastScape T coupling is only meaningful for a "
                             "two-dimensional ASPECT model."));
      AssertThrow(t_coupling_width > 0.0,
                  ExcMessage("FastScape T coupling width must be positive."));
      if constexpr (dim == 2)
        {
          const auto *box = dynamic_cast<const GeometryModel::Box<2> *>(
                              &this->get_geometry_model());
          const auto *two_merged_boxes =
            dynamic_cast<const GeometryModel::TwoMergedBoxes<2> *>(
              &this->get_geometry_model());
          const auto *chunk = dynamic_cast<const GeometryModel::Chunk<2> *>(
                                &this->get_geometry_model());
          const auto *shell =
            dynamic_cast<const GeometryModel::SphericalShell<2> *>(
              &this->get_geometry_model());
          Assert(box != nullptr || two_merged_boxes != nullptr ||
                 chunk != nullptr || shell != nullptr,
                 ExcInternalError());
          const bool curved_stem = chunk != nullptr || shell != nullptr;
          spherical_geometry = curved_stem;
          this->set_surface_transfer_options(surface_transfer_scheme,
                                             surface_transfer_neighbors,
                                             curved_stem);
          if (Utilities::MPI::this_mpi_process(
                this->get_mpi_communicator()) != 0)
            return;

          Point<2> origin;
          Point<2> extents;
          double outer_radius = 0.0;
          double minimum_angle = 0.0;
          double angle_range = 0.0;
          if (box != nullptr || two_merged_boxes != nullptr)
            {
              origin = box != nullptr ? box->get_origin()
                       : two_merged_boxes->get_origin();
              extents = box != nullptr ? box->get_extents()
                        : two_merged_boxes->get_extents();
            }
          else if (chunk != nullptr)
            {
              outer_radius = chunk->outer_radius();
              minimum_angle = chunk->west_longitude();
              angle_range = chunk->longitude_range();
              extents[0] = outer_radius * angle_range;
            }
          else
            {
              outer_radius = shell->outer_radius();
              angle_range = shell->opening_angle()
                            * numbers::PI / 180.0;
              extents[0] = outer_radius * angle_range;
            }
          const unsigned int refinement_factor =
            Utilities::pow(2, surface_refinement);
          const unsigned int nx = t_coupling_cell_size > 0.0
                                  ? std::max(1u, static_cast<unsigned int>(
                                               std::ceil(extents[0] / t_coupling_cell_size)))
                                  : box_repetitions * refinement_factor;
          const double dx = extents[0] / nx;
          const unsigned int ny = std::max(
                                    1u, static_cast<unsigned int>(std::lround(t_coupling_width / dx)));
          const double represented_width = ny * dx;

          const Point<2> lower(origin[0], -0.5 * represented_width);
          const Point<2> upper(origin[0] + extents[0],
                               0.5 * represented_width);
          Triangulation<2,2> planar_surface_mesh;
          GridGenerator::subdivided_hyper_rectangle(
            planar_surface_mesh, std::vector<unsigned int> {nx, ny},
            lower, upper, true);
          GridGenerator::flatten_triangulation(planar_surface_mesh,
                                               t_surface_mesh);

          std::vector<GridTools::PeriodicFacePair<
          typename Triangulation<2,3>::cell_iterator>> periodicity;
          for (const unsigned int direction : periodic_surface_dimensions)
            GridTools::collect_periodic_faces(t_surface_mesh,
                                              2*direction,
                                              2*direction+1,
                                              direction,
                                              periodicity);
          if (!periodicity.empty())
            t_surface_mesh.add_periodicity(periodicity);

          fastscape_points.resize(nx);
          fastscape_point_areas.assign(nx, dx);
          for (unsigned int column = 0; column < nx; ++column)
            {
              const double longitudinal_distance = (column + 0.5) * dx;
              if (!curved_stem)
                fastscape_points[column] =
                  Point<2>(origin[0] + longitudinal_distance,
                           origin[1] + extents[1]);
              else
                {
                  const std::array<double,2> natural =
                  {
                    {
                      outer_radius,
                      minimum_angle + longitudinal_distance / outer_radius
                    }
                  };
                  fastscape_points[column] =
                    this->get_geometry_model().natural_to_cartesian_coordinates(
                      natural);
                }
            }

          t_surface_points.clear();
          t_surface_point_areas.clear();
          t_column_of_cell.clear();
          t_cells_by_column.assign(nx, {});
          std::vector<Point<2>> climate_coordinates;
          xt::xarray<double> initial_elevation =
            xt::zeros<double>({t_surface_mesh.n_active_cells()});
          unsigned int index = 0;
          for (const auto &cell : t_surface_mesh.active_cell_iterators())
            {
              const Point<3> point = cell->center();
              t_surface_points.push_back(point);
              t_surface_point_areas.push_back(cell->measure());
              climate_coordinates.emplace_back(point[0], point[1]);
              const unsigned int column = std::min(
                                            nx - 1,
                                            static_cast<unsigned int>(
                                              std::floor((point[0] - origin[0]) / dx)));
              t_column_of_cell.push_back(column);
              t_cells_by_column[column].push_back(index);

              const Point<1> surface_coordinate = curved_stem
                                                  ? Point<1>(minimum_angle +
                                                             (point[0] - origin[0]) / outer_radius)
                                                  : Point<1>(point[0]);
              initial_elevation[index] =
                this->get_initial_topography_model().value(surface_coordinate);
              if (initial_relief != 0.0)
                initial_elevation[index] += initial_relief
                                            * std::cos(2.0 * numbers::PI
                                                       * (point[0] - origin[0]) / extents[0])
                                            * std::cos(2.0 * numbers::PI
                                                       * (point[1] - lower[1]) / represented_width);
              ++index;
            }

          t_landscape->initialize(t_surface_mesh,
                                  false,
                                  initial_elevation,
                                  incision_rate,
                                  drainage_area_exponent,
                                  slope_exponent,
                                  flow_routing_method,
                                  flow_routing_connectivity,
                                  multiple_flow_slope_exponent,
                                  nonlinear_tolerance,
                                  marine_sediment_transport_coefficient,
                                  marine_sediment_porosity,
                                  marine_transport_depth_scale,
                                  limit_marine_deposition_to_available_accommodation,
                                  maximum_marine_deposition_above_sea_level,
                                  use_sea_level_as_drainage_base_level,
                                  restrict_ocean_to_largest_connected_component,
                                  submarine_river_incision_factor,
                                  open_marine_sediment_boundary,
                                  lithology_names,
                                  lithology_probabilities,
                                  lithology_erodibility_factors,
                                  lithology_random_seed);
          t_spatial_erosion_strength->initialize(spatial_erosion_strength_file,
                                                 climate_coordinates);
          t_spatial_surface_runoff->initialize(spatial_surface_runoff_file,
                                               climate_coordinates);
          t_spatial_ice_thickness->initialize(spatial_ice_thickness_file,
                                              climate_coordinates);
          t_spatial_basal_ice_velocity->initialize(
            spatial_basal_ice_velocity_file, climate_coordinates);
          t_surface_results->initialize(initial_elevation);

          regional_delayed_ice_load_displacement.assign(
            initial_elevation.size(), 0.0);
          regional_total_ice_load_displacement.assign(
            initial_elevation.size(), 0.0);
          regional_ice_load_velocity.assign(initial_elevation.size(), 0.0);
          regional_ice_load_state_is_initialized = false;

          const double inward_offset =
            1e-6 * this->get_geometry_model().length_scale();
          const double maximum_landscape_height =
            std::max(std::abs(*std::min_element(initial_elevation.begin(),
                                                initial_elevation.end())),
                     std::abs(*std::max_element(initial_elevation.begin(),
                                                initial_elevation.end())));
          const double sampling_depth =
            2.0 * maximum_landscape_height + inward_offset;
          for (Point<2> &point : fastscape_points)
            point[1] -= sampling_depth;

          this->get_pcout()
              << "   FastScape T-coupled surface grid: " << nx << " x " << ny
              << " cells, represented width " << represented_width << " m, "
              << t_surface_points.size() << " landscape cells ("
              << (curved_stem ? "curved" : "Cartesian")
              << " ASPECT stem)." << std::endl;
        }
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
    bool
    FastscapeCpp<dim>::should_write_surface_results() const
    {
      if (output_interval > 0.0)
        {
          if (last_output_time <
              this->get_parameters().start_time - output_interval)
            last_output_time = this->get_parameters().start_time;

          if (this->get_time() >= last_output_time + output_interval ||
              this->get_time() == this->get_end_time())
            {
              const double magic =
                1.0 + 2.0 * std::numeric_limits<double>::epsilon();
              last_output_time +=
                std::floor((this->get_time() - last_output_time) /
                           output_interval * magic) *
                output_interval / magic;
              return true;
            }
          return false;
        }

      return result_interval > 0 &&
             this->get_timestep_number() % result_interval == 0;
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

      const double aspect_dt_years = this->get_timestep() / year_in_seconds;
      const double current_sea_level =
        use_sea_level_function
        ? sea_level_function.value(Point<1>())
        : sea_level;

      if (use_t_coupling)
        {
          Assert(dim == 2, ExcInternalError());
          Assert(t_landscape && t_spatial_erosion_strength &&
                 t_spatial_surface_runoff && t_spatial_ice_thickness &&
                 t_spatial_basal_ice_velocity && t_surface_results,
                 ExcInternalError());
          AssertDimension(solution_at_points.size(), t_cells_by_column.size());

          xt::xarray<double> uplift_rate =
            xt::zeros<double>(t_landscape->get_elevation().shape());
          std::vector<Tensor<1,3>> tangential_velocity(
            t_landscape->get_elevation().size());
          for (unsigned int cell = 0; cell < t_column_of_cell.size(); ++cell)
            {
              const unsigned int column = t_column_of_cell[cell];
              Tensor<1,dim> material_velocity;
              for (unsigned int d = 0; d < dim; ++d)
                material_velocity[d] = solution_at_points[column]
                                       [this->introspection().component_indices.velocities[d]];
              const Tensor<1,dim> normal =
                outward_direction(this->evaluation_points[column]);
              const double normal_material_velocity =
                material_velocity * normal;
              uplift_rate[cell] = apply_normal_material_velocity
                                  ? normal_material_velocity * year_in_seconds
                                  : 0.0;
              tangential_velocity[cell][0] =
                (material_velocity - normal_material_velocity * normal)[0]
                * year_in_seconds;
            }

          const auto effective_ice = effective_ice_fields(
                                       t_landscape->get_elevation(),
                                       t_spatial_ice_thickness->get_values(),
                                       t_spatial_basal_ice_velocity->get_values());
          const std::vector<double> regional_load_velocity =
            update_regional_ice_load_response(
              aspect_dt_years,
              effective_ice.first);
          AssertDimension(regional_load_velocity.size(), uplift_rate.size());
          for (unsigned int i = 0; i < uplift_rate.size(); ++i)
            uplift_rate[i] += regional_load_velocity[i];

          const typename FastscapeLandscape<2,3>::StepResult landscape_step =
            t_landscape->advance(
              uplift_rate,
              tangential_velocity,
              aspect_dt_years,
              landscape_steps_per_geodynamic_step,
              maximum_landscape_step_years,
              current_sea_level,
              t_spatial_erosion_strength->get_values(),
              t_spatial_surface_runoff->get_values(),
              effective_ice.first,
              effective_ice.second,
              glacial_erosion_mode,
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
              routed_glacier_spread_erosion,
              routed_glacier_terminate_at_sea_level,
              glacial_erosion_coefficient,
              glacial_velocity_exponent,
              glacial_ice_thickness_scale,
              minimum_ice_thickness,
              restrict_glacial_erosion_to_grounded_ice,
              minimum_glacial_erosion_elevation,
              ice_density,
              seawater_density,
              advect_surface_state,
              maximum_surface_advection_courant,
              hillslope_diffusion_coefficient,
              submarine_hillslope_diffusion_coefficient >= 0.0
              ? submarine_hillslope_diffusion_coefficient
              : hillslope_diffusion_coefficient,
              maximum_marine_sediment_transport_courant,
              maximum_hillslope_diffusion_courant);

          for (unsigned int column = 0; column < result.size(); ++column)
            {
              Assert(!t_cells_by_column[column].empty(), ExcInternalError());
              double elevation_change = 0.0;
              if (t_coupling_reduction == "center")
                {
                  const auto center = *std::min_element(
                                        t_cells_by_column[column].begin(),
                                        t_cells_by_column[column].end(),
                                        [this](const unsigned int a, const unsigned int b)
                  {
                    return std::abs(t_surface_points[a][1]) <
                           std::abs(t_surface_points[b][1]);
                  });
                  elevation_change = t_landscape->get_elevation()[center]
                                     - landscape_step.previous_elevation[center];
                }
              else
                {
                  double total_area = 0.0;
                  for (const unsigned int cell : t_cells_by_column[column])
                    {
                      elevation_change += t_surface_point_areas[cell]
                                          * (t_landscape->get_elevation()[cell]
                                             - landscape_step.previous_elevation[cell]);
                      total_area += t_surface_point_areas[cell];
                    }
                  elevation_change /= total_area;
                }
              result[column] = elevation_change / this->get_timestep()
                               * outward_direction(this->evaluation_points[column]);
            }

          t_surface_results->write_budget(
            this->get_output_directory(),
            this->get_timestep_number(),
            this->get_time() / year_in_seconds,
            current_sea_level,
            landscape_step.eroded_volume,
            landscape_step.fluvial_eroded_volume,
            landscape_step.glacial_eroded_volume,
            landscape_step.exported_sediment_flux,
            landscape_step.accommodation_limited_exported_sediment_flux,
            landscape_step.coastal_sediment_flux,
            landscape_step.deposited_sediment_volume,
            landscape_step.stored_sediment_volume,
            *t_landscape);

          if (should_write_surface_results())
            t_surface_results->write(
              this->get_output_directory(),
              this->get_timestep_number(),
              this->get_time() / year_in_seconds,
              write_visualization_results,
              output_drainage_diagnostics,
              t_surface_mesh,
              t_surface_points,
              *t_landscape,
              t_spatial_erosion_strength->get_values(),
              t_spatial_surface_runoff->get_values(),
              effective_ice.first,
              effective_ice.second,
              regional_total_ice_load_displacement,
              regional_ice_load_velocity);
          return result;
        }

      Assert(landscape && spatial_erosion_strength &&
             spatial_surface_runoff && spatial_ice_thickness &&
             spatial_basal_ice_velocity && surface_results,
             ExcInternalError());
      AssertDimension(solution_at_points.size(),
                      landscape->get_elevation().size());

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
          uplift_rate[i] = apply_normal_material_velocity
                           ? normal_material_velocity * year_in_seconds
                           : 0.0;
          tangential_velocity[i] =
            (material_velocity -
             normal_material_velocity * surface_normal) * year_in_seconds;
        }

      const auto effective_ice = effective_ice_fields(
                                   landscape->get_elevation(),
                                   spatial_ice_thickness->get_values(),
                                   spatial_basal_ice_velocity->get_values());
      const std::vector<double> regional_load_velocity =
        update_regional_ice_load_response(
          aspect_dt_years, effective_ice.first);
      AssertDimension(regional_load_velocity.size(), uplift_rate.size());
      for (unsigned int i = 0; i < uplift_rate.size(); ++i)
        uplift_rate[i] += regional_load_velocity[i];

      const typename FastscapeLandscape<dim-1,dim>::StepResult landscape_step =
        landscape->advance(
          uplift_rate,
          tangential_velocity,
          aspect_dt_years,
          landscape_steps_per_geodynamic_step,
          maximum_landscape_step_years,
          current_sea_level,
          spatial_erosion_strength->get_values(),
          spatial_surface_runoff->get_values(),
          effective_ice.first,
          effective_ice.second,
          glacial_erosion_mode,
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
          routed_glacier_spread_erosion,
          routed_glacier_terminate_at_sea_level,
          glacial_erosion_coefficient,
          glacial_velocity_exponent,
          glacial_ice_thickness_scale,
          minimum_ice_thickness,
          restrict_glacial_erosion_to_grounded_ice,
          minimum_glacial_erosion_elevation,
          ice_density,
          seawater_density,
          advect_surface_state,
          maximum_surface_advection_courant,
          hillslope_diffusion_coefficient,
          submarine_hillslope_diffusion_coefficient >= 0.0
          ? submarine_hillslope_diffusion_coefficient
          : hillslope_diffusion_coefficient,
          maximum_marine_sediment_transport_courant,
          maximum_hillslope_diffusion_courant);

      for (unsigned int i = 0; i < result.size(); ++i)
        {
          const double normal_velocity =
            (landscape->get_elevation()[i]
             - landscape_step.previous_elevation[i])
            / this->get_timestep();
          result[i] = normal_velocity * outward_direction(this->evaluation_points[i]);
        }

      surface_results->write_budget(
        this->get_output_directory(),
        this->get_timestep_number(),
        this->get_time() / year_in_seconds,
        current_sea_level,
        landscape_step.eroded_volume,
        landscape_step.fluvial_eroded_volume,
        landscape_step.glacial_eroded_volume,
        landscape_step.exported_sediment_flux,
        landscape_step.accommodation_limited_exported_sediment_flux,
        landscape_step.coastal_sediment_flux,
        landscape_step.deposited_sediment_volume,
        landscape_step.stored_sediment_volume,
        *landscape);

      if (should_write_surface_results())
        surface_results->write(
          this->get_output_directory(),
          this->get_timestep_number(),
          this->get_time() / year_in_seconds,
          write_visualization_results,
          output_drainage_diagnostics,
          surface_mesh,
          fastscape_points,
          *landscape,
          spatial_erosion_strength->get_values(),
          spatial_surface_runoff->get_values(),
          effective_ice.first,
          effective_ice.second,
          regional_total_ice_load_displacement,
          regional_ice_load_velocity);

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

      std::vector<double> elevation_values;
      std::vector<double> sediment_thickness_values;
      if (use_t_coupling)
        {
          elevation_values.assign(t_landscape->get_elevation().begin(),
                                  t_landscape->get_elevation().end());
          sediment_thickness_values.assign(
            t_landscape->get_sediment_thickness().begin(),
            t_landscape->get_sediment_thickness().end());
        }
      else
        {
          elevation_values.assign(landscape->get_elevation().begin(),
                                  landscape->get_elevation().end());
          sediment_thickness_values.assign(
            landscape->get_sediment_thickness().begin(),
            landscape->get_sediment_thickness().end());
        }
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
      if (use_t_coupling)
        for (const auto &field :
             t_landscape->get_sediment_thickness_by_lithology())
          sediment_by_lithology.emplace_back(field.begin(), field.end());
      else
        for (const auto &field :
             landscape->get_sediment_thickness_by_lithology())
          sediment_by_lithology.emplace_back(field.begin(), field.end());
      std::ostringstream lithology_stream;
      {
        aspect::oarchive archive(lithology_stream);
        archive << (use_t_coupling
                    ? t_landscape->get_lithology_names()
                    : landscape->get_lithology_names());
        archive << (use_t_coupling
                    ? t_landscape->get_bedrock_lithology()
                    : landscape->get_bedrock_lithology());
        archive << sediment_by_lithology;
      }
      status_strings["FastscapeLithologyProvenance"] =
        lithology_stream.str();

      std::ostringstream output_state_stream;
      {
        aspect::oarchive archive(output_state_stream);
        archive << (use_t_coupling
                    ? t_surface_results->get_reference_elevation()
                    : surface_results->get_reference_elevation());
        archive << (use_t_coupling
                    ? t_surface_results->get_output_history()
                    : surface_results->get_output_history());
      }
      status_strings["FastscapeSurfaceOutputState"] =
        output_state_stream.str();

      std::ostringstream output_timing_stream;
      {
        aspect::oarchive archive(output_timing_stream);
        archive << last_output_time;
      }
      status_strings["FastscapeSurfaceOutputTiming"] =
        output_timing_stream.str();

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

      std::ostringstream regional_load_stream;
      {
        aspect::oarchive archive(regional_load_stream);
        archive << regional_delayed_ice_load_displacement;
        archive << regional_total_ice_load_displacement;
        archive << regional_ice_load_velocity;
        archive << regional_ice_load_state_is_initialized;
      }
      status_strings["FastscapeRegionalIceLoadResponse"] =
        regional_load_stream.str();
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
      if (use_t_coupling)
        t_landscape->set_elevation(elevation_values);
      else
        landscape->set_elevation(elevation_values);

      const auto sediment_state =
        status_strings.find("FastscapeMarineSediment");
      if (sediment_state != status_strings.end())
        {
          std::vector<double> sediment_thickness_values;
          std::istringstream sediment_stream(sediment_state->second);
          aspect::iarchive sediment_archive(sediment_stream);
          sediment_archive >> sediment_thickness_values;
          if (use_t_coupling)
            t_landscape->set_sediment_thickness(sediment_thickness_values);
          else
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
          const std::vector<std::string> &active_lithology_names =
            use_t_coupling ? t_landscape->get_lithology_names()
            : landscape->get_lithology_names();
          AssertThrow(stored_names == active_lithology_names,
                      ExcMessage("The lithology names in the checkpoint differ "
                                 "from the current parameter file."));
          if (use_t_coupling)
            t_landscape->set_lithology_state(stored_bedrock_lithology,
                                             stored_sediment_thickness);
          else
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
          if (use_t_coupling)
            t_surface_results->restore_output_state(reference_elevation,
                                                    output_history);
          else
            surface_results->restore_output_state(reference_elevation,
                                                  output_history);
        }

      const auto output_timing_state =
        status_strings.find("FastscapeSurfaceOutputTiming");
      if (output_timing_state != status_strings.end())
        {
          std::istringstream output_timing_stream(output_timing_state->second);
          aspect::iarchive output_timing_archive(output_timing_stream);
          output_timing_archive >> last_output_time;
        }
      else if (output_interval > 0.0 && output_state != status_strings.end())
        {
          std::vector<double> reference_elevation;
          std::vector<std::pair<double,std::string>> output_history;
          std::istringstream output_state_stream(output_state->second);
          aspect::iarchive output_state_archive(output_state_stream);
          output_state_archive >> reference_elevation;
          output_state_archive >> output_history;
          if (!output_history.empty())
            {
              const double last_actual_output_time =
                output_history.back().first * year_in_seconds;
              const double elapsed_time =
                std::max(0.0, last_actual_output_time -
                         this->get_parameters().start_time);
              const double magic =
                1.0 + 2.0 * std::numeric_limits<double>::epsilon();
              last_output_time = this->get_parameters().start_time +
                                 std::floor(elapsed_time / output_interval * magic) *
                                 output_interval / magic;
            }
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

      const auto regional_load_state =
        status_strings.find("FastscapeRegionalIceLoadResponse");
      if (regional_load_state != status_strings.end())
        {
          std::istringstream regional_load_stream(regional_load_state->second);
          aspect::iarchive archive(regional_load_stream);
          archive >> regional_delayed_ice_load_displacement;
          archive >> regional_total_ice_load_displacement;
          archive >> regional_ice_load_velocity;
          archive >> regional_ice_load_state_is_initialized;
          const std::size_t active_size [[maybe_unused]] = use_t_coupling
                                                           ? t_landscape->get_elevation().size()
                                                           : landscape->get_elevation().size();
          AssertDimension(regional_delayed_ice_load_displacement.size(),
                          active_size);
          AssertDimension(regional_total_ice_load_displacement.size(),
                          active_size);
          AssertDimension(regional_ice_load_velocity.size(), active_size);
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
        prm.declare_entry("Use T coupling in 2d", "false", Patterns::Bool(),
                          "For a two-dimensional ASPECT Box model, evolve a "
                          "two-dimensional x-y landscape rather than a line. "
                          "ASPECT velocities are replicated across y and the "
                          "landscape response is reduced back to the ASPECT "
                          "section. This is analogous to the established "
                          "Fortran FastScape T coupling.");
        prm.declare_entry("T coupling width", "100000", Patterns::Double(0),
                          "Requested out-of-plane landscape width in meters. "
                          "The actual width is rounded to an integer number "
                          "of square cells with dy equal to dx.");
        prm.declare_entry("T coupling cell size", "0", Patterns::Double(0),
                          "Requested maximum cell size in meters for a 2-D "
                          "T-coupled landscape. A positive value determines "
                          "the number of x cells directly and constructs "
                          "square cells. Zero retains Box repetitions times "
                          "two to the Surface refinement level for backward "
                          "compatibility.");
        prm.declare_entry("T coupling reduction", "average",
                          Patterns::Selection("average|center"),
                          "Reduce the two-dimensional landscape response to "
                          "the ASPECT section using an area-weighted average "
                          "across y or the center row.");
        prm.declare_entry("Periodic surface dimensions", "",
                          Patterns::List(Patterns::Integer(0)),
                          "Comma-separated zero-based horizontal coordinate "
                          "directions to connect periodically on a box surface. "
                          "For a three-dimensional Cartesian model, use 1 to "
                          "make the transverse y direction periodic. The default "
                          "keeps all surface boundaries open/fixed as before.");
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
        prm.declare_entry("Flow routing method", "single flow",
                          Patterns::Selection("single flow|multiple flow"),
                          "Route water, sediment, and routed ice discharge to "
                          "one steepest-descent receiver or partition them "
                          "among all downslope neighbors. Multiple flow retains "
                          "single-flow routing only while resolving depressions, "
                          "then rebuilds the final graph with multiple receivers.");
        prm.declare_entry("Flow routing connectivity", "face",
                          Patterns::Selection("face|face and vertex"),
                          "Use only cells sharing a face or also cells sharing "
                          "a vertex as possible flow receivers. On a regular "
                          "two-dimensional quadrilateral surface these choices "
                          "correspond to D4 and D8 connectivity, respectively. "
                          "Periodic diagonal neighbors are included.");
        prm.declare_entry("Multiple flow slope exponent", "1.0",
                          Patterns::Double(0),
                          "Exponent p used to partition multiple-direction flow "
                          "in proportion to local slope raised to p. This entry "
                          "is ignored for single-flow routing. Larger values "
                          "concentrate flow more strongly along the steepest paths.");
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
        prm.declare_entry("Glacial ice thickness scale", "0",
                          Patterns::Double(0),
                          "Ice-thickness saturation scale H* in meters. If "
                          "positive, glacial erosion is multiplied by "
                          "1-exp(-H/H*), matching the Fortran coupling. Zero "
                          "retains the earlier thickness-independent law.");
        prm.declare_entry("Minimum ice thickness", "1",
                          Patterns::Double(0),
                          "Minimum ice thickness in meters required for "
                          "glacial erosion.");
        prm.declare_entry("Glacial erosion mode", "prescribed fields",
                          Patterns::Selection("prescribed fields|routed glacier"),
                          "Use externally supplied ice thickness and basal "
                          "velocity, or construct a low-cost glacier proxy "
                          "by routing elevation-dependent ice mass balance "
                          "over the evolving FastScape drainage graph.");
        prm.declare_entry("Routed glacier ELA", "1500",
                          Patterns::Double(),
                          "Equilibrium-line altitude in meters for the "
                          "routed-glacier proxy.");
        prm.declare_entry("Routed glacier accumulation gradient", "1e-3",
                          Patterns::Double(0),
                          "Ice-equivalent accumulation-rate increase in "
                          "meters per year per meter above the ELA.");
        prm.declare_entry("Routed glacier maximum accumulation", "2",
                          Patterns::Double(0),
                          "Maximum local ice accumulation in meters per year.");
        prm.declare_entry("Routed glacier ablation gradient", "2e-3",
                          Patterns::Double(0),
                          "Ice-equivalent ablation-rate increase in meters "
                          "per year per meter below the ELA.");
        prm.declare_entry("Routed glacier maximum ablation", "5",
                          Patterns::Double(0),
                          "Maximum local ice ablation in meters per year.");
        prm.declare_entry("Routed glacier minimum discharge", "1e6",
                          Patterns::Double(0),
                          "Minimum routed ice discharge in cubic meters per "
                          "year required to define a glacier corridor.");
        prm.declare_entry("Routed glacier reference discharge", "1e8",
                          Patterns::Double(0),
                          "Reference ice discharge in cubic meters per year "
                          "for the empirical width and thickness laws.");
        prm.declare_entry("Routed glacier reference width", "3000",
                          Patterns::Double(0),
                          "Glacier width in meters at the reference discharge.");
        prm.declare_entry("Routed glacier width exponent", "0.3",
                          Patterns::Double(0),
                          "Exponent in W=Wref*(Q/Qref)^p.");
        prm.declare_entry("Routed glacier minimum width", "1000",
                          Patterns::Double(0),
                          "Minimum empirical glacier width in meters.");
        prm.declare_entry("Routed glacier maximum width", "8000",
                          Patterns::Double(0),
                          "Maximum empirical glacier width in meters.");
        prm.declare_entry("Routed glacier reference thickness", "300",
                          Patterns::Double(0),
                          "Ice thickness in meters at the reference discharge "
                          "and a surface slope of 0.05.");
        prm.declare_entry("Routed glacier thickness exponent", "0.2",
                          Patterns::Double(0),
                          "Discharge exponent in the empirical ice-thickness law.");
        prm.declare_entry("Routed glacier minimum thickness", "50",
                          Patterns::Double(0),
                          "Minimum empirical ice thickness in meters.");
        prm.declare_entry("Routed glacier maximum thickness", "1200",
                          Patterns::Double(0),
                          "Maximum empirical ice thickness in meters.");
        prm.declare_entry("Routed glacier thickness slope exponent", "0.2",
                          Patterns::Double(0),
                          "Exponent applied to 0.05/max(surface_slope,1e-3) "
                          "in the empirical ice-thickness law.");
        prm.declare_entry("Spread routed glacier erosion over width", "true",
                          Patterns::Bool(),
                          "Distribute centerline erosion with a conservative "
                          "Gaussian graph-distance kernel whose diameter is "
                          "the empirical glacier width. Disable this option "
                          "to retain centerline-only erosion.");
        prm.declare_entry("Terminate routed glacier at sea level", "true",
                          Patterns::Bool(),
                          "Remove routed ice when it first reaches a cell at "
                          "or below sea level, representing coastal calving "
                          "and preventing ice discharge from following a flat "
                          "numerical ocean boundary. Disable this option for "
                          "grounded marine or ice-shelf experiments.");
        prm.declare_entry("Ice distribution mode", "spatial field",
                          Patterns::Selection(
                            "spatial field|elevation threshold|"
                            "spatial field and elevation threshold"),
                          "Choose prescribed climate-model fields, a PRM "
                          "elevation threshold, or climate fields masked by "
                          "the elevation threshold.");
        prm.declare_entry("Glacial elevation threshold", "1000",
                          Patterns::Double(),
                          "Surface elevation in meters above which the "
                          "elevation-threshold ice distribution is active.");
        prm.declare_entry("Glacial elevation transition width", "0",
                          Patterns::Double(0),
                          "Optional elevation interval in meters over which "
                          "ice thickness and basal velocity increase linearly "
                          "from zero above the threshold. Zero applies a hard "
                          "threshold.");
        prm.declare_entry("Elevation threshold ice thickness", "300",
                          Patterns::Double(0),
                          "Ice thickness in meters above the elevation "
                          "threshold when Ice distribution mode is elevation "
                          "threshold.");
        prm.declare_entry("Elevation threshold basal ice velocity", "0.1",
                          Patterns::Double(0),
                          "Basal sliding velocity in meters per year above "
                          "the elevation threshold when Ice distribution "
                          "mode is elevation threshold.");
        prm.declare_entry("Restrict glacial erosion to grounded ice", "true",
                          Patterns::Bool(),
                          "Suppress basal erosion where prescribed ice would "
                          "float according to rho_ice H < rho_water "
                          "(sea_level-bed_elevation).");
        prm.declare_entry("Minimum glacial erosion elevation", "-1e99",
                          Patterns::Double(),
                          "Optional lower bed-elevation limit in meters for "
                          "glacial erosion. The default leaves grounded marine "
                          "ice possible; use zero to forbid submarine erosion.");
        prm.declare_entry("Seawater density", "1028",
                          Patterns::Double(0),
                          "Seawater density in kilograms per cubic meter used "
                          "by the ice-flotation test.");
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
        prm.declare_entry("Maximum marine sediment transport Courant number", "0",
                          Patterns::Double(0),
                          "Maximum explicit Courant number used for marine "
                          "sediment diffusion. A positive value activates "
                          "automatic transport substeps; zero preserves the "
                          "single-step behavior used by earlier models.");
        prm.declare_entry("Limit marine deposition to available accommodation",
                          "false",
                          Patterns::Bool(),
                          "Limit coastal deposition to the pore-corrected "
                          "volume available below sea level plus Maximum "
                          "marine deposition above sea level. Sediment that "
                          "does not fit is conservatively counted as exported. "
                          "This prevents a drainage base-level cell from being "
                          "filled far above sea level in a single landscape step.");
        prm.declare_entry("Maximum marine deposition above sea level", "0",
                          Patterns::Double(0),
                          "Maximum elevation in meters above the current sea "
                          "level to which direct coastal sediment deposition "
                          "may fill a drainage base-level cell when the "
                          "accommodation limit is enabled.");
        prm.declare_entry("Use sea level as drainage base level", "false",
                          Patterns::Bool(),
                          "Stop river routing in connected cells at or below sea "
                          "level and deliver their sediment to the marine transport "
                          "model. This is always done on a closed spherical surface; "
                          "the default preserves the previous box behavior.");
        prm.declare_entry("Restrict ocean to largest connected water body",
                          "true",
                          Patterns::Bool(),
                          "Treat only the largest face-connected group of "
                          "below-sea-level cells as the global ocean. This "
                          "prevents disconnected inland depressions from "
                          "receiving marine sediment.");
        prm.declare_entry("Submarine river incision factor", "1",
                          Patterns::Double(0),
                          "Multiplier applied to stream-power incision below sea "
                          "level. Set this to zero when submarine channels are not "
                          "explicitly represented; one preserves the earlier behavior.");
        prm.declare_entry("Submarine hillslope diffusion coefficient", "-1",
                          Patterns::Double(-1),
                          "Hillslope diffusivity used between two ocean cells. A "
                          "negative value inherits Hillslope diffusion coefficient; "
                          "zero disables generic hillslope diffusion below sea level.");
        prm.declare_entry("Open marine sediment boundary", "none",
                          Patterns::Selection(
                            "none|minimum x|maximum x|all nonperiodic"),
                          "Absorbing box-surface boundary through which mobile marine "
                          "sediment leaves the landscape. Periodic sides can never "
                          "export sediment. The default closed basin preserves the "
                          "previous behavior.");
        prm.declare_entry("Advect surface state", "false",
                          Patterns::Bool(),
                          "Advect bedrock elevation and mobile sediment over "
                          "the fixed FastScape grid using the tangential "
                          "ASPECT material velocity.");
        prm.declare_entry("Apply normal material velocity", "true",
                          Patterns::Bool(),
                          "Apply the normal component of ASPECT's material "
                          "velocity as uplift or subsidence. Disable this only "
                          "to isolate tangential transport in an advection "
                          "benchmark; regional ice-load motion remains active.");
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
        prm.declare_entry("Enable regional ice load response", "false",
                          Patterns::Bool(),
                          "Add local bedrock subsidence and rebound caused by "
                          "the prescribed ice thickness to the uplift passed "
                          "through FastScape. This reduced response is intended "
                          "for regional box models and is disabled by default.");
        prm.declare_entry("Regional compensation density", "3300",
                          Patterns::Double(0),
                          "Density in kilograms per cubic meter that converts "
                          "ice mass per unit area into the local equilibrium "
                          "bedrock displacement.");
        prm.declare_entry("Regional immediate response fraction", "0",
                          Patterns::Double(0, 1),
                          "Fraction of the local equilibrium displacement "
                          "applied immediately. The remaining fraction follows "
                          "the regional ice-load relaxation time.");
        prm.declare_entry("Regional ice load relaxation time", "10000",
                          Patterns::Double(0),
                          "Time in years over which delayed regional subsidence "
                          "or rebound approaches local isostatic equilibrium. "
                          "Zero applies the delayed response immediately.");
        prm.declare_entry("Initialize regional ice load in equilibrium", "true",
                          Patterns::Bool(),
                          "Treat the initial prescribed ice as an already "
                          "compensated reference state. Disable this when the "
                          "initial ice load is emplaced at model start.");
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
                          "files when 'Time between graphical output' is zero. "
                          "Zero disables step-based result files.");
        prm.declare_entry("Time between graphical output", "0",
                          Patterns::Double(0),
                          "Simulation time between landscape surface result "
                          "files. If years are used in output, this value is "
                          "interpreted in years; otherwise it is interpreted "
                          "in seconds. A positive value takes precedence over "
                          "'Result interval'. Zero retains the legacy "
                          "step-based output schedule. Sediment-budget rows "
                          "are written every geodynamic timestep regardless "
                          "of this value.");
        prm.declare_entry("Write visualization results", "true",
                          Patterns::Bool(),
                          "Write cell-based elevation, erosion, drainage-area, "
                          "and sediment-flux visualization files.");
        prm.declare_entry("Output drainage diagnostics", "false",
                          Patterns::Bool(),
                          "Add runoff-weighted and geometric drainage areas, "
                          "dominant basin and outlet identifiers, primary "
                          "receiver diagnostics, and drainage/coastal outlet "
                          "masks to surface CSV and visualization output. "
                          "Also write drainage_basin_budget.csv. The legacy "
                          "drainage_area field is always retained.");
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
        use_t_coupling = prm.get_bool("Use T coupling in 2d");
        t_coupling_width = prm.get_double("T coupling width");
        t_coupling_cell_size = prm.get_double("T coupling cell size");
        t_coupling_reduction = prm.get("T coupling reduction");
        AssertThrow(!use_t_coupling || dim == 2,
                    ExcMessage("Use T coupling in 2d can only be enabled in a "
                               "two-dimensional ASPECT model."));
        periodic_surface_dimensions.clear();
        const std::string periodic_dimensions =
          prm.get("Periodic surface dimensions");
        if (!periodic_dimensions.empty())
          periodic_surface_dimensions = Utilities::string_to_unsigned_int(
                                          Utilities::split_string_list(periodic_dimensions));
        for (const unsigned int direction : periodic_surface_dimensions)
          AssertThrow(direction < (use_t_coupling ? 2u : dim-1),
                      ExcMessage("Periodic surface dimensions must be smaller "
                                 "than the number of horizontal dimensions."));
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
        flow_routing_method = prm.get("Flow routing method");
        flow_routing_connectivity =
          prm.get("Flow routing connectivity");
        multiple_flow_slope_exponent =
          prm.get_double("Multiple flow slope exponent");
        AssertThrow(flow_routing_method != "multiple flow" ||
                    std::abs(slope_exponent - 1.0) <=
                    std::numeric_limits<double>::epsilon(),
                    ExcMessage("FastScape multiple-flow routing currently "
                               "requires Slope exponent = 1.0."));
        nonlinear_tolerance = prm.get_double("Nonlinear tolerance");
        glacial_erosion_coefficient =
          prm.get_double("Glacial erosion coefficient");
        glacial_velocity_exponent =
          prm.get_double("Glacial velocity exponent");
        glacial_ice_thickness_scale =
          prm.get_double("Glacial ice thickness scale");
        AssertThrow(glacial_velocity_exponent > 0.0,
                    ExcMessage("Glacial velocity exponent must be positive."));
        minimum_ice_thickness =
          prm.get_double("Minimum ice thickness");
        glacial_erosion_mode = prm.get("Glacial erosion mode");
        routed_glacier_ela = prm.get_double("Routed glacier ELA");
        routed_glacier_accumulation_gradient =
          prm.get_double("Routed glacier accumulation gradient");
        routed_glacier_maximum_accumulation =
          prm.get_double("Routed glacier maximum accumulation");
        routed_glacier_ablation_gradient =
          prm.get_double("Routed glacier ablation gradient");
        routed_glacier_maximum_ablation =
          prm.get_double("Routed glacier maximum ablation");
        routed_glacier_minimum_discharge =
          prm.get_double("Routed glacier minimum discharge");
        routed_glacier_reference_discharge =
          prm.get_double("Routed glacier reference discharge");
        routed_glacier_reference_width =
          prm.get_double("Routed glacier reference width");
        routed_glacier_width_exponent =
          prm.get_double("Routed glacier width exponent");
        routed_glacier_minimum_width =
          prm.get_double("Routed glacier minimum width");
        routed_glacier_maximum_width =
          prm.get_double("Routed glacier maximum width");
        routed_glacier_reference_thickness =
          prm.get_double("Routed glacier reference thickness");
        routed_glacier_thickness_exponent =
          prm.get_double("Routed glacier thickness exponent");
        routed_glacier_minimum_thickness =
          prm.get_double("Routed glacier minimum thickness");
        routed_glacier_maximum_thickness =
          prm.get_double("Routed glacier maximum thickness");
        routed_glacier_thickness_slope_exponent =
          prm.get_double("Routed glacier thickness slope exponent");
        routed_glacier_spread_erosion =
          prm.get_bool("Spread routed glacier erosion over width");
        routed_glacier_terminate_at_sea_level =
          prm.get_bool("Terminate routed glacier at sea level");
        AssertThrow(routed_glacier_reference_discharge > 0.0,
                    ExcMessage("Routed glacier reference discharge must be positive."));
        AssertThrow(routed_glacier_minimum_width <=
                    routed_glacier_maximum_width,
                    ExcMessage("Routed glacier minimum width must not exceed its maximum width."));
        AssertThrow(routed_glacier_minimum_thickness <=
                    routed_glacier_maximum_thickness,
                    ExcMessage("Routed glacier minimum thickness must not exceed its maximum thickness."));
        ice_distribution_mode = prm.get("Ice distribution mode");
        glacial_elevation_threshold =
          prm.get_double("Glacial elevation threshold");
        glacial_elevation_transition_width =
          prm.get_double("Glacial elevation transition width");
        elevation_threshold_ice_thickness =
          prm.get_double("Elevation threshold ice thickness");
        elevation_threshold_basal_ice_velocity =
          prm.get_double("Elevation threshold basal ice velocity");
        restrict_glacial_erosion_to_grounded_ice =
          prm.get_bool("Restrict glacial erosion to grounded ice");
        minimum_glacial_erosion_elevation =
          prm.get_double("Minimum glacial erosion elevation");
        seawater_density = prm.get_double("Seawater density");
        AssertThrow(seawater_density > 0.0,
                    ExcMessage("Seawater density must be positive."));
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
        maximum_marine_sediment_transport_courant =
          prm.get_double("Maximum marine sediment transport Courant number");
        limit_marine_deposition_to_available_accommodation =
          prm.get_bool("Limit marine deposition to available accommodation");
        maximum_marine_deposition_above_sea_level =
          prm.get_double("Maximum marine deposition above sea level");
        use_sea_level_as_drainage_base_level =
          prm.get_bool("Use sea level as drainage base level");
        restrict_ocean_to_largest_connected_component =
          prm.get_bool("Restrict ocean to largest connected water body");
        submarine_river_incision_factor =
          prm.get_double("Submarine river incision factor");
        submarine_hillslope_diffusion_coefficient =
          prm.get_double("Submarine hillslope diffusion coefficient");
        open_marine_sediment_boundary =
          prm.get("Open marine sediment boundary");
        advect_surface_state =
          prm.get_bool("Advect surface state");
        apply_normal_material_velocity =
          prm.get_bool("Apply normal material velocity");
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
        regional_ice_load_response_enabled =
          prm.get_bool("Enable regional ice load response");
        regional_compensation_density =
          prm.get_double("Regional compensation density");
        AssertThrow(regional_compensation_density > 0.0,
                    ExcMessage("Regional compensation density must be positive."));
        regional_immediate_response_fraction =
          prm.get_double("Regional immediate response fraction");
        regional_ice_load_relaxation_time =
          prm.get_double("Regional ice load relaxation time");
        initialize_regional_ice_load_in_equilibrium =
          prm.get_bool("Initialize regional ice load in equilibrium");
        true_polar_wander_enabled =
          prm.get_bool("Enable true polar wander");
        include_ice_load_in_true_polar_wander =
          prm.get_bool("Include ice load in true polar wander");
        degree_two_self_gravity_enabled =
          prm.get_bool("Enable degree two self gravity");
        AssertThrow(!(regional_ice_load_response_enabled &&
                      degree_two_self_gravity_enabled),
                    ExcMessage("The regional ice-load response and the global "
                               "degree-two self-gravity correction cannot be "
                               "enabled together because this would count the "
                               "solid-Earth load response twice."));
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
        output_interval = prm.get_double("Time between graphical output");
        if (this->convert_output_to_years())
          output_interval *= year_in_seconds;
        write_visualization_results =
          prm.get_bool("Write visualization results");
        output_drainage_diagnostics =
          prm.get_bool("Output drainage diagnostics");
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
