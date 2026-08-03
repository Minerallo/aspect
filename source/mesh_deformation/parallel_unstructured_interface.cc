/*
  Copyright (C) 2014 - 2026 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/parallel_unstructured_interface.h>
#include <aspect/simulator_signals.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools_evaluate.h>
#include <deal.II/numerics/rtree.h>

#include <boost/geometry/index/predicates.hpp>


namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    void
    ParallelUnstructuredInterface<dim>::
    set_surface_transfer_options(const std::string &scheme,
                                 const unsigned int neighbors,
                                 const bool normalize_coordinates)
    {
      if (scheme == "nearest")
        surface_transfer_scheme = SurfaceTransferScheme::nearest;
      else if (scheme == "weighted")
        surface_transfer_scheme = SurfaceTransferScheme::weighted;
      else if (scheme == "conservative")
        surface_transfer_scheme = SurfaceTransferScheme::conservative;
      else
        AssertThrow(false,
                    ExcMessage("Unknown surface transfer scheme '" + scheme + "'."));

      surface_transfer_neighbors = std::max(1u, neighbors);
      normalize_transfer_coordinates = normalize_coordinates;
    }



    template <int dim>
    void
    ParallelUnstructuredInterface<dim>::
    set_evaluation_point_areas(const std::vector<double> &areas)
    {
      evaluation_point_areas = areas;
    }

    template <int dim>
    void
    ParallelUnstructuredInterface<dim>::
    compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                             AffineConstraints<double> &mesh_velocity_constraints,
                                             const std::set<types::boundary_id> &boundary_ids) const
    {
      // First compute a (global) vector that has the correct ASPECT solution at all boundary nodes.
      const std::vector<std::vector<double>> aspect_surface_solution = evaluate_aspect_solution_at_points();

      // Extract the velocity from the external tool. This is implemented in the derived class.
      const std::vector<Tensor<1,dim>> external_surface_velocities
        = compute_updated_velocities_at_points(aspect_surface_solution);

      // Take the velocities computed above and constrain the DoF's on the surface of ASPECT's mesh.
      const LinearAlgebra::Vector v_interpolated
        = interpolate_external_velocities_to_surface_support_points(external_surface_velocities);

      // Create a ghosted vector that contains the interpolated velocities. This is necessary
      // because the constraints are set on all locally relevant DoFs (shared DoFs), not just
      // the locally owned ones.
      const DoFHandler<dim> &mesh_dof_handler = this->get_mesh_deformation_handler().get_mesh_deformation_dof_handler();
      const IndexSet mesh_locally_relevant = DoFTools::extract_locally_relevant_dofs (mesh_dof_handler);
      LinearAlgebra::Vector v_interpolated_ghosted(mesh_dof_handler.locally_owned_dofs(),
                                                   mesh_locally_relevant,
                                                   this->get_mpi_communicator());
      v_interpolated_ghosted = v_interpolated;

      // Set the constraints. For this, loop over all boundary DoFs
      // and if a boundary DoF is locally owned, create a constraint.
      // We later make that consistent across processors to ensure we
      // also know about the locally relevant DoFs constraints.
      // Now insert the relevant part of the solution into the mesh constraints
      const IndexSet constrained_dofs =
        DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                        ComponentMask(dim, true),
                                        boundary_ids);

      for (const types::global_dof_index index : constrained_dofs)
        {
          if (mesh_velocity_constraints.can_store_line(index))
            if (mesh_velocity_constraints.is_constrained(index)==false)
              {
#if DEAL_II_VERSION_GTE(9,6,0)
                mesh_velocity_constraints.add_constraint(index,
                                                         {},
                                                         v_interpolated_ghosted(index));
#else
                mesh_velocity_constraints.add_line(index);
                mesh_velocity_constraints.set_inhomogeneity(index, v_interpolated_ghosted(index));
#endif
              }
        }
    }


    template <int dim>
    void
    ParallelUnstructuredInterface<dim>::
    set_evaluation_points (const std::vector<Point<dim>> &evaluation_points)
    {
      // First, save a copy of the points at which we need the solution,
      // among other reasons so that we can track that input arguments
      // for later function calls describe the same number of points.
      this->evaluation_points = evaluation_points;

      // Set up RemotePointEvaluation. The evaluation points are given in reference coordinates,
      // so we need to use a simple mapping instead of the MappingQEulerian stored in the Simulator. The latter
      // would produce the deformed mesh. We pick a high order mapping for curved meshes (sphere, shell),
      // in order to find evaluation points on the surface of the sphere/shell. For cartesian meshes, we
      // can use a Q1 mapping.
      static MappingQ<dim> mapping(this->get_geometry_model().has_curved_elements() ? 4 : 1);
      remote_point_evaluator = std::make_unique<Utilities::MPI::RemotePointEvaluation<dim, dim>>();
      remote_point_evaluator->reinit(this->evaluation_points, this->get_triangulation(), mapping);

      if (!remote_point_evaluator->all_points_found())
        {
          this->get_pcout() << "WARNING: not all evaluation points were found inside the domain!" << std::endl;
          this->get_pcout() << "Evaluation points not found:" << std::endl;
          for (unsigned int p=0; p<evaluation_points.size(); ++p)
            {
              if (!remote_point_evaluator->point_found(p))
                {
                  this->get_pcout() << "Point " << p << ": " << evaluation_points[p] << std::endl;
                }
            }
        }

      // Create a global mapping from external evaluation points to every
      // ASPECT surface support point. This must not be limited to points that
      // happen to lie in the same volume cell: an independent external mesh
      // may be coarser or have a different topology.
      {
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
        const DoFHandler<dim> &mesh_dof_handler = this->get_mesh_deformation_handler().get_mesh_deformation_dof_handler();
        const auto boundary_ids = this->get_mesh_deformation_boundary_indicators();
        const unsigned int dofs_per_cell = mesh_dof_handler.get_fe().dofs_per_cell;
        const std::vector<Point<dim>> &unit_support_points =
          mesh_dof_handler.get_fe().get_unit_support_points();
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        std::vector<SurfaceSupportPointData> local_support_points;

        for (const auto &cell : mesh_dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              cell->get_dof_indices(local_dof_indices);
              for (const unsigned int face : cell->face_indices())
                if (cell->face(face)->at_boundary() &&
                    boundary_ids.find(cell->face(face)->boundary_id()) != boundary_ids.end())
                  {
                    const unsigned int coordinate = face / 2;
                    const double side = face % 2;
                    unsigned int scalar_support_points_on_face = 0;
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      if (mesh_dof_handler.get_fe().system_to_component_index(j).first == 0
                          && std::abs(unit_support_points[j][coordinate] - side) < 1e-12)
                        ++scalar_support_points_on_face;

                    Assert(scalar_support_points_on_face > 0, ExcInternalError());
                    const double nodal_area =
                      cell->face(face)->measure() / scalar_support_points_on_face;
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      if (std::abs(unit_support_points[j][coordinate] - side) < 1e-12)
                        local_support_points.push_back(
                          {local_dof_indices[j],
                           mesh_dof_handler.get_fe().system_to_component_index(j).first,
                           mapping.transform_unit_to_real_cell(cell, unit_support_points[j]),
                           nodal_area});
                  }
            }

        const auto gathered_support_points =
          Utilities::MPI::gather(this->get_mpi_communicator(),
                                 local_support_points,
                                 0);
        gathered_evaluation_points =
          Utilities::MPI::gather(this->get_mpi_communicator(),
                                 this->evaluation_points,
                                 0);
        gathered_evaluation_point_areas =
          Utilities::MPI::gather(this->get_mpi_communicator(),
                                 evaluation_point_areas,
                                 0);

        map_dof_to_eval_point.clear();
        if (my_rank == 0)
          {
            struct SourcePoint
            {
              Point<dim> point;
              unsigned int rank;
              unsigned int index;
            };

            std::vector<SourcePoint> sources;
            std::vector<std::pair<Point<dim>,unsigned int>> tree_entries;
            for (unsigned int rank = 0; rank < gathered_evaluation_points.size(); ++rank)
              for (unsigned int index = 0;
                   index < gathered_evaluation_points[rank].size();
                   ++index)
                {
                  Point<dim> point = gathered_evaluation_points[rank][index];
                  if (normalize_transfer_coordinates && point.norm() > 0.0)
                    point /= point.norm();
                  sources.push_back({point, rank, index});
                  tree_entries.emplace_back(point, sources.size()-1);
                }

            AssertThrow(!sources.empty(),
                        ExcMessage("No external surface evaluation points were provided."));
            const auto tree = pack_rtree(tree_entries);
            // Merge repeated contributions to nodal control areas, including
            // support points shared by surface faces and MPI subdomains.
            std::map<types::global_dof_index,SurfaceSupportPointData> targets;
            for (const auto &rank_points : gathered_support_points)
              for (const auto &target : rank_points)
                {
                  const auto position = targets.find(target.dof_index);
                  if (position == targets.end())
                    targets.emplace(target.dof_index, target);
                  else
                    position->second.area += target.area;
                }

            namespace bgi = boost::geometry::index;
            const unsigned int requested_neighbors =
              (surface_transfer_scheme == SurfaceTransferScheme::nearest
               ? 1
               : surface_transfer_neighbors);
            const unsigned int n_neighbors =
              std::min<unsigned int>(requested_neighbors, sources.size());

            for (const auto &[dof_index, target] : targets)
              {
                Point<dim> target_point = target.point;
                if (normalize_transfer_coordinates && target_point.norm() > 0.0)
                  target_point /= target_point.norm();

                std::vector<std::pair<Point<dim>,unsigned int>> nearest;
                tree.query(bgi::nearest(target_point, n_neighbors),
                           std::back_inserter(nearest));
                std::vector<double> weights(nearest.size(), 0.0);
                unsigned int exact_neighbor = numbers::invalid_unsigned_int;
                double weight_sum = 0.0;
                for (unsigned int i = 0; i < nearest.size(); ++i)
                  {
                    const double distance_squared =
                      target_point.distance_square(nearest[i].first);
                    if (distance_squared < 1e-28)
                      exact_neighbor = i;
                    else
                      {
                        weights[i] = 1.0 / distance_squared;
                        weight_sum += weights[i];
                      }
                  }
                if (exact_neighbor != numbers::invalid_unsigned_int)
                  {
                    std::fill(weights.begin(), weights.end(), 0.0);
                    weights[exact_neighbor] = 1.0;
                  }
                else
                  for (double &weight : weights)
                    weight /= weight_sum;

                const double normal_component =
                  target.point.norm() > 0.0
                  ? target.point[target.component] / target.point.norm()
                  : 0.0;
                for (unsigned int i = 0; i < nearest.size(); ++i)
                  {
                    const SourcePoint &source = sources[nearest[i].second];
                    map_dof_to_eval_point.push_back(
                      {dof_index, source.rank, source.index, target.component,
                       target_point.distance_square(source.point), weights[i],
                       normal_component, target.area});
                  }
              }
          }
      }

      // Finally, also ensure that upon mesh refinement, all of the
      // information set herein is invalidated:
      this->get_signals().pre_refinement_store_user_data
      .connect([this](typename parallel::distributed::Triangulation<dim> &)
      {
        this->evaluation_points.clear();
        this->evaluation_point_areas.clear();
        this->gathered_evaluation_points.clear();
        this->gathered_evaluation_point_areas.clear();
        this->map_dof_to_eval_point.clear();
        this->remote_point_evaluator.reset();
      });
    }



    template <int dim>
    std::vector<std::vector<double>>
    ParallelUnstructuredInterface<dim>::
    evaluate_aspect_solution_at_points () const
    {
      Assert (remote_point_evaluator != nullptr,
              ExcMessage("You can only call this function if you have previously "
                         "set the evaluation points by calling set_evaluation_points(), "
                         "and if the evaluator has not been invalidated by a mesh "
                         "refinement step."));

      // All components are evaluated (velocity, pressure, temperature, and N compositional fields).
      const unsigned int n_components = this->introspection().n_components;
      std::vector<std::vector<double>> solution_at_points (evaluation_points.size(), std::vector<double>(n_components, 0.0));

      // VectorTools::point_values can evaluate N components at a time, but this is a template argument and not a
      // runtime argument. For now, we just evaluate them one component at a time. Of course it would be more
      // efficient to branch and evaluate up to K at a time (for a reasonable number of K, say 10). Maybe something
      // to implement directly into deal.II.
      for (unsigned int c=0; c<n_components; ++c)
        {
          const std::vector<double> values = VectorTools::point_values<1>(*this->remote_point_evaluator,
                                                                          this->get_dof_handler(),
                                                                          this->get_solution(),
                                                                          dealii::VectorTools::EvaluationFlags::avg,
                                                                          c);
          for (unsigned int i=0; i<evaluation_points.size(); ++i)
            solution_at_points[i][c] = values[i];
        }

      return solution_at_points;
    }



    template <int dim>
    LinearAlgebra::Vector
    ParallelUnstructuredInterface<dim>::
    interpolate_external_velocities_to_surface_support_points (const std::vector<Tensor<1,dim>> &velocities) const
    {
      Assert (remote_point_evaluator != nullptr,
              ExcMessage("You can only call this function if you have previously "
                         "set the evaluation points by calling set_evaluation_points(), "
                         "and if the evaluator has not been invalidated by a mesh "
                         "refinement step."));
      AssertDimension(velocities.size(), evaluation_points.size());


      // Create the output vector.
      const DoFHandler<dim> &mesh_dof_handler = this->get_mesh_deformation_handler().get_mesh_deformation_dof_handler();
      LinearAlgebra::Vector vector_with_surface_velocities(mesh_dof_handler.locally_owned_dofs(),
                                                           this->get_mpi_communicator());

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      const auto gathered_velocities =
        Utilities::MPI::gather(this->get_mpi_communicator(),
                               velocities,
                               0);

      if (my_rank == 0)
        {
          std::map<types::global_dof_index,double> transferred_values;
          if (normalize_transfer_coordinates)
            {
              // Interpolate scalar normal speeds, then reconstruct Cartesian
              // components along each target support point's own radial
              // direction. Interpolating source Cartesian components directly
              // introduces artificial tangential motion on a sphere.
              std::map<types::global_dof_index,double> normal_speeds;
              for (const auto &entry : map_dof_to_eval_point)
                {
                  const Point<dim> &source_point =
                    gathered_evaluation_points[entry.evaluation_point_rank]
                                               [entry.evaluation_point_index];
                  Tensor<1,dim> source_normal;
                  for (unsigned int d = 0; d < dim; ++d)
                    source_normal[d] = source_point[d] / source_point.norm();
                  normal_speeds[entry.dof_index] +=
                    entry.weight
                    * (gathered_velocities[entry.evaluation_point_rank]
                                            [entry.evaluation_point_index]
                       * source_normal);
                }

              double positive_scale = 1.0;
              double negative_scale = 1.0;
              if (surface_transfer_scheme == SurfaceTransferScheme::conservative)
                {
                  double source_area = 0.0;
                  double source_positive = 0.0;
                  double source_negative = 0.0;
                  for (unsigned int rank = 0;
                       rank < gathered_evaluation_points.size();
                       ++rank)
                    {
                      AssertThrow(gathered_evaluation_point_areas[rank].size()
                                  == gathered_evaluation_points[rank].size(),
                                  ExcMessage("The conservative surface transfer requires "
                                             "one positive area for every external "
                                             "evaluation point."));
                      for (unsigned int i = 0;
                           i < gathered_evaluation_points[rank].size();
                           ++i)
                        {
                          const Point<dim> &point =
                            gathered_evaluation_points[rank][i];
                          Tensor<1,dim> normal;
                          for (unsigned int d = 0; d < dim; ++d)
                            normal[d] = point[d] / point.norm();
                          const double value =
                            gathered_velocities[rank][i] * normal;
                          const double area =
                            gathered_evaluation_point_areas[rank][i];
                          AssertThrow(area >= 0.0,
                                      ExcMessage("External surface point areas "
                                                 "must be nonnegative."));
                          source_area += area;
                          source_positive += area * std::max(value, 0.0);
                          source_negative += area * std::max(-value, 0.0);
                        }
                    }

                  double target_area = 0.0;
                  double target_positive = 0.0;
                  double target_negative = 0.0;
                  std::set<types::global_dof_index> counted_target_dofs;
                  for (const auto &entry : map_dof_to_eval_point)
                    if (entry.component == 0
                        && counted_target_dofs.insert(entry.dof_index).second)
                      {
                        const double value = normal_speeds[entry.dof_index];
                        target_area += entry.target_area;
                        target_positive += entry.target_area * std::max(value, 0.0);
                        target_negative += entry.target_area * std::max(-value, 0.0);
                      }

                  AssertThrow(source_area > 0.0 && target_area > 0.0,
                              ExcMessage("The conservative surface transfer requires "
                                         "positive source and target areas."));
                  const double source_area_normalization = target_area / source_area;
                  source_positive *= source_area_normalization;
                  source_negative *= source_area_normalization;
                  if (source_positive > 0.0)
                    AssertThrow(target_positive > 0.0,
                                ExcMessage("The target stencil lost all positive "
                                           "normal surface motion."));
                  if (source_negative > 0.0)
                    AssertThrow(target_negative > 0.0,
                                ExcMessage("The target stencil lost all negative "
                                           "normal surface motion."));
                  positive_scale =
                    source_positive > 0.0 ? source_positive / target_positive : 0.0;
                  negative_scale =
                    source_negative > 0.0 ? source_negative / target_negative : 0.0;
                }

              for (const auto &entry : map_dof_to_eval_point)
                {
                  double normal_speed = normal_speeds[entry.dof_index];
                  normal_speed *= normal_speed >= 0.0
                                  ? positive_scale
                                  : negative_scale;
                  transferred_values[entry.dof_index] =
                    normal_speed * entry.target_normal_component;
                }
            }
          else
            {
              for (const auto &entry : map_dof_to_eval_point)
                transferred_values[entry.dof_index] +=
                  entry.weight
                  * gathered_velocities[entry.evaluation_point_rank]
                                         [entry.evaluation_point_index]
                                         [entry.component];

              if (surface_transfer_scheme == SurfaceTransferScheme::conservative)
                for (unsigned int component = 0; component < dim; ++component)
                  {
                    double source_area = 0.0;
                    double source_positive = 0.0;
                    double source_negative = 0.0;
                    for (unsigned int rank = 0;
                         rank < gathered_evaluation_points.size();
                         ++rank)
                      {
                        AssertThrow(gathered_evaluation_point_areas[rank].size()
                                    == gathered_evaluation_points[rank].size(),
                                    ExcMessage("The conservative surface transfer requires "
                                               "one positive area for every external "
                                               "evaluation point."));
                        for (unsigned int i = 0;
                             i < gathered_evaluation_points[rank].size();
                             ++i)
                          {
                            const double area =
                              gathered_evaluation_point_areas[rank][i];
                            AssertThrow(area >= 0.0,
                                        ExcMessage("External surface point areas "
                                                   "must be nonnegative."));
                            const double value =
                              gathered_velocities[rank][i][component];
                            source_area += area;
                            source_positive += area * std::max(value, 0.0);
                            source_negative += area * std::max(-value, 0.0);
                          }
                      }

                    double target_area = 0.0;
                    double target_positive = 0.0;
                    double target_negative = 0.0;
                    std::set<types::global_dof_index> counted_target_dofs;
                    for (const auto &entry : map_dof_to_eval_point)
                      if (entry.component == component
                          && counted_target_dofs.insert(entry.dof_index).second)
                        {
                          const double value = transferred_values[entry.dof_index];
                          target_area += entry.target_area;
                          target_positive += entry.target_area * std::max(value, 0.0);
                          target_negative += entry.target_area * std::max(-value, 0.0);
                        }

                    AssertThrow(source_area > 0.0 && target_area > 0.0,
                                ExcMessage("The conservative surface transfer requires "
                                           "positive source and target areas."));
                    const double source_area_normalization = target_area / source_area;
                    source_positive *= source_area_normalization;
                    source_negative *= source_area_normalization;
                    if (source_positive > 0.0)
                      AssertThrow(target_positive > 0.0,
                                  ExcMessage("The target stencil lost all positive "
                                             "surface motion."));
                    if (source_negative > 0.0)
                      AssertThrow(target_negative > 0.0,
                                  ExcMessage("The target stencil lost all negative "
                                             "surface motion."));
                    const double positive_scale =
                      source_positive > 0.0 ? source_positive / target_positive : 0.0;
                    const double negative_scale =
                      source_negative > 0.0 ? source_negative / target_negative : 0.0;

                    counted_target_dofs.clear();
                    for (const auto &entry : map_dof_to_eval_point)
                      if (entry.component == component
                          && counted_target_dofs.insert(entry.dof_index).second)
                        {
                          double &value = transferred_values[entry.dof_index];
                          value *= value >= 0.0 ? positive_scale : negative_scale;
                        }
                  }
            }

          for (const auto &[dof_index, value] : transferred_values)
            vector_with_surface_velocities[dof_index] = value;
        }

      vector_with_surface_velocities.compress(VectorOperation::insert);

      return vector_with_surface_velocities;
    }
  }



  namespace MeshDeformation
  {
#define INSTANTIATE(dim) \
  template class ParallelUnstructuredInterface<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
