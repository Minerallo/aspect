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

#ifndef _aspect_mesh_deformation_fastscape_cpp_grid_h
#define _aspect_mesh_deformation_fastscape_cpp_grid_h

#include <aspect/config.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <deal.II/grid/tria.h>
#include <deal.II/base/tensor.h>

#include <fastscapelib/grid/base.hpp>
#include <fastscapelib/version.hpp>

#include <array>
#include <unordered_map>
#include <vector>

#if __has_include(<xtensor/core/xtensor_config.hpp>)
#  include <xtensor/core/xtensor_config.hpp>
#  include <xtensor/containers/xarray.hpp>
#  include <xtensor/containers/xtensor.hpp>
#else
#  include <xtensor/xtensor_config.hpp>
#  include <xtensor/xarray.hpp>
#  include <xtensor/xtensor.hpp>
#endif

#if FASTSCAPELIB_VERSION_MAJOR * 10000 + FASTSCAPELIB_VERSION_MINOR * 100 + FASTSCAPELIB_VERSION_PATCH <= 202
#  define ASPECT_FASTSCAPELIB_LEGACY_GRID_API 1
#  include <fastscapelib/utils/xtensor_utils.hpp>
#else
#  define ASPECT_FASTSCAPELIB_LEGACY_GRID_API 0
#  include <fastscapelib/utils/containers.hpp>
#endif


namespace fastscapelib
{
  /**
   * An adapter that exposes the active cells of a serial deal.II surface
   * triangulation as an unstructured FastScape grid. One FastScape node
   * corresponds to one surface cell. Unlike the adapter prototyped in
   * geodynamics/aspect#6389, this implementation retains the individual cell
   * areas. This is important on a spherical surface, where the cells are not
   * exactly equal-area.
   */
  template <class TriangulationType, class Selector>
  class dealii_surface_grid;

  template <class TriangulationType, class Selector>
  struct grid_inner_types<dealii_surface_grid<TriangulationType, Selector>>
  {
    static constexpr bool is_structured = false;
    static constexpr bool is_uniform = false;
    using grid_data_type = double;
#if ASPECT_FASTSCAPELIB_LEGACY_GRID_API
    using xt_selector = Selector;
    static constexpr std::size_t xt_ndims = 1;
#else
    using container_selector = Selector;
    static constexpr std::size_t container_ndims = 1;
#endif
    static constexpr uint8_t n_neighbors_max = 8;
    using neighbors_cache_type = neighbors_no_cache<n_neighbors_max>;
  };


  template <class TriangulationType, class Selector = xt_selector>
  class dealii_surface_grid
    : public grid<dealii_surface_grid<TriangulationType, Selector>>
  {
    public:
      using self_type = dealii_surface_grid<TriangulationType, Selector>;
      using base_type = grid<self_type>;
      using inner_types = grid_inner_types<self_type>;
      using grid_data_type = typename base_type::grid_data_type;
#if ASPECT_FASTSCAPELIB_LEGACY_GRID_API
      using xt_selector = typename base_type::xt_selector;
      using container_type = xt_tensor_t<xt_selector, grid_data_type, 1>;
#else
      using container_selector = typename base_type::container_selector;
      using container_type =
        fixed_shape_container_t<container_selector, grid_data_type, 1>;
#endif
      using nodes_status_type = typename base_type::nodes_status_type;
      using size_type = typename base_type::size_type;
      using shape_type = typename base_type::shape_type;
      using point_type = dealii::Point<TriangulationType::space_dimension>;
      using direction_type = dealii::Tensor<1,TriangulationType::space_dimension>;
      using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
      using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

      /**
       * Build the adapter. If @p closed_surface is false, cells touching the
       * boundary are fixed-value/base-level nodes. Closed surfaces (global
       * spheres) receive their base levels dynamically from the coupling
       * plugin.
       */
      dealii_surface_grid(TriangulationType &triangulation,
                          const bool closed_surface,
                          const bool include_vertex_neighbors = false);

      size_type
      number_of_cell_neighbors(const size_type index) const;

      size_type
      cell_neighbor(const size_type index,
                    const size_type neighbor_number) const;

      grid_data_type
      cell_neighbor_distance(const size_type index,
                             const size_type neighbor_number) const;

      grid_data_type
      cell_shared_face_measure(const size_type index,
                               const size_type neighbor_number) const;

      /** Return whether a cell touches a non-periodic exterior boundary. */
      bool
      cell_touches_boundary(const size_type index,
                            const unsigned int coordinate_direction,
                            const bool upper_side) const;

      direction_type
      cell_neighbor_direction(const size_type index,
                              const size_type neighbor_number) const;

    protected:
      TriangulationType &triangulation;
      const bool closed_surface;
      shape_type m_shape;
      size_type m_size;
      container_type node_areas;
      nodes_status_type m_nodes_status;
      std::vector<point_type> cell_centers;
      std::vector<size_type> neighbors_count;
      std::vector<neighbors_indices_impl_type> neighbors_indices;
      std::vector<neighbors_distances_impl_type> neighbors_distances;
      std::vector<neighbors_distances_impl_type> shared_face_measures;
      std::vector<std::array<bool,2*TriangulationType::dimension>> boundary_sides;
      std::vector<std::array<direction_type,inner_types::n_neighbors_max>> neighbor_directions;

      container_type nodes_areas_impl() const;
      grid_data_type nodes_areas_impl(const size_type &index) const noexcept;
      size_type neighbors_count_impl(const size_type &index) const;
      void neighbors_indices_impl(neighbors_indices_impl_type &neighbors,
                                  const size_type &index) const;
      const neighbors_distances_impl_type &
      neighbors_distances_impl(const size_type &index) const;

      friend class grid<self_type>;
  };


  template <class TriangulationType, class Selector>
  dealii_surface_grid<TriangulationType, Selector>::
  dealii_surface_grid(TriangulationType &input_triangulation,
                      const bool closed_surface,
                      const bool include_vertex_neighbors)
    :
    base_type(0),
    triangulation(input_triangulation),
    closed_surface(closed_surface),
    m_size(triangulation.n_active_cells()),
    node_areas(container_type::from_shape({m_size})),
             m_nodes_status(nodes_status_type::from_shape({m_size})),
             cell_centers(m_size),
             neighbors_count(m_size),
             neighbors_indices(m_size),
             neighbors_distances(m_size),
             shared_face_measures(m_size),
             boundary_sides(m_size),
             neighbor_directions(m_size)
  {
    m_shape = {{static_cast<typename shape_type::value_type>(m_size)}};
    std::fill(m_nodes_status.begin(), m_nodes_status.end(), node_status::core);

    std::unordered_map<unsigned int, size_type> global_to_local;
    std::unordered_map<unsigned int, std::vector<size_type>> vertex_to_cells;
    size_type local_index = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        global_to_local[cell->global_active_cell_index()] = local_index;
        if (include_vertex_neighbors)
          for (const unsigned int vertex : cell->vertex_indices())
            vertex_to_cells[cell->vertex_index(vertex)].push_back(local_index);
        node_areas[local_index] = cell->measure();
        cell_centers[local_index] = cell->center();

        boundary_sides[local_index].fill(false);
        if (!closed_surface)
          for (const unsigned int face : cell->face_indices())
            if (cell->at_boundary(face) && !cell->has_periodic_neighbor(face))
              {
                m_nodes_status[local_index] = node_status::fixed_value;
                const unsigned int boundary_id = cell->face(face)->boundary_id();
                if (boundary_id < 2*TriangulationType::dimension)
                  boundary_sides[local_index][boundary_id] = true;
              }
        ++local_index;
      }

    local_index = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        const dealii::Point<TriangulationType::space_dimension> center = cell->center();
        size_type count = 0;

        const auto add_neighbor = [&] (const size_type neighbor_index,
                                       const grid_data_type distance,
                                       const direction_type &direction,
                                       const grid_data_type shared_measure)
        {
          if (neighbor_index == local_index)
            return;
          for (size_type existing = 0; existing < count; ++existing)
            if (neighbors_indices[local_index][existing] == neighbor_index)
              return;

          AssertThrow(count < inner_types::n_neighbors_max,
                      dealii::ExcMessage(
                        "FastScape vertex connectivity produced more than "
                        "eight neighbors for a surface cell."));
          neighbors_indices[local_index][count] = neighbor_index;
          neighbors_distances[local_index][count] = distance;
          neighbor_directions[local_index][count] = direction;
          shared_face_measures[local_index][count] = shared_measure;
          ++count;
        };

        for (const unsigned int face : cell->face_indices())
          if (!cell->at_boundary(face) || cell->has_periodic_neighbor(face))
            {
              const bool periodic = cell->has_periodic_neighbor(face);
              const auto neighbor = cell->neighbor_or_periodic_neighbor(face);
              const auto entry = global_to_local.find(neighbor->global_active_cell_index());
              AssertThrow(entry != global_to_local.end(),
                          dealii::ExcMessage("FastScape surface-grid neighbor was not found."));
              grid_data_type distance;
              direction_type direction;
              if (periodic)
                {
                  distance =
                    center.distance(cell->face(face)->center()) +
                    neighbor->center().distance(
                      neighbor->face(cell->periodic_neighbor_of_periodic_neighbor(face))->center());
                  direction = cell->face(face)->center() - center;
                }
              else
                {
                  distance = center.distance(neighbor->center());
                  direction = neighbor->center() - center;
                }
              grid_data_type shared_measure;
              if constexpr (TriangulationType::dimension == 1)
                shared_measure = 1.0;
              else
                shared_measure = cell->face(face)->measure();
              add_neighbor(entry->second, distance, direction, shared_measure);
            }

        if (include_vertex_neighbors)
          {
            for (const unsigned int vertex : cell->vertex_indices())
              for (const size_type neighbor_index :
                   vertex_to_cells[cell->vertex_index(vertex)])
                add_neighbor(neighbor_index,
                             center.distance(cell_centers[neighbor_index]),
                             cell_centers[neighbor_index] - center,
                             0.0);

            // Vertices on opposite periodic boundaries are not identified in
            // the triangulation. Add the two diagonal cells across each
            // periodic face explicitly, including doubly-periodic corners.
            if constexpr (TriangulationType::dimension == 2)
              for (const unsigned int face : cell->face_indices())
                if (cell->has_periodic_neighbor(face))
                  {
                    const auto periodic_neighbor =
                      cell->neighbor_or_periodic_neighbor(face);
                    const unsigned int opposite_face =
                      cell->periodic_neighbor_of_periodic_neighbor(face);
                    const grid_data_type periodic_distance =
                      center.distance(cell->face(face)->center()) +
                      periodic_neighbor->center().distance(
                        periodic_neighbor->face(opposite_face)->center());
                    direction_type periodic_direction =
                      cell->face(face)->center() - center;
                    periodic_direction *=
                      periodic_distance / periodic_direction.norm();

                    for (const unsigned int transverse_face :
                         periodic_neighbor->face_indices())
                      if (transverse_face / 2 != opposite_face / 2 &&
                          (!periodic_neighbor->at_boundary(transverse_face) ||
                           periodic_neighbor->has_periodic_neighbor(transverse_face)))
                        {
                          const auto diagonal_neighbor =
                            periodic_neighbor->neighbor_or_periodic_neighbor(
                              transverse_face);
                          const auto diagonal_entry = global_to_local.find(
                                                        diagonal_neighbor->global_active_cell_index());
                          AssertThrow(diagonal_entry != global_to_local.end(),
                                      dealii::ExcMessage(
                                        "FastScape periodic diagonal neighbor "
                                        "was not found."));

                          direction_type transverse_direction;
                          grid_data_type transverse_distance;
                          if (periodic_neighbor->has_periodic_neighbor(
                                transverse_face))
                            {
                              const unsigned int transverse_opposite_face =
                                periodic_neighbor->
                                periodic_neighbor_of_periodic_neighbor(
                                  transverse_face);
                              transverse_distance =
                                periodic_neighbor->center().distance(
                                  periodic_neighbor->face(
                                    transverse_face)->center()) +
                                diagonal_neighbor->center().distance(
                                  diagonal_neighbor->face(
                                    transverse_opposite_face)->center());
                              transverse_direction =
                                periodic_neighbor->face(
                                  transverse_face)->center() -
                                periodic_neighbor->center();
                              transverse_direction *= transverse_distance /
                                                      transverse_direction.norm();
                            }
                          else
                            {
                              transverse_direction =
                                diagonal_neighbor->center() -
                                periodic_neighbor->center();
                              transverse_distance =
                                transverse_direction.norm();
                            }

                          const direction_type diagonal_direction =
                            periodic_direction + transverse_direction;
                          add_neighbor(diagonal_entry->second,
                                       diagonal_direction.norm(),
                                       diagonal_direction,
                                       0.0);
                        }
                  }
          }

        neighbors_count[local_index] = count;
        ++local_index;
      }
  }


  template <class TriangulationType, class Selector>
  bool
  dealii_surface_grid<TriangulationType, Selector>::
  cell_touches_boundary(const size_type index,
                        const unsigned int coordinate_direction,
                        const bool upper_side) const
  {
    AssertIndexRange(index, m_size);
    AssertIndexRange(coordinate_direction, TriangulationType::dimension);
    return boundary_sides[index][2*coordinate_direction + (upper_side ? 1 : 0)];
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  number_of_cell_neighbors(const size_type index) const
  -> size_type
  {
    AssertIndexRange(index, m_size);
    return neighbors_count[index];
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  cell_neighbor(const size_type index,
                const size_type neighbor_number) const
  -> size_type
  {
    AssertIndexRange(index, m_size);
    AssertIndexRange(neighbor_number, neighbors_count[index]);
    return neighbors_indices[index][neighbor_number];
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  cell_neighbor_distance(const size_type index,
                         const size_type neighbor_number) const
  -> grid_data_type
  {
    AssertIndexRange(index, m_size);
    AssertIndexRange(neighbor_number, neighbors_count[index]);
    return neighbors_distances[index][neighbor_number];
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  cell_shared_face_measure(const size_type index,
                           const size_type neighbor_number) const
  -> grid_data_type
  {
    AssertIndexRange(index, m_size);
    AssertIndexRange(neighbor_number, neighbors_count[index]);
    return shared_face_measures[index][neighbor_number];
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  cell_neighbor_direction(const size_type index,
                          const size_type neighbor_number) const
  -> direction_type
  {
    AssertIndexRange(index, m_size);
    AssertIndexRange(neighbor_number, neighbors_count[index]);
    direction_type direction = neighbor_directions[index][neighbor_number];

    if (closed_surface)
      {
        const size_type neighbor_index =
          neighbors_indices[index][neighbor_number];
        direction_type surface_normal =
          cell_centers[index] + cell_centers[neighbor_index];
        const double normal_norm = surface_normal.norm();
        if (normal_norm > 0.0)
          {
            surface_normal /= normal_norm;
            direction -= (direction * surface_normal) * surface_normal;
          }
      }

    const double direction_norm = direction.norm();
    AssertThrow(direction_norm > 0.0,
                dealii::ExcMessage("Degenerate FastScape cell-neighbor direction."));
    return direction / direction_norm;
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::nodes_areas_impl() const
  -> container_type
  {
    return node_areas;
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  nodes_areas_impl(const size_type &index) const noexcept
  -> grid_data_type
  {
    return node_areas[index];
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  neighbors_count_impl(const size_type &index) const
  -> size_type
  {
    return neighbors_count[index];
  }


  template <class TriangulationType, class Selector>
  void
  dealii_surface_grid<TriangulationType, Selector>::
  neighbors_indices_impl(neighbors_indices_impl_type &output,
                         const size_type &index) const
  {
    for (size_type i = 0; i < neighbors_count[index]; ++i)
      output[i] = neighbors_indices[index][i];
  }


  template <class TriangulationType, class Selector>
  auto
  dealii_surface_grid<TriangulationType, Selector>::
  neighbors_distances_impl(const size_type &index) const
  -> const neighbors_distances_impl_type &
  {
    return neighbors_distances[index];
  }
}

#undef ASPECT_FASTSCAPELIB_LEGACY_GRID_API

#endif
#endif
