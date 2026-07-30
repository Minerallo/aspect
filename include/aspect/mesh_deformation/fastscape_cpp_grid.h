/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#ifndef _aspect_mesh_deformation_fastscape_cpp_grid_h
#define _aspect_mesh_deformation_fastscape_cpp_grid_h

#include <aspect/config.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#include <deal.II/grid/tria.h>
#include <deal.II/base/tensor.h>

#include <fastscapelib/grid/base.hpp>
#include <fastscapelib/version.hpp>

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
                          const bool closed_surface);

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
                      const bool closed_surface)
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
    shared_face_measures(m_size)
  {
    m_shape = {{static_cast<typename shape_type::value_type>(m_size)}};
    std::fill(m_nodes_status.begin(), m_nodes_status.end(), node_status::core);

    std::unordered_map<unsigned int, size_type> global_to_local;
    size_type local_index = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        global_to_local[cell->global_active_cell_index()] = local_index;
        node_areas[local_index] = cell->measure();
        cell_centers[local_index] = cell->center();

        if (!closed_surface)
          for (const unsigned int face : cell->face_indices())
            if (cell->at_boundary(face))
              {
                m_nodes_status[local_index] = node_status::fixed_value;
                break;
              }
        ++local_index;
      }

    local_index = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        const dealii::Point<TriangulationType::space_dimension> center = cell->center();
        size_type count = 0;

        for (const unsigned int face : cell->face_indices())
          if (!cell->at_boundary(face))
            {
              const auto neighbor = cell->neighbor(face);
              const auto entry = global_to_local.find(neighbor->global_active_cell_index());
              AssertThrow(entry != global_to_local.end(),
                          dealii::ExcMessage("FastScape surface-grid neighbor was not found."));
              AssertIndexRange(count, inner_types::n_neighbors_max);
              neighbors_indices[local_index][count] = entry->second;
              neighbors_distances[local_index][count] = center.distance(neighbor->center());
              if constexpr (TriangulationType::dimension == 1)
                shared_face_measures[local_index][count] = 1.0;
              else
                shared_face_measures[local_index][count] =
                  cell->face(face)->measure();
              ++count;
            }

        neighbors_count[local_index] = count;
        ++local_index;
      }
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
    const size_type neighbor_index =
      neighbors_indices[index][neighbor_number];
    direction_type direction =
      cell_centers[neighbor_index] - cell_centers[index];

    if (closed_surface)
      {
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
