/*
 Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/global.h>
#include <aspect/simulator.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/volume_of_fluid/utilities.h>

#include <array>

namespace aspect
{
  template <>
  void VolumeOfFluidHandler<2>::update_volume_of_fluid_normals (const VolumeOfFluidField<2> &field,
                                                                LinearAlgebra::BlockVector &solution)
  {
    const unsigned int dim = 2;
    const unsigned int max_degree = 1;

    LinearAlgebra::BlockVector initial_solution;

    this->get_computing_timer().enter_subsection("Reconstruct VolumeOfFluid interfaces");

    initial_solution.reinit(sim.system_rhs, false);

    // Boundary reference
    typename DoFHandler<dim>::active_cell_iterator endc =
      this->get_dof_handler().end ();

    // Number of cells in the local reconstruction stencil
    const unsigned int n_cells_local_stencil = 9;

    Vector<double> local_volume_of_fluids (n_cells_local_stencil);
    std::vector<Point<dim>> stencil_unit_cell_centers (n_cells_local_stencil);
    std::vector<typename DoFHandler<dim>::active_cell_iterator> neighbor_cells(n_cells_local_stencil);

    const unsigned int n_candidate_normals_per_dim = 3; // Named variable for number of candidate sums for each dimension
    std::vector<double> strip_sums (dim * n_candidate_normals_per_dim);

    // Named variable for number of candidate interface normal vectors for the reconstruction
    const unsigned int n_candidate_normals = dim*n_candidate_normals_per_dim+1;
    std::vector<Tensor<1, dim, double>> normals (n_candidate_normals);
    std::vector<double> d_vals (n_candidate_normals);
    std::vector<double> errs (n_candidate_normals);

    // Variables to do volume calculations

    const QGauss<dim> quadrature(max_degree);

    std::vector<double> xFEM_values(quadrature.size());

    const FiniteElement<dim> &system_fe = this->get_fe();

    FEValues<dim> fevalues(this->get_mapping(), system_fe, quadrature,
                           update_JxW_values);

    // Normal holding vars
    Point<dim> reconstruction_stencil_unit_cell_center;
    Tensor<1, dim, double> normal;
    double d;

    for (unsigned int i=0; i<dim; ++i)
      reconstruction_stencil_unit_cell_center[i] = 0.5;

    std::vector<types::global_dof_index> cell_dof_indices (system_fe.dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (system_fe.dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = field.volume_fraction;
    const unsigned int volume_of_fluid_c_index = volume_of_fluid_var.first_component_index;
    const unsigned int volume_of_fluid_ind
      = this->get_fe().component_to_system_index(volume_of_fluid_c_index, 0);

    const FEVariable<dim> &volume_of_fluidN_var = field.reconstruction;
    const unsigned int volume_of_fluidN_c_index = volume_of_fluidN_var.first_component_index;
    const unsigned int volume_of_fluidN_blockidx = volume_of_fluidN_var.block_index;

    const FEVariable<dim> &volume_of_fluidLS_var = field.level_set;
    const unsigned int volume_of_fluidLS_c_index = volume_of_fluidLS_var.first_component_index;
    const unsigned int n_volume_of_fluidLS_dofs = volume_of_fluidLS_var.fe->dofs_per_cell;
    const unsigned int volume_of_fluidLS_blockidx = volume_of_fluidLS_var.block_index;

    //Iterate over cells
    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        // Obtain data for this cell and neighbors
        cell->get_dof_indices (local_dof_indices);
        const double cell_volume_of_fluid = solution(local_dof_indices[volume_of_fluid_ind]);

        normal[0] = 0.0;
        normal[1] = 0.0;
        d = -1.0;

        if (cell_volume_of_fluid > 1.0 - volume_fraction_threshold)
          {
            d = 1.0;
            initial_solution(local_dof_indices[volume_of_fluid_ind]) = 1.0;
          }
        else if (cell_volume_of_fluid < volume_fraction_threshold)
          {
            d = -1.0;
            initial_solution(local_dof_indices[volume_of_fluid_ind]) = 0.0;
          }
        else
          {
            initial_solution(local_dof_indices[volume_of_fluid_ind]) = cell_volume_of_fluid;

            // Get references to neighboring cells to build stencil references
            //
            // Due to mesh structure, we need to obtain the cells by
            // considering the pattern of neighboring cells
            //
            // Indices used for the stencil references follow the
            // pattern
            //
            // 6 7 8
            // 3 4 5
            // 0 1 2
            const unsigned int stencil_side_cell_count = 3;
            for (unsigned int i = 0; i < stencil_side_cell_count; ++i)
              {
                typename DoFHandler<dim>::active_cell_iterator cen; // holding variable for cell at column center
                if (i == 0 || i == 2)
                  {
                    // Not on center column, so obtain the center cell of the appropriate row
                    const unsigned int neighbor_no = (i/2);
                    const typename DoFHandler<dim>::face_iterator face = cell->face (neighbor_no);
                    if ((face->at_boundary() && !cell->has_periodic_neighbor(neighbor_no)) ||
                        face->has_children())
                      cen = endc;
                    else
                      {
                        const typename DoFHandler<dim>::cell_iterator neighbor =
                          cell->neighbor_or_periodic_neighbor(neighbor_no);
                        if (neighbor->level() == cell->level() &&
                            neighbor->is_active())
                          cen = neighbor;
                        else
                          cen = endc;
                      }
                  }
                else
                  {
                    // On center column, so current cell is center of column
                    cen = cell;
                  }

                for (unsigned int j = 0; j < stencil_side_cell_count; ++j)
                  {
                    typename DoFHandler<dim>::active_cell_iterator curr; // Variable for cell at current stencil location
                    if (cen == endc)
                      {
                        // Current column center is not in mesh, so assume that
                        // the desired piece of the stencil is also outside.
                        curr = endc;
                      }
                    else
                      {
                        // Current column center exists, so obtain correct cell reference
                        if (j == 0 || j == 2)
                          {
                            //
                            const unsigned int neighbor_no = 2+(j/2);
                            const typename DoFHandler<dim>::face_iterator face = cen->face (neighbor_no);
                            if ((face->at_boundary() && !cen->has_periodic_neighbor(neighbor_no)) ||
                                face->has_children())
                              curr = endc;
                            else
                              {
                                const typename DoFHandler<dim>::cell_iterator neighbor =
                                  cen->neighbor_or_periodic_neighbor(neighbor_no);
                                if (neighbor->level() == cell->level() &&
                                    neighbor->is_active())
                                  curr = neighbor;
                                else
                                  curr = endc;
                              }
                          }
                        else
                          {
                            // Current stencil reference is column center column center
                            curr = cen;
                          }
                      }
                    if (curr != endc)
                      {
                        // Cell reference is valid, so get data
                        curr->get_dof_indices (cell_dof_indices);
                        stencil_unit_cell_centers[3 * j + i] = Point<dim> (-1.0 + i,
                                                                           -1.0 + j);
                      }
                    else
                      {
                        // Cell reference is invalid, so replace with current
                        // cell to reduce branching complexity in the later
                        // algorithm
                        cell->get_dof_indices (cell_dof_indices);
                        stencil_unit_cell_centers[3 * j + i] = Point<dim> (0.0,
                                                                           0.0);
                      }
                    local_volume_of_fluids (3 * j + i) = solution (cell_dof_indices[volume_of_fluid_ind]);
                    neighbor_cells[3 * j + i] = curr;
                  }
              }

            // Gather cell strip sums
            //
            // Sums are indexed as
            // [n_sums_per_dim * parallel_dim + ind]
            const unsigned int n_sums_per_dim = 3;
            for (unsigned int i = 0; i < dim * n_candidate_normals_per_dim; ++i)
              strip_sums[i] = 0.0;

            for (unsigned int i = 0; i < stencil_side_cell_count; ++i)
              {
                for (unsigned int j = 0; j < stencil_side_cell_count; ++j)
                  {
                    strip_sums[n_sums_per_dim * 0 + i] += local_volume_of_fluids (stencil_side_cell_count * j + i);
                    strip_sums[n_sums_per_dim * 1 + j] += local_volume_of_fluids (stencil_side_cell_count * j + i);
                  }
              }

            // Calculate normal vectors for the 6 candidates from the efficient
            // least squares approach
            //
            // Labeling the sums of the perpendicular strips as
            // 0 1 2
            // L C R
            //
            // For each dimension consider the interface normals implied by the
            // 3 divided differences
            //
            // L-C/1
            // C-R/2
            // L-R/1
            //
            for (unsigned int di = 0; di < dim; ++di)
              {
                // Get index other dimension
                unsigned int di2 = (di == 0) ? 1 : 0;
                for (unsigned int i = 0; i < 3; ++i)
                  {
                    normals[n_candidate_normals_per_dim * di + i][di] = 0.0;
                    normals[n_candidate_normals_per_dim * di + i][di2] = 0.0;
                    if (i==0 || i == 2)
                      {
                        // Positive sum in difference is L
                        normals[n_candidate_normals_per_dim * di + i][di] += strip_sums[n_sums_per_dim * di + 0];
                        normals[3 * di + i][di2] += 1.0;
                      }
                    else
                      {
                        // Positive sum in difference is C
                        normals[3 * di + i][di] += strip_sums[3 * di + 1];
                        normals[3 * di + i][di2] += 0.0;
                      }
                    if (i == 0)
                      {
                        // Negative sum in difference is C
                        normals[3 * di + i][di] -= strip_sums[3 * di + 1];
                        normals[3 * di + i][di2] += 0.0;
                      }
                    else
                      {
                        // Negative sum in difference is R
                        normals[3 * di + i][di] -= strip_sums[3 * di + 2];
                        normals[3 * di + i][di2] += 1.0;
                      }

                    if (strip_sums[3 * di2 + 2] > strip_sums[3 * di2 + 0])
                      {
                        // There is more fluid in in area above the interface on the stencil,
                        // so flip normal direction
                        normals[3 * di + i][di2] *= -1.0;
                      }
                  }
              }

            // Add time extrapolated local normal as candidate
            // this is not expected to be the best candidate in general, but
            // should result in exact reconstruction for linear interface
            // translations
            // Inclusion of this candidate will not reduce accuracy due to it
            // only being selected if it produces a better interface
            // approximation than the ELS candidates. Note that this will
            // render linear translation problems less dependent on the
            // interface reconstruction, so other tests will also be necessary.
            for (unsigned int i=0; i<dim; ++i)
              normals[6][i] = solution(local_dof_indices[system_fe
                                                         .component_to_system_index(volume_of_fluidN_c_index+i, 0)]);

            // If candidate normal too small, remove from consideration
            if (normals[6]*normals[6]< volume_fraction_threshold)
              {
                normals[6][0] = 0;
                normals[6][1] = 0;
              }

            unsigned int index_of_best_normal = 0;
            {
              fevalues.reinit(cell);
              const std::vector<double> weights = fevalues.get_JxW_values();

              double cell_vol = 0.0;
              for (const double weight : weights)
                {
                  cell_vol+=weight;
                }
              for (unsigned int nind = 0; nind < n_candidate_normals; ++nind)
                {
                  errs[nind] = 0.0;
                  const double normal_norm = normals[nind].norm_square();

                  if (normal_norm > volume_fraction_threshold) // If candidate normal too small set error to maximum
                    {
                      d_vals[nind] = VolumeOfFluid::Utilities::compute_interface_location_newton<dim> (
                                       max_degree,
                                       normals[nind],
                                       cell_volume_of_fluid,
                                       cell_vol,
                                       volume_of_fluid_reconstruct_epsilon,
                                       quadrature.get_points(), weights);
                    }
                  else
                    {
                      errs[nind] = 9.0;
                    }
                }
            }

            for (unsigned int i = 0; i < n_cells_local_stencil; ++i)
              {
                if (neighbor_cells[i] == endc)
                  {
                    continue;
                  }

                fevalues.reinit(neighbor_cells[i]);

                const std::vector<double> weights = fevalues.get_JxW_values();

                double cell_vol = 0.0;
                for (const double weight : weights)
                  {
                    cell_vol+=weight;
                  }

                for (unsigned int nind = 0; nind < n_candidate_normals; ++nind)
                  {
                    const double normal_norm = normals[nind]*normals[nind];

                    if (normal_norm > volume_fraction_threshold) // If candidate normal too small skip as set to max already
                      {
                        double dot = 0.0;
                        for (unsigned int di = 0; di < dim; ++di)
                          dot += normals[nind][di] * stencil_unit_cell_centers[i][di];
                        const double n_volume_of_fluid = VolumeOfFluid::Utilities::compute_fluid_volume<dim> (max_degree, normals[nind], d_vals[nind]-dot,
                                                         quadrature.get_points(), weights)/cell_vol;
                        const double cell_err = local_volume_of_fluids (i) - n_volume_of_fluid;
                        errs[nind] += cell_err * cell_err;
                      }
                  }
              }

            for (unsigned int nind = 0; nind < n_candidate_normals; ++nind)
              {
                if (errs[index_of_best_normal] >= errs[nind])
                  index_of_best_normal = nind;
              }

            normal = normals[index_of_best_normal];
            d = d_vals[index_of_best_normal];
          }

        for (unsigned int i=0; i<dim; ++i)
          initial_solution (local_dof_indices[system_fe
                                              .component_to_system_index(volume_of_fluidN_c_index+i, 0)]) = normal[i];

        initial_solution (local_dof_indices[system_fe
                                            .component_to_system_index(volume_of_fluidN_c_index+dim, 0)]) = d;

        for (unsigned int i=0; i<n_volume_of_fluidLS_dofs; ++i)
          {
            // Recenter unit cell on origin
            Tensor<1, dim, double> recentered_support_point = volume_of_fluidLS_var.fe->unit_support_point(i)-reconstruction_stencil_unit_cell_center;
            initial_solution (local_dof_indices[system_fe
                                                .component_to_system_index(volume_of_fluidLS_c_index, i)])
              = d-recentered_support_point*normal;
          }
      }

    initial_solution.compress(VectorOperation::insert);

    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    solution.block(volume_of_fluidN_blockidx) = initial_solution.block(volume_of_fluidN_blockidx);
    solution.block(volume_of_fluidLS_blockidx) = initial_solution.block(volume_of_fluidLS_blockidx);

    this->get_computing_timer().leave_subsection("Reconstruct VolumeOfFluid interfaces");
  }


  template <>
  void VolumeOfFluidHandler<3>::update_volume_of_fluid_normals (const VolumeOfFluidField<3> &field,
                                                                LinearAlgebra::BlockVector &solution)
  {
    const unsigned int dim = 3;
    const unsigned int stencil_side_cell_count = 3;
    const unsigned int n_cells_local_stencil = 27;
    const unsigned int n_candidate_normals = 28;

    LinearAlgebra::BlockVector initial_solution;
    initial_solution.reinit(sim.system_rhs, false);

    this->get_computing_timer().enter_subsection("Reconstruct VolumeOfFluid interfaces");

    const typename DoFHandler<dim>::active_cell_iterator endc =
      this->get_dof_handler().end();

    Vector<double> local_volume_of_fluids(n_cells_local_stencil);
    std::vector<Point<dim>> stencil_unit_cell_centers(n_cells_local_stencil);
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    neighbor_cells(n_cells_local_stencil);
    std::vector<Tensor<1, dim>> normals(n_candidate_normals);
    std::vector<double> d_values(n_candidate_normals);
    std::vector<double> errors(n_candidate_normals);

    Point<dim> reconstruction_stencil_unit_cell_center;
    for (unsigned int direction = 0; direction < dim; ++direction)
      reconstruction_stencil_unit_cell_center[direction] = 0.5;

    const FiniteElement<dim> &system_fe = this->get_fe();
    std::vector<types::global_dof_index> cell_dof_indices(system_fe.dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(system_fe.dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = field.volume_fraction;
    const unsigned int volume_of_fluid_index =
      system_fe.component_to_system_index(volume_of_fluid_var.first_component_index, 0);

    const FEVariable<dim> &reconstruction_var = field.reconstruction;
    const unsigned int reconstruction_component = reconstruction_var.first_component_index;
    const unsigned int reconstruction_block = reconstruction_var.block_index;

    const FEVariable<dim> &level_set_var = field.level_set;
    const unsigned int level_set_component = level_set_var.first_component_index;
    const unsigned int n_level_set_dofs = level_set_var.fe->dofs_per_cell;
    const unsigned int level_set_block = level_set_var.block_index;

    const auto stencil_index = [] (const unsigned int x,
                                   const unsigned int y,
                                   const unsigned int z)
    {
      return x + stencil_side_cell_count * (y + stencil_side_cell_count * z);
    };

    for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        cell->get_dof_indices(local_dof_indices);
        const double cell_volume_of_fluid = solution(local_dof_indices[volume_of_fluid_index]);

        Tensor<1, dim> normal;
        double interface_location = -1.0;

        if (cell_volume_of_fluid > 1.0 - volume_fraction_threshold)
          interface_location = 1.0;
        else if (cell_volume_of_fluid >= volume_fraction_threshold)
          {
            for (unsigned int z = 0; z < stencil_side_cell_count; ++z)
              for (unsigned int y = 0; y < stencil_side_cell_count; ++y)
                for (unsigned int x = 0; x < stencil_side_cell_count; ++x)
                  {
                    const std::array<unsigned int, dim> stencil_coordinates = {{x, y, z}};
                    typename DoFHandler<dim>::active_cell_iterator current_cell = cell;
                    bool valid_neighbor = true;

                    for (unsigned int direction = 0; direction < dim; ++direction)
                      if (stencil_coordinates[direction] != 1)
                        {
                          const unsigned int face_number =
                            2 * direction + (stencil_coordinates[direction] == 2 ? 1 : 0);
                          const typename DoFHandler<dim>::face_iterator face =
                            current_cell->face(face_number);

                          if ((face->at_boundary()
                               && !current_cell->has_periodic_neighbor(face_number))
                              || face->has_children())
                            {
                              valid_neighbor = false;
                              break;
                            }

                          const typename DoFHandler<dim>::cell_iterator neighbor =
                            current_cell->neighbor_or_periodic_neighbor(face_number);
                          if (neighbor->level() != cell->level() || !neighbor->is_active())
                            {
                              valid_neighbor = false;
                              break;
                            }
                          current_cell = neighbor;
                        }

                    const unsigned int index = stencil_index(x, y, z);
                    if (valid_neighbor)
                      {
                        current_cell->get_dof_indices(cell_dof_indices);
                        stencil_unit_cell_centers[index] =
                          Point<dim>(static_cast<double>(x) - 1.0,
                                     static_cast<double>(y) - 1.0,
                                     static_cast<double>(z) - 1.0);
                        neighbor_cells[index] = current_cell;
                      }
                    else
                      {
                        cell->get_dof_indices(cell_dof_indices);
                        stencil_unit_cell_centers[index] = Point<dim>();
                        neighbor_cells[index] = endc;
                      }

                    local_volume_of_fluids(index) =
                      solution(cell_dof_indices[volume_of_fluid_index]);
                  }

            unsigned int candidate_index = 0;
            for (unsigned int height_direction = 0;
                 height_direction < dim;
                 ++height_direction)
              {
                std::array<unsigned int, 2> tangent_directions;
                unsigned int tangent_index = 0;
                for (unsigned int direction = 0; direction < dim; ++direction)
                  if (direction != height_direction)
                    tangent_directions[tangent_index++] = direction;

                std::array<std::array<double, 3>, 3> column_sums = {{}};
                std::array<double, 3> tangent_profile_0 = {{0.0, 0.0, 0.0}};
                std::array<double, 3> tangent_profile_1 = {{0.0, 0.0, 0.0}};
                double negative_height_plane_sum = 0.0;
                double positive_height_plane_sum = 0.0;

                for (unsigned int tangent_1 = 0; tangent_1 < 3; ++tangent_1)
                  for (unsigned int tangent_0 = 0; tangent_0 < 3; ++tangent_0)
                    {
                      for (unsigned int height = 0; height < 3; ++height)
                        {
                          std::array<unsigned int, dim> coordinates = {{0, 0, 0}};
                          coordinates[height_direction] = height;
                          coordinates[tangent_directions[0]] = tangent_0;
                          coordinates[tangent_directions[1]] = tangent_1;
                          const double fluid = local_volume_of_fluids(
                                                 stencil_index(coordinates[0],
                                                               coordinates[1],
                                                               coordinates[2]));
                          column_sums[tangent_0][tangent_1] += fluid;
                          if (height == 0)
                            negative_height_plane_sum += fluid;
                          else if (height == 2)
                            positive_height_plane_sum += fluid;
                        }
                    }

                for (unsigned int position = 0; position < 3; ++position)
                  for (unsigned int transverse_position = 0;
                       transverse_position < 3;
                       ++transverse_position)
                    {
                      tangent_profile_0[position] +=
                        column_sums[position][transverse_position] / 3.0;
                      tangent_profile_1[position] +=
                        column_sums[transverse_position][position] / 3.0;
                    }

                const std::array<double, 3> tangent_differences_0 =
                {{tangent_profile_0[0] - tangent_profile_0[1],
                  tangent_profile_0[1] - tangent_profile_0[2],
                  0.5 * (tangent_profile_0[0] - tangent_profile_0[2])}};
                const std::array<double, 3> tangent_differences_1 =
                {{tangent_profile_1[0] - tangent_profile_1[1],
                  tangent_profile_1[1] - tangent_profile_1[2],
                  0.5 * (tangent_profile_1[0] - tangent_profile_1[2])}};
                const double height_component =
                  (positive_height_plane_sum > negative_height_plane_sum ? -1.0 : 1.0);

                for (unsigned int difference_1 = 0; difference_1 < 3; ++difference_1)
                  for (unsigned int difference_0 = 0; difference_0 < 3; ++difference_0)
                    {
                      normals[candidate_index][height_direction] = height_component;
                      normals[candidate_index][tangent_directions[0]] =
                        tangent_differences_0[difference_0];
                      normals[candidate_index][tangent_directions[1]] =
                        tangent_differences_1[difference_1];
                      ++candidate_index;
                    }
              }

            Assert(candidate_index == n_candidate_normals - 1, ExcInternalError());

            for (unsigned int direction = 0; direction < dim; ++direction)
              normals.back()[direction] =
                solution(local_dof_indices[system_fe.component_to_system_index(
                           reconstruction_component + direction, 0)]);

            for (unsigned int candidate = 0; candidate < n_candidate_normals; ++candidate)
              {
                errors[candidate] = 0.0;
                if (normals[candidate].norm_square() > volume_fraction_threshold)
                  d_values[candidate] =
                    VolumeOfFluid::Utilities::compute_interface_location(normals[candidate],
                                                                         cell_volume_of_fluid);
                else
                  errors[candidate] = static_cast<double>(n_cells_local_stencil);
              }

            for (unsigned int neighbor_index = 0;
                 neighbor_index < n_cells_local_stencil;
                 ++neighbor_index)
              {
                if (neighbor_cells[neighbor_index] == endc)
                  continue;

                for (unsigned int candidate = 0;
                     candidate < n_candidate_normals;
                     ++candidate)
                  if (normals[candidate].norm_square() > volume_fraction_threshold)
                    {
                      const double offset =
                        normals[candidate] * stencil_unit_cell_centers[neighbor_index];
                      const double reconstructed_fraction =
                        VolumeOfFluid::Utilities::compute_fluid_fraction(
                          normals[candidate], d_values[candidate] - offset);
                      const double cell_error =
                        local_volume_of_fluids(neighbor_index) - reconstructed_fraction;
                      errors[candidate] += cell_error * cell_error;
                    }
              }

            const unsigned int best_candidate =
              std::distance(errors.begin(), std::min_element(errors.begin(), errors.end()));
            normal = normals[best_candidate];
            interface_location = d_values[best_candidate];
          }

        for (unsigned int direction = 0; direction < dim; ++direction)
          initial_solution(local_dof_indices[system_fe.component_to_system_index(
                             reconstruction_component + direction, 0)]) = normal[direction];
        initial_solution(local_dof_indices[system_fe.component_to_system_index(
                           reconstruction_component + dim, 0)]) = interface_location;

        for (unsigned int i = 0; i < n_level_set_dofs; ++i)
          {
            const Tensor<1, dim> recentered_support_point =
              level_set_var.fe->unit_support_point(i)
              - reconstruction_stencil_unit_cell_center;
            initial_solution(local_dof_indices[system_fe.component_to_system_index(
                               level_set_component, i)]) =
              interface_location - recentered_support_point * normal;
          }
      }

    initial_solution.compress(VectorOperation::insert);
    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    solution.block(reconstruction_block) = initial_solution.block(reconstruction_block);
    solution.block(level_set_block) = initial_solution.block(level_set_block);

    this->get_computing_timer().leave_subsection("Reconstruct VolumeOfFluid interfaces");
  }

  template <>
  void VolumeOfFluidHandler<2>::update_volume_of_fluid_composition (const typename Simulator<2>::AdvectionField &composition_field,
                                                                    const VolumeOfFluidField<2> &volume_of_fluid_field,
                                                                    LinearAlgebra::BlockVector &solution)
  {
    const unsigned int dim = 2;

    LinearAlgebra::BlockVector initial_solution;

    this->get_computing_timer().enter_subsection("Compute VolumeOfFluid compositions");

    initial_solution.reinit(sim.system_rhs, false);

    // Normal holding vars
    Point<dim> reconstruction_stencil_unit_cell_center;

    for (unsigned int i=0; i<dim; ++i)
      reconstruction_stencil_unit_cell_center[i] = 0.5;

    const FiniteElement<dim> &system_fe = this->get_fe();

    std::vector<types::global_dof_index> local_dof_indices (system_fe.dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = volume_of_fluid_field.volume_fraction;
    const unsigned int volume_of_fluid_c_index = volume_of_fluid_var.first_component_index;
    const unsigned int volume_of_fluid_ind
      = system_fe.component_to_system_index(volume_of_fluid_c_index, 0);

    const FEVariable<dim> &volume_of_fluidN_var = volume_of_fluid_field.reconstruction;
    const unsigned int volume_of_fluidN_c_index = volume_of_fluidN_var.first_component_index;

    const unsigned int base_element = composition_field.base_element(this->introspection());
    const std::vector<Point<dim>> support_points = system_fe.base_element(base_element).get_unit_support_points();

    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned())
          continue;

        cell->get_dof_indices (local_dof_indices);
        const double cell_volume_of_fluid = solution(local_dof_indices[volume_of_fluid_ind]);

        Tensor<1, dim, double> normal;

        for (unsigned int i=0; i<dim; ++i)
          normal[i] = solution(local_dof_indices[system_fe
                                                 .component_to_system_index(volume_of_fluidN_c_index+i, 0)]);

        // Compute L1 normalized normal vector
        // (Use of L1 normalization makes calculating the correct correction
        // factor much simpler to retain $0\leq C\leq 1$ bound much simpler)
        double normal_l1_norm = 0.0;
        for (unsigned int i=0; i<dim; ++i)
          {
            normal_l1_norm += std::abs(normal[i]);
          }
        //Calculate correct factor to retain vol frac and [0,1] bound
        double correction_factor = (normal_l1_norm<volume_fraction_threshold)
                                   ?
                                   0.0
                                   :
                                   2.0*(0.5-std::abs(cell_volume_of_fluid-0.5))/normal_l1_norm;
        for (unsigned int i=0; i<system_fe.base_element(base_element).dofs_per_cell; ++i)
          {
            const unsigned int system_local_dof
              = system_fe.component_to_system_index(composition_field.component_index(sim.introspection),
                                                    /*dof index within component*/i);

            Tensor<1, dim, double> recentered_support_point = support_points[i]-reconstruction_stencil_unit_cell_center;

            const double value = cell_volume_of_fluid - correction_factor*(recentered_support_point*normal);

            initial_solution(local_dof_indices[system_local_dof]) = value;
          }

      }

    initial_solution.compress(VectorOperation::insert);

    const unsigned int blockidx = composition_field.block_index(this->introspection());
    solution.block(blockidx) = initial_solution.block(blockidx);

    this->get_computing_timer().leave_subsection("Compute VolumeOfFluid compositions");
  }



  template <>
  void VolumeOfFluidHandler<3>::update_volume_of_fluid_composition (const Simulator<3>::AdvectionField &composition_field,
                                                                    const VolumeOfFluidField<3> &volume_of_fluid_field,
                                                                    LinearAlgebra::BlockVector &solution)
  {
    const unsigned int dim = 3;

    LinearAlgebra::BlockVector initial_solution;
    this->get_computing_timer().enter_subsection("Compute VolumeOfFluid compositions");
    initial_solution.reinit(sim.system_rhs, false);

    Point<dim> reconstruction_stencil_unit_cell_center;
    for (unsigned int direction = 0; direction < dim; ++direction)
      reconstruction_stencil_unit_cell_center[direction] = 0.5;

    const FiniteElement<dim> &system_fe = this->get_fe();
    std::vector<types::global_dof_index> local_dof_indices(system_fe.dofs_per_cell);

    const unsigned int volume_of_fluid_index =
      system_fe.component_to_system_index(
        volume_of_fluid_field.volume_fraction.first_component_index, 0);
    const unsigned int reconstruction_component =
      volume_of_fluid_field.reconstruction.first_component_index;

    const unsigned int base_element = composition_field.base_element(this->introspection());
    const std::vector<Point<dim>> support_points =
      system_fe.base_element(base_element).get_unit_support_points();

    for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        cell->get_dof_indices(local_dof_indices);
        const double cell_volume_of_fluid = solution(local_dof_indices[volume_of_fluid_index]);

        Tensor<1, dim> normal;
        double normal_l1_norm = 0.0;
        for (unsigned int direction = 0; direction < dim; ++direction)
          {
            normal[direction] =
              solution(local_dof_indices[system_fe.component_to_system_index(
                         reconstruction_component + direction, 0)]);
            normal_l1_norm += std::abs(normal[direction]);
          }

        const double correction_factor =
          (normal_l1_norm < volume_fraction_threshold
           ? 0.0
           : 2.0 * (0.5 - std::abs(cell_volume_of_fluid - 0.5)) / normal_l1_norm);

        for (unsigned int i = 0;
             i < system_fe.base_element(base_element).dofs_per_cell;
             ++i)
          {
            const unsigned int system_local_dof =
              system_fe.component_to_system_index(
                composition_field.component_index(sim.introspection), i);
            const Tensor<1, dim> recentered_support_point =
              support_points[i] - reconstruction_stencil_unit_cell_center;
            initial_solution(local_dof_indices[system_local_dof]) =
              cell_volume_of_fluid
              - correction_factor * (recentered_support_point * normal);
          }
      }

    initial_solution.compress(VectorOperation::insert);
    const unsigned int block_index = composition_field.block_index(this->introspection());
    solution.block(block_index) = initial_solution.block(block_index);
    this->get_computing_timer().leave_subsection("Compute VolumeOfFluid compositions");
  }
}
