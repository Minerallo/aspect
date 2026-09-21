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

#include "common.h"

#include <aspect/config.h>

#ifdef ASPECT_WITH_FASTSCAPELIB

#  include <aspect/mesh_deformation/fastscape_cpp_grid.h>

#  include <deal.II/grid/grid_generator.h>
#  include <deal.II/grid/grid_tools.h>
#  include <deal.II/grid/tria.h>


TEST_CASE("Fastscapelib deal.II surface grid")
{
  dealii::Triangulation<2> triangulation;
  dealii::GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  fastscapelib::dealii_surface_grid<dealii::Triangulation<2>> grid(
    triangulation,
    false);

  REQUIRE(grid.size() == 4);

  for (unsigned int cell = 0; cell < grid.size(); ++cell)
    {
      REQUIRE(grid.number_of_cell_neighbors(cell) == 2);
      REQUIRE(grid.nodes_areas(cell) == Approx(0.25));

      for (unsigned int neighbor = 0;
           neighbor < grid.number_of_cell_neighbors(cell);
           ++neighbor)
        {
          REQUIRE(grid.cell_neighbor_distance(cell, neighbor) == Approx(0.5));
          REQUIRE(grid.cell_shared_face_measure(cell, neighbor) == Approx(0.5));
          REQUIRE(grid.cell_neighbor_direction(cell, neighbor).norm()
                  == Approx(1.0));
        }
    }
}


TEST_CASE("Fastscapelib deal.II D8 surface grid")
{
  dealii::Triangulation<2> triangulation;
  dealii::GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  fastscapelib::dealii_surface_grid<dealii::Triangulation<2>> grid(
    triangulation,
    false,
    true);

  REQUIRE(grid.size() == 4);
  for (unsigned int cell = 0; cell < grid.size(); ++cell)
    {
      REQUIRE(grid.number_of_cell_neighbors(cell) == 3);
      unsigned int diagonal_neighbors = 0;
      for (unsigned int neighbor = 0;
           neighbor < grid.number_of_cell_neighbors(cell);
           ++neighbor)
        {
          if (grid.cell_shared_face_measure(cell, neighbor) == Approx(0.0))
            {
              ++diagonal_neighbors;
              REQUIRE(grid.cell_neighbor_distance(cell, neighbor)
                      == Approx(std::sqrt(0.5)));
            }
          REQUIRE(grid.cell_neighbor_direction(cell, neighbor).norm()
                  == Approx(1.0));
        }
      REQUIRE(diagonal_neighbors == 1);
    }
}


TEST_CASE("Fastscapelib deal.II periodic surface grid")
{
  dealii::Triangulation<2> triangulation;
  dealii::GridGenerator::subdivided_hyper_rectangle(
    triangulation,
    std::vector<unsigned int>({4,3}),
    dealii::Point<2>(),
    dealii::Point<2>(1.0,1.0),
    true);

  std::vector<dealii::GridTools::PeriodicFacePair<
  dealii::Triangulation<2>::cell_iterator>> periodicity;
  dealii::GridTools::collect_periodic_faces(triangulation,
                                            2,
                                            3,
                                            1,
                                            periodicity);
  triangulation.add_periodicity(periodicity);

  fastscapelib::dealii_surface_grid<dealii::Triangulation<2>> grid(
    triangulation,
    false);

  REQUIRE(grid.size() == 12);
  for (unsigned int cell = 0; cell < grid.size(); ++cell)
    {
      REQUIRE_FALSE(grid.cell_touches_boundary(cell, 1, false));
      REQUIRE_FALSE(grid.cell_touches_boundary(cell, 1, true));
      REQUIRE(grid.number_of_cell_neighbors(cell) >= 3);
      REQUIRE(grid.number_of_cell_neighbors(cell) <= 4);

      for (unsigned int neighbor = 0;
           neighbor < grid.number_of_cell_neighbors(cell);
           ++neighbor)
        {
          REQUIRE(grid.cell_neighbor_distance(cell, neighbor) > 0.0);
          REQUIRE(grid.cell_neighbor_direction(cell, neighbor).norm()
                  == Approx(1.0));
        }
    }
}


TEST_CASE("Fastscapelib deal.II periodic D8 surface grid")
{
  dealii::Triangulation<2> triangulation;
  dealii::GridGenerator::subdivided_hyper_rectangle(
    triangulation,
    std::vector<unsigned int>({4,3}),
    dealii::Point<2>(),
    dealii::Point<2>(1.0,1.0),
    true);

  std::vector<dealii::GridTools::PeriodicFacePair<
  dealii::Triangulation<2>::cell_iterator>> periodicity;
  dealii::GridTools::collect_periodic_faces(triangulation,
                                            2,
                                            3,
                                            1,
                                            periodicity);
  triangulation.add_periodicity(periodicity);

  fastscapelib::dealii_surface_grid<dealii::Triangulation<2>> grid(
    triangulation,
    false,
    true);

  unsigned int cells_with_five_neighbors = 0;
  unsigned int cells_with_eight_neighbors = 0;
  for (unsigned int cell = 0; cell < grid.size(); ++cell)
    {
      REQUIRE_FALSE(grid.cell_touches_boundary(cell, 1, false));
      REQUIRE_FALSE(grid.cell_touches_boundary(cell, 1, true));
      if (grid.number_of_cell_neighbors(cell) == 5)
        ++cells_with_five_neighbors;
      else if (grid.number_of_cell_neighbors(cell) == 8)
        ++cells_with_eight_neighbors;
      else
        FAIL("Periodic D8 cell has an unexpected neighbor count");

      for (unsigned int neighbor = 0;
           neighbor < grid.number_of_cell_neighbors(cell);
           ++neighbor)
        {
          REQUIRE(grid.cell_neighbor_distance(cell, neighbor) > 0.0);
          REQUIRE(grid.cell_neighbor_direction(cell, neighbor).norm()
                  == Approx(1.0));
          for (unsigned int other = neighbor + 1;
               other < grid.number_of_cell_neighbors(cell);
               ++other)
            REQUIRE(grid.cell_neighbor(cell, neighbor) !=
                    grid.cell_neighbor(cell, other));
        }
    }
  REQUIRE(cells_with_five_neighbors == 6);
  REQUIRE(cells_with_eight_neighbors == 6);
}

#endif
