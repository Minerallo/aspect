/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/simulator/solver/stokes_matrix_free.h>
#include <aspect/simulator/solver/matrix_free_operators.h>
#include <aspect/simulator/solver/stokes_matrix_free_local_smoothing.h>
#include <aspect/simulator/solver/stokes_matrix_free_global_coarsening.h>

namespace aspect
{
  template <int dim>
  std::unique_ptr<StokesMatrixFreeHandler<dim>> create_matrix_free_solver(Simulator<dim> &simulator, const Parameters<dim> &parameters)
  {
    if (parameters.stokes_gmg_type == Parameters<dim>::StokesGMGType::local_smoothing)
      {
        switch (parameters.stokes_velocity_degree)
          {
            case 2:
              return std::make_unique<StokesMatrixFreeHandlerLocalSmoothingImplementation<dim,2>>(simulator, parameters);
              break;
            case 3:
              return std::make_unique<StokesMatrixFreeHandlerLocalSmoothingImplementation<dim,3>>(simulator, parameters);
              break;
            default:
              AssertThrow(false, ExcMessage("The finite element degree for the Stokes system you selected is not supported yet."));
          }
      }
    else if (parameters.stokes_gmg_type == Parameters<dim>::StokesGMGType::global_coarsening)
      {
        switch (parameters.stokes_velocity_degree)
          {
            case 2:
              return std::make_unique<StokesMatrixFreeHandlerGlobalCoarseningImplementation<dim,2>>(simulator, parameters);
              break;
            case 3:
              return std::make_unique<StokesMatrixFreeHandlerGlobalCoarseningImplementation<dim,3>>(simulator, parameters);
              break;
            default:
              AssertThrow(false, ExcMessage("The finite element degree for the Stokes system you selected is not supported yet."));
          }
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }

  template <int dim>
  void StokesMatrixFreeHandler<dim>::declare_parameters(ParameterHandler &prm)
  {
    StokesMatrixFreeHandlerLocalSmoothingImplementation<dim,2>::declare_parameters(prm);
    StokesMatrixFreeHandlerGlobalCoarseningImplementation<dim,2>::declare_parameters(prm);
  }

}

// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template std::unique_ptr<StokesMatrixFreeHandler<dim>> create_matrix_free_solver(Simulator<dim> &, const Parameters<dim> &); \
  template class StokesMatrixFreeHandler<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
