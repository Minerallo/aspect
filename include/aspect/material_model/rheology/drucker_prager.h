/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_drucker_prager_h
#define _aspect_material_model_rheology_drucker_prager_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>

#include <aspect/material_model/utilities.h>
#include <aspect/simulator_access.h>
namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * A data structure with all input file parameters relevant to
       * Drucker-Prager plasticity.
       */
      struct DruckerPragerParameters
      {
        /**
         * Internal friction angle (phi) for the current composition and phase
         */
        double angle_internal_friction;

        /**
         * Dilation angle (psi) for the current composition and phase
         */
        double angle_dilation;

        /**
         * Cohesion for the current composition and phase
         */
        double cohesion;

        /**
         * Prefactor for the yield stress for the current composition and phase
         */
        double yield_stress_prefactor;

        /**
         * Limit maximum yield stress from drucker prager yield criterion.
         */
        double max_yield_stress;

        /**
         * Constructor. Initializes all values to NaN.
         */
        DruckerPragerParameters();
      };

      template <int dim>
      class DruckerPrager : public ::aspect::SimulatorAccess<dim>
      {
        public:
          DruckerPrager();
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers, this is considered the number of phases
           * for each compositional field (plus possibly a background field)
           * and this number will be checked against the parsed parameters.
           *
           * @param [in] prm The ParameterHandler to read from.
           * @param expected_n_phases_per_composition Optional list of number of phases.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition = nullptr);

          /**
           * Compute the parameters for the Drucker Prager plasticity.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          const DruckerPragerParameters
          compute_drucker_prager_parameters (const unsigned int composition,
                                             const std::vector<double> &phase_function_values = std::vector<double>(),
                                             const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the plastic yield stress based on the Drucker Prager yield criterion.
           *
           * @deprecated This function is deprecated, use the other function with the same name
           * that takes DruckerPragerParameters as input.
           */
          DEAL_II_DEPRECATED
          double
          compute_yield_stress (const double cohesion,
                                const double angle_internal_friction,
                                const double pressure,
                                const double max_yield_stress) const;

          /**
           * Compute the plastic yield stress based on the Drucker Prager yield criterion for
           * the given pressure, and drucker prager parameters.
           *
           * @param pressure The current pressure.
           * @param p The parameters of the drucker prager plasticity. This can either be
           * manually created, or computed by compute_drucker_prager_parameters().
           *
           * @return The yield stress for the given conditions.
           */
          double
          compute_yield_stress (const double pressure,
                                const DruckerPragerParameters &p) const;

          /**
           * Compute the apparent viscosity using the yield stress and effective strain rate.
           * If the non_yielding_viscosity is not infinite
           * (i.e., if there are other rheological elements accommodating strain), the returned
           * value is the effective composite viscosity, not the pure "plastic" viscosity.
           *
           * @deprecated This function is deprecated, use the other function with the same name
           * that takes DruckerPragerParameters as input.
           */
          DEAL_II_DEPRECATED
          double
          compute_viscosity (const double cohesion,
                             const double angle_internal_friction,
                             const double pressure,
                             const double effective_strain_rate,
                             const double max_yield_stress,
                             const double non_yielding_viscosity = std::numeric_limits<double>::infinity()) const;

          /**
           * Compute the apparent viscosity using the yield stress and effective strain rate.
           * If the non_yielding_viscosity is not infinite
           * (i.e., if there are other rheological elements accommodating strain), the returned
           * value is the effective composite viscosity, not the pure "plastic" viscosity.
           */
          double
          compute_viscosity (const double pressure,
                             const double effective_strain_rate,
                             const DruckerPragerParameters &p,
                             const double non_yielding_viscosity = std::numeric_limits<double>::infinity()) const;

          /**
           * Compute the strain rate and first stress derivative
           * as a function of stress based on the damped Drucker-Prager flow law.
           */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double stress,
                                              const double pressure,
                                              const DruckerPragerParameters &p) const;

          /**
           * Compute the derivative of the plastic viscosity with respect to pressure.
           */
          double
          compute_derivative (const double angle_internal_friction,
                              const double effective_strain_rate) const;

          /**
           * Compute the LHS and RHS dilation terms for the Stokes system.
           * LHS: $\bar\alpha\alpha / \eta^{ve}$;
           * RHS: $\bar(2\eta^{ve}\varepsilon^{eff} - k) / \eta^{ve}$.
           * Here $\alpha$ and $\bar\alpha$ correspond to the friction angle
           * and the dilation angle, respectively, $k$ is cohesion,
           * $\eta^{ve}$ is the non-yielding viscosity, and $\varepsilon^{eff}$
           * is the effective viscosity.
           */
          std::pair<double,double>
          compute_dilation_terms_for_stokes_system(const DruckerPragerParameters &drucker_prager_parameters,
                                                   const double non_yielding_viscosity,
                                                   const double effective_strain_rate) const;

        private:

          /**
           * The Drucker-Prager rheology is a simple plastic model
           * that yields at a stress of
           * (6.0 * cohesion * cos_phi + 6.0 * pressure * sin_phi) / (sqrt(3.0) * (3.0 + sin_phi))
           * in 3D or
           * cohesion * cos_phi + pressure * sin_phi in 2D.
           * Phi is an angle of internal friction, that is
           * input by the user in degrees, but stored as radians.
           */
          std::vector<double> angles_internal_friction;

          /**
           * The plastic potential of the Drucker-Prager model is given by
           * tau_II - 6.0 * pressure * sin_psi / (sqrt(3.0) * (3.0 + sin_psi)) in 3D or
           * tau_II - pressure * sin_psi in 2D, where tau_II is the second invariant
           * of the deviatoric stress.
           * Psi is an angle of dilation, that is input by the user in degrees,
           * but stored as radians. Note that the dilation angle must not be larger
           * than the internal friction angle.
           */
          std::vector<double> angles_dilation;

          /**
           * The cohesion is provided and stored in Pa.
           */
          std::vector<double> cohesions;

          /**
           *The prefactors for the yield stress.
           */
          std::vector<double> yield_stress_prefactors;

          /**
           * The yield stress is limited to a constant value, stored in Pa.
           */
          std::vector<double> max_yield_stresses;

          /**
           * Whether to add a plastic damper in the computation
           * of the drucker prager plastic viscosity.
           */
          bool use_plastic_damper;

          /**
           * Viscosity of a damper used to stabilize plasticity
           */
          double damper_viscosity;
      };
    }
  }
}
#endif
