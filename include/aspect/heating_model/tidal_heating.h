/*
  Copyright (C) 2025 by the authors of the ASPECT code.

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


#ifndef _aspect_heating_model_tidal_heating_h
#define _aspect_heating_model_tidal_heating_h

#include <aspect/heating_model/interface.h>

namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements tidal heating.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class TidalHeating : public Interface<dim>
    {
      public:
        /**
         * Return the tidal heating terms.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        /**
         * Specify which material model outputs the heating model requires
         * for computing the heating terms.
         */
        MaterialModel::MaterialProperties::Property
        get_required_properties () const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Parameters used for tidal heating (H), which is defined using the following
         * Equation is from Tobie et al. (2003) (https://doi.org/10.1029/2003JE002099)
         * H = 2*(viscosity)*(time-averaged tidal strain rate)^2/(1+((viscosity)*(tidal frequency)/(elastic shear modulus))^2))
         * viscosity (Pa s) = viscosity calculated by the selected material in ASPECT
         * time-averaged strain rate = constant_tidal_strain_rate
         * tidal frequency = tidal_frequency
         * elastic shear modulus = elastic_shear_modulus
         */
        double tidal_frequency;
        double elastic_shear_modulus;
        double constant_tidal_strain_rate;
    };
  }
}


#endif
