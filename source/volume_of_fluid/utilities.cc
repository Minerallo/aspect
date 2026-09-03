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

#include <aspect/volume_of_fluid/utilities.h>

#include <array>

namespace aspect
{
  namespace VolumeOfFluid
  {
    namespace Utilities
    {
      double compute_fluid_fraction (const Tensor<1, 2> normal,
                                     const double d)
      {
        const int dim = 2;

        //Get 1-Norm
        double norm1 = 0.0;
        double max = 0.0;
        for (unsigned int i = 0; i < dim; ++i)
          {
            const double normal_component = std::abs (normal[i]);
            norm1 += normal_component;
            max = (max < normal_component) ? normal_component : max;
          }

        //Obtain volume
        if (d <= -0.5*norm1)
          {
            return 0.0;
          }
        if (d >= 0.5*norm1)
          {
            return 1.0;
          }
        const double dtest = d / norm1; // Threshold value for changes in computation behavior
        // Normalized parameter indicating the "slope" (positive and finite)
        // Chosen due to resulting in simple formulas
        // Equal to the absolute value of the smaller vector entry in the 1-norm normalized normal vector
        const double mpos = 1.0 - max/norm1;
        if (dtest < mpos - 0.5)
          {
            return (dtest + 0.5) * (dtest + 0.5) / (2.0*mpos * (1.0 - mpos));
          }
        if (dtest > 0.5 - mpos)
          {
            return 1.0 - (dtest - 0.5) * (dtest - 0.5) / (2.0*mpos * (1.0 - mpos));
          }
        return 0.5 + dtest / (1.0 - mpos);
      }



      double compute_interface_location (const Tensor<1, 2> normal,
                                         const double vol)
      {
        const int dim = 2;

        //Get 1-Norm
        double norm1 = 0.0;
        double max = 0.0;
        for (unsigned int i = 0; i < dim; ++i)
          {
            double normal_component = std::abs (normal[i]);
            norm1 += normal_component;
            max = (max < normal_component) ? normal_component : max;
          }

        // Normalized parameter indicating the "slope" (positive and finite)
        // Chosen due to resulting in simple formulas
        // Equal to the absolute value of the smaller vector entry in the 1-norm normalized normal vector
        const double mpos = (norm1==0.0)?0.0:(1.0-max/norm1);
        norm1 = (norm1==0.0)? 1.0:norm1;

        // Obtain correct RHS term for interface position
        if (vol <= 0.0)
          {
            return -0.5 * norm1;
          }
        else if (vol >= 1.0)
          {
            return 0.5 * norm1;
          }
        else if (vol < 0.5 * mpos / (1 - mpos))
          {
            return norm1 * (-0.5 + std::sqrt (2.0*vol * mpos * (1 - mpos)));
          }
        else if (vol > 1.0 - 0.5 * mpos / (1 - mpos))
          {
            return norm1 * (0.5 - std::sqrt (2.0*(1.0 - vol) * mpos * (1 - mpos)));
          }
        else
          {
            return norm1 * (1 - mpos) * (vol - 0.5);
          }
      }



      double compute_fluid_fraction (const Tensor<1, 3> normal,
                                     const double d)
      {
        // Shift the centered unit cell [-1/2,1/2]^3 to [0,1]^3 and use
        // inclusion-exclusion for the volume of a clipped box. This is the
        // three-dimensional counterpart of the analytic 2d implementation
        // above, but remains well-defined when one or more normal components
        // vanish.
        std::array<double, 3> coefficients;
        unsigned int n_nonzero_components = 0;
        for (unsigned int direction = 0; direction < 3; ++direction)
          if (std::abs(normal[direction]) > std::numeric_limits<double>::epsilon())
            coefficients[n_nonzero_components++] = std::abs(normal[direction]);

        if (n_nonzero_components == 0)
          return (d < 0.0 ? 0.0 : (d > 0.0 ? 1.0 : 0.5));

        double normal_l1_norm = 0.0;
        for (unsigned int direction = 0; direction < n_nonzero_components; ++direction)
          normal_l1_norm += coefficients[direction];

        const double shifted_interface_location = d + 0.5 * normal_l1_norm;
        if (shifted_interface_location <= 0.0)
          return 0.0;
        if (shifted_interface_location >= normal_l1_norm)
          return 1.0;

        if (n_nonzero_components == 1)
          return shifted_interface_location / coefficients[0];

        if (n_nonzero_components == 2)
          {
            Tensor<1, 2> reduced_normal;
            reduced_normal[0] = coefficients[0];
            reduced_normal[1] = coefficients[1];
            return compute_fluid_fraction(reduced_normal, d);
          }

        const auto positive_cube = [] (const double value)
        {
          const double positive_part = std::max(0.0, value);
          return positive_part * positive_part * positive_part;
        };

        double volume = positive_cube(shifted_interface_location);
        for (unsigned int i = 0; i < 3; ++i)
          volume -= positive_cube(shifted_interface_location - coefficients[i]);
        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = i + 1; j < 3; ++j)
            volume += positive_cube(shifted_interface_location
                                    - coefficients[i] - coefficients[j]);
        volume -= positive_cube(shifted_interface_location - normal_l1_norm);

        volume /= 6.0 * coefficients[0] * coefficients[1] * coefficients[2];
        return std::max(0.0, std::min(1.0, volume));
      }



      double compute_interface_location (const Tensor<1, 3, double> normal,
                                         const double vol)
      {
        double normal_l1_norm = 0.0;
        for (unsigned int direction = 0; direction < 3; ++direction)
          normal_l1_norm += std::abs(normal[direction]);

        if (normal_l1_norm == 0.0)
          return 0.0;
        if (vol <= 0.0)
          return -0.5 * normal_l1_norm;
        if (vol >= 1.0)
          return 0.5 * normal_l1_norm;

        // The clipped volume is continuous and monotone in d. A bounded
        // bisection is inexpensive compared with interface reconstruction and
        // avoids the singular special cases of the closed-form inverse.
        double lower_bound = -0.5 * normal_l1_norm;
        double upper_bound = 0.5 * normal_l1_norm;
        for (unsigned int iteration = 0; iteration < 64; ++iteration)
          {
            const double midpoint = 0.5 * (lower_bound + upper_bound);
            if (compute_fluid_fraction(normal, midpoint) < vol)
              lower_bound = midpoint;
            else
              upper_bound = midpoint;
          }
        return 0.5 * (lower_bound + upper_bound);
      }



      void xFEM_Heaviside(const unsigned int degree,
                          const Tensor<1, 2> normal,
                          const double d,
                          const std::vector<Point<2>> &points,
                          std::vector<double> &values)
      {
        const int basis_count=4;
        std::vector<double> coeffs(basis_count);

        const double n_xp = std::fabs(normal[0]), n_yp = std::fabs(normal[1]);
        const double sign_n_x = (((normal[0]) > 0) - ((normal[0]) < 0)),
                     sign_n_y = (((normal[1]) > 0) - ((normal[1]) < 0));

        const double norm1 = n_xp + n_yp;
        const double triangle_break = 0.5*std::fabs(n_xp-n_yp);

        const unsigned int max_degree = 1;

        AssertThrow(degree<=max_degree,
                    ExcMessage("Cannot generate xFEM polynomials. Only implemented for degrees<2."));

        // The formulas below calculate the correct coefficients for a given
        // basis in order to form a polynomial $f$ which will satisfy $\int
        // fpdx=\int pH(d-n\cdot x)dx$ for all polynomials $p$ less than or
        // equal to the given degree
        //
        // The functions for the correct values were calculated and exported using sympy
        if (d<-0.5*norm1)
          {
            for (unsigned int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d>0.5*norm1)
          {
            // Full cell
            coeffs[0] = 1.0;
            for (unsigned int i =1; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (norm1< 1e-7)
          {
            coeffs[0] = 0.5;
            for (unsigned int i =1; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d<=-triangle_break)
          {
            //Triangle
            const double d_n = d + 0.5*norm1;
            coeffs[0]=0.5*d_n*d_n/(n_xp*n_yp); // 1
            coeffs[1]=d_n*d_n*(d_n - 1.5*n_yp)/(n_xp*n_yp*n_yp)*sign_n_y; // 2*y - 1
            coeffs[2]=d_n*d_n*(d_n - 1.5*n_xp)/(n_xp*n_xp*n_yp)*sign_n_x; // 2*x - 1
            coeffs[3]=1.5*d_n*d_n*(d_n*d_n - 2.0*d_n*n_xp - 2.0*d_n*n_yp + 3.0*n_xp*n_yp)/
                      (n_xp*n_xp*n_yp*n_yp)*sign_n_x*sign_n_y; // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_xp<n_yp)
          {
            //Trapezoid X
            coeffs[0]=(d + 0.5*n_yp)/n_yp; // 1
            coeffs[1]=0.25*(12.0*d*d + n_xp*n_xp - 3.0*n_yp*n_yp)/(n_yp*n_yp)*sign_n_y; // 2*y - 1
            coeffs[2]=-0.5*n_xp/n_yp*sign_n_x; // 2*x - 1
            coeffs[3]=-3.0*d*n_xp/(n_yp*n_yp)*sign_n_x*sign_n_y; // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_yp<n_xp)
          {
            //Trapezoid Y
            coeffs[0]=(d + 0.5*n_xp)/n_xp; // 1
            coeffs[1]=-0.5*n_yp/n_xp*sign_n_y; // 2*y - 1
            coeffs[2]=0.25*(12.0*dealii::Utilities::fixed_power<2>(d) - 3.0*n_xp*n_xp + n_yp*n_yp)/(n_xp*n_xp)*sign_n_x; // 2*x - 1
            coeffs[3]=-3.0*d*n_yp/(n_xp*n_xp)*sign_n_x*sign_n_y; // (2*x - 1)*(2*y - 1)
          }
        else
          {
            //ITriangle
            const double d_nn = 0.5*norm1-d;
            coeffs[0]=1.0L-0.5L*d_nn*d_nn/(n_xp*n_yp); // 1
            coeffs[1]=0.5L*(d_nn*d_nn)*sign_n_y*(2*d_nn - 3*n_yp)/(n_xp*(n_yp*n_yp)); // 2*y - 1
            coeffs[2]=0.5L*(d_nn*d_nn)*sign_n_x*(2*d_nn - 3*n_xp)/((n_xp*n_xp)*n_yp); // 2*x - 1
            coeffs[3]=1.5L*(d_nn*d_nn)*sign_n_x*sign_n_y*(-(d_nn*d_nn) + 2*d_nn*n_xp + 2*d_nn*n_yp - 3*n_xp*n_yp)/((n_xp*n_xp)*(n_yp*n_yp)); // (2*x - 1)*(2*y - 1)
          }

        // Calculate the correct values at the provided quadrature points by
        // multiplying coefficients by the basis polynomials.
        for (unsigned int i = 0; i<points.size(); ++i)
          {
            const Point<2> point = points[i];
            const double x = point[0], y = point[1];
            values[i] = coeffs[0];
            if (degree>=1)
              {
                values[i] += coeffs[1]*(2.0*y-1.0) +
                             coeffs[2]*(2.0*x-1.0) +
                             coeffs[3]*(2.0*x - 1)*(2.0*y - 1.0);
              }
          }
      }



      void xFEM_Heaviside(const unsigned int /*degree*/,
                          const Tensor<1, 3> /*normal*/,
                          const double /*d*/,
                          const std::vector<Point<3>> &/*points*/,
                          std::vector<double> &/*values*/)
      {
        AssertThrow(false, ExcNotImplemented());
      }



      void xFEM_Heaviside_derivative_d(const unsigned int degree,
                                       const Tensor<1, 2> normal,
                                       const double d,
                                       const std::vector<Point<2>> &points,
                                       std::vector<double> &values)
      {
        const int basis_count=4;
        std::vector<double> coeffs(basis_count);

        const double n_xp = std::fabs(normal[0]), n_yp = std::fabs(normal[1]);
        const double sign_n_x = (((normal[0]) > 0) - ((normal[0]) < 0)),
                     sign_n_y = (((normal[1]) > 0) - ((normal[1]) < 0));

        const double norm1 = n_xp + n_yp;
        const double triangle_break = 0.5L*std::fabs(n_xp-n_yp);

        const int max_degree = 1;

        AssertThrow(degree<=max_degree,
                    ExcMessage("Cannot generate xFEM polynomials are only functional for degrees<2."));


        // The formulas below calculate the correct coefficients for a given
        // basis in order to form a polynomial $f$ which will satisfy $\int
        // fpdx=\int pH(d-n\cdot x)dx$ for all polynomials $p$ less than or
        // equal to the given degree
        //
        // The functions for the correct values were calculated and exported using sympy
        if (d<-0.5*norm1)
          {
            for (int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d>0.5*norm1)
          {
            // Full cell
            for (int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (norm1<1e-7)
          {
            for (int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d<=-triangle_break)
          {
            //D Triangle
            const double d_n = d + 0.5*norm1;
            coeffs[0]=d_n/(n_xp*n_yp); // 1
            coeffs[1]=3*d_n*sign_n_y*(d_n - n_yp)/(n_xp*(n_yp*n_yp)); // 2*y - 1
            coeffs[2]=3*d_n*sign_n_x*(d_n - n_xp)/((n_xp*n_xp)*n_yp); // 2*x - 1
            coeffs[3]=3*d_n*sign_n_x*sign_n_y*(2*(d_n*d_n) - 3*d_n*n_xp - 3*d_n*n_yp + 3*n_xp*n_yp)/((n_xp*n_xp)*(n_yp*n_yp)); // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_xp<n_yp)
          {
            //D Trapezoid X
            coeffs[0]=1.0/n_yp; // 1
            coeffs[1]=6*d*sign_n_y/(n_yp*n_yp); // 2*y - 1
            coeffs[2]=0; // 2*x - 1
            coeffs[3]=-3*n_xp*sign_n_x*sign_n_y/(n_yp*n_yp); // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_yp<n_xp)
          {
            //D Trapezoid Y
            coeffs[0]=1.0/n_xp; // 1
            coeffs[1]=0; // 2*y - 1
            coeffs[2]=6*d*sign_n_x/(n_xp*n_xp); // 2*x - 1
            coeffs[3]=-3*n_yp*sign_n_x*sign_n_y/(n_xp*n_xp); // (2*x - 1)*(2*y - 1)
          }
        else
          {
            //D ITriangle
            const double d_nn = 0.5*norm1-d;
            coeffs[0]=d_nn/(n_xp*n_yp); // 1
            coeffs[1]=3*d_nn*sign_n_y*(-d_nn + n_yp)/(n_xp*(n_yp*n_yp)); // 2*y - 1
            coeffs[2]=3*d_nn*sign_n_x*(-d_nn + n_xp)/((n_xp*n_xp)*n_yp); // 2*x - 1
            coeffs[3]=3*d_nn*sign_n_x*sign_n_y*(2*(d_nn*d_nn) - 3*d_nn*n_xp - 3*d_nn*n_yp + 3*n_xp*n_yp)/((n_xp*n_xp)*(n_yp*n_yp)); // (2*x - 1)*(2*y - 1)
          }


        // Calculate the correct values at the provided quadrature points by
        // multiplying coefficients by the basis polynomials.
        for (unsigned int i = 0; i<points.size(); ++i)
          {
            const Point<2> point = points[i];
            const double x = point[0], y = point[1];
            values[i] = coeffs[0];
            if (degree>=1)
              {
                values[i] += coeffs[1]*(2*y-1.0) +
                             coeffs[2]*(2*x-1.0) +
                             coeffs[3]*(2*x - 1)*(2*y - 1);
              }
          }
      }



      void xFEM_Heaviside_derivative_d(const unsigned int /*degree*/,
                                       const Tensor<1, 3> /*normal*/,
                                       const double /*d*/,
                                       const std::vector<Point<3>> &/*points*/,
                                       std::vector<double> &/*values*/)
      {
        AssertThrow(false, ExcNotImplemented());
      }



      template <int dim>
      double compute_interface_location_newton(const unsigned int degree,
                                               const Tensor<1, dim> normal,
                                               const double volume_fraction,
                                               const double vol,
                                               const double epsilon,
                                               const std::vector<Point<dim>> &points,
                                               const std::vector<double> &weights)
      {
        double norm1=0.0;
        for (int i=0; i<dim; ++i)
          norm1+=std::fabs(normal[i]);
        double d_l=-0.5L*norm1, d_h=0.5L*norm1;
        double f_l=0.0, f_h=1.0;
        double d_guess= d_l + (volume_fraction-f_l)*(d_h-d_l)/(f_h-f_l);

        std::vector<double> f_values(points.size());
        std::vector<double> df_values(points.size());

        for (int iter=0; iter<40; ++iter)
          {
            xFEM_Heaviside(degree, normal, d_guess, points, f_values);
            xFEM_Heaviside_derivative_d(degree, normal, d_guess, points, df_values);

            double f_guess=0.0;
            double df_guess=0.0;
            for (unsigned int i=0; i<points.size(); ++i)
              {
                const double factor = weights[i]/vol;
                f_guess  += f_values[i]*factor;
                df_guess += df_values[i]*factor;
              }

            // Break if within tolerance
            if (std::fabs(f_guess-volume_fraction)<epsilon)
              {
                break;
              }

            if (volume_fraction<f_guess)
              {
                d_h = d_guess;
                f_h = f_guess;
              }
            else
              {
                d_l = d_guess;
                f_l = f_guess;
              }

            if (std::fabs(df_guess)<epsilon)
              {
                d_guess = (volume_fraction-f_l)/(f_h-f_l)*(d_h-d_l);
              }
            else
              {
                d_guess += (volume_fraction-f_guess)/(df_guess);

                if (d_guess < d_l || d_guess > d_h)
                  {
                    d_guess = d_l + (volume_fraction-f_l)*(d_h-d_l)/(f_h-f_l);
                  }
              }
          }

        return d_guess;
      }

      template <int dim>
      double compute_fluid_volume(const unsigned int degree,
                                  const Tensor<1, dim> normal,
                                  const double d,
                                  const std::vector<Point<dim>> &points,
                                  const std::vector<double> &weights)
      {
        std::vector<double> f_values(points.size());

        xFEM_Heaviside(degree, normal, d, points, f_values);

        double fluid_volume=0.0;
        for (unsigned int i=0; i<points.size(); ++i)
          fluid_volume += f_values[i]*weights[i];

        return fluid_volume;
      }

      template <int dim>
      double calculate_volume_flux(const unsigned int dir,
                                   const double time_direction_derivative,
                                   const Tensor<1, dim> normal,
                                   const double d_face)
      {
        Tensor<1, dim> i_normal;
        double i_d;

        // Get d value at center of "aperture" (cell face cross timestep)
        i_d = d_face+0.5*time_direction_derivative;
        // Get normal vector on "aperture" by replacing the appropriate component with the time_direction_derivative
        i_normal = normal;
        i_normal[dir] = time_direction_derivative;

        return compute_fluid_fraction (i_normal, i_d);
      }
    }
  }
}

namespace aspect
{
  namespace VolumeOfFluid
  {
    namespace Utilities
    {
#define INSTANTIATE(dim) \
  template double calculate_volume_flux<dim>(const unsigned int dir, \
                                             const double time_direction_derivative, \
                                             const Tensor<1, dim> normal, \
                                             const double d); \
  template double compute_interface_location_newton<dim>(const unsigned int degree, \
                                                         const Tensor<1, dim> normal, \
                                                         const double volume_fraction, \
                                                         const double vol, \
                                                         const double epsilon, \
                                                         const std::vector<Point<dim>> &points, \
                                                         const std::vector<double> &weights); \
  template double compute_fluid_volume<dim>(const unsigned int degree, \
                                            const Tensor<1, dim> normal,\
                                            const double d,\
                                            const std::vector<Point<dim>> &points,\
                                            const std::vector<double> &weights);

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
