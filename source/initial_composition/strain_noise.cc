/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/utilities.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/function_lib.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class Strain_noise : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Strain_noise ();

        /**
         * Initialization function.
         */
        void
        initialize ();

        /**
         * Return the initial composition as a function of position and number
         * of compositional field. The composition varies as a Gaussian distribution
         * around a user-defined set of line-segments.
         */
        virtual
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
    
        // Box with lithosphere 2d parameters     
        double x_extent_2d;
        double y_extent_2d_bl;


        /**
         * The number of the compositional field representing the strain.
         */
        unsigned int strain_composition_number;

        /**
         * Whether or not the domain is a cartesian box.
         */
        bool cartesian_domain = true;

        /**
         * The value of the seed for the random number generator
         */
        double seed;

        /**
         * The maximum amplitude of the Gaussian amplitude of the noise.
         */
        double A;

        /**
         * The standard deviation of the Gaussian amplitude of the noise.
         */
        double sigma;

        /**
         * The depth around and halfwidth with which the noise is smoothed out.
         */
        double strain_depth;
        double strain_halfwidth;

        /**
         * The list of line segments consisting of two 2d coordinates per segment (even in 2d).
         * The segments represent the Strain_noise axis.
         */
        std::vector<std::array<Point<2>,2 > > point_list;

        /**
         * A table with random noise for the
         * second invariant of the strain.
         */
        std::array<unsigned int,dim> grid_intervals;
        Functions::InterpolatedUniformGridData<dim> *interpolate_noise;
    };
  }
}


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    Strain_noise<dim>::Strain_noise ()
    {}

    template <int dim>
    void
    Strain_noise<dim>::initialize ()
    {
      AssertThrow(dynamic_cast<const MaterialModel::ViscoPlastic<dim> *>(&this->get_material_model()) != NULL,
                  ExcMessage("This initial condition only makes sense in combination with the visco_plastic material model."));

      // From shear_bands.cc
      Point<dim> extents_min, extents_max;
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        size_idx[d] = grid_intervals[d]+1;

      Table<dim,double> white_noise;
      white_noise.TableBase<dim,double>::reinit(size_idx);
      std::array<std::pair<double,double>,dim> grid_extents;


          const GeometryModel::Box<dim> *geometry_model
            = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());

            grid_extents[0].first = 0;
            grid_extents[1].first = 0;
            grid_extents[0].second = x_extent_2d;
            grid_extents[1].second = y_extent_2d_bl;


      // use a fixed number as seed for random generator
      // this is important if we run the code on more than 1 processor
      std::srand(seed);

      TableIndices<dim> idx;

      for (unsigned int i=0; i<white_noise.size()[0]; ++i)
        {
          idx[0] = i;
          for (unsigned int j=0; j<white_noise.size()[1]; ++j)
            {
              idx[1] = j;
                white_noise(idx) = ((std::rand() % 10000) / 10000.0);
            }
        }

      interpolate_noise = new Functions::InterpolatedUniformGridData<dim> (grid_extents,
                                                                           grid_intervals,
                                                                           white_noise);
    }

    template <int dim>
    double
    Strain_noise<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      // If n_comp does not represent the strain field,
      // return 0 right away.
      //std::cout<<"name composition :"<<this->introspection().name_for_compositional_index(strain_composition_number)<<std::endl;  
      if (n_comp != strain_composition_number)
        return 0.0;

      // Initiate distance with large value
      double distance_to_Strain_noise_axis = 1e23;
      double temp_distance = 0;

      // For spherical geometries we need to reorder the coordinates
      Point<dim> natural_coords = position;

      // Loop over all line segments
      for (unsigned int i_segments = 0; i_segments < point_list.size(); ++i_segments)
        {
          temp_distance = std::abs(natural_coords[0]-point_list[i_segments][0][0]);
          // Get the minimum distance
          distance_to_Strain_noise_axis = std::min(distance_to_Strain_noise_axis, temp_distance);
        }

      // Smoothing of noise with depth
      const double depth_smoothing = 0.5 * (1.0 - std::tanh((this->get_geometry_model().depth(position) - strain_depth) / strain_halfwidth));
      // Smoothing of noise with lateral distance to the Strain_noise axis
      const double noise_amplitude = A * std::exp((-std::pow(distance_to_Strain_noise_axis,2)/(2.0*std::pow(sigma,2)))) * depth_smoothing;
      // Add randomness
      return noise_amplitude * interpolate_noise->value(natural_coords);
    }

    template <int dim>
    void
    Strain_noise<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Strain_noise");
        {
          prm.enter_subsection("Box with lithosphere 2d");
          {
          prm.declare_entry("X extent in 2d", "2208000",
                            Patterns::Double(),
                            "Set a X extent in meters.");  
          prm.declare_entry("Y extent in 2d", "100000",
                            Patterns::Double(),
                            "Y extent when used with 2D");                                                
          }prm.leave_subsection();
            
          
          prm.declare_entry ("Random number generator seed", "0",
                             Patterns::Double (0),
                             "The value of the seed used for the random number generator. "
                             "Units: none.");
          prm.declare_entry ("Standard deviation of Gaussian noise amplitude distribution", "20000",
                             Patterns::Double (0),
                             "The standard deviation of the Gaussian distribution of the amplitude of the strain noise. "
                             "Note that this parameter is taken to be the same for all Strain_noise segments. "
                             "Units: $m$ or degrees.");
          prm.declare_entry ("Maximum amplitude of Gaussian noise amplitude distribution", "0.2",
                             Patterns::Double (0),
                             "The amplitude of the Gaussian distribution of the amplitude of the strain noise. "
                             "Note that this parameter is taken to be the same for all Strain_noise segments. "
                             "Units: none.");
          prm.declare_entry ("Depth around which Gaussian noise is smoothed out", "40000",
                             Patterns::Double (0),
                             "The depth around which smoothing out of the strain noise starts with a hyperbolic tangent. "
                             "Note that this parameter is taken to be the same for all Strain_noise segments. "
                             "Units: $m$.");
          prm.declare_entry ("Halfwidth with which Gaussian noise is smoothed out in depth", "40000",
                             Patterns::Double (0),
                             "The halfwidth with which smoothing out of the strain noise is done with a hyperbolic tangent. "
                             "Note that this parameter is taken to be the same for all Strain_noise segments. "
                             "Units: $m$.");
          prm.declare_entry ("Grid intervals for noise X or radius", "25",
                             Patterns::Integer (0),
                             "Grid intervals in X (cartesian domain) or radial (spherical) direction for the white noise "
                             "added to the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
          prm.declare_entry ("Grid intervals for noise Y or longitude", "25",
                             Patterns::Integer (0),
                             "Grid intervals in Y (cartesian domain) or longitude (spherical) direction for the white noise "
                             "added to the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
          prm.declare_entry ("Grid intervals for noise Z or latitude", "25",
                             Patterns::Integer (0),
                             "Grid intervals in Z (cartesian domain) or latitude (spherical) direction for the white noise "
                             "added to the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
          prm.declare_entry ("Strain_noise axis line segments",
                             "",
                             Patterns::Anything(),
                             "Set the line segments that represent the Strain_noise axis. In 3d each segment is made up of "
                             "two points that represent horizontal coordinates (x,y) or (lon,lat). "
                             "The exact format for the point list describing the segments is "
                             "\"x1,y1>x2,y2;x2,y2>x3,y3;x4,y4>x5,y5\". In 2d, a segment is made up by 1 horizontal "
                             "x or longitude coordinate: \"x1;x2;x3\". Note that the segments can be connected "
                             "or isolated. The units of the coordinates are "
                             "dependent on the geometry model. In the box model they are in meters, in the "
                             "chunks they are in degrees.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Strain_noise<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Check that there is a compositional field called strain and retrieve its index
      
      if(this->introspection().compositional_name_exists("plastic_strain"))
                  strain_composition_number = this->introspection().compositional_index_for_name("plastic_strain");
      else if(this->introspection().compositional_name_exists("viscous_strain"))
                  strain_composition_number = this->introspection().compositional_index_for_name("viscous_strain");
      else if(this->introspection().compositional_name_exists("total_strain"))
                  strain_composition_number = this->introspection().compositional_index_for_name("total_strain");
      else
          AssertThrow(false, ExcMessage("This plugin requires a compositional strain field (plastic_strain, viscous_strain, or total_strain). "));

      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Strain_noise");
        {
         
            sigma                = prm.get_double ("Standard deviation of Gaussian noise amplitude distribution");
            A                    = prm.get_double ("Maximum amplitude of Gaussian noise amplitude distribution");
            seed                 = prm.get_double ("Random number generator seed");
            strain_depth         = prm.get_double ("Depth around which Gaussian noise is smoothed out");
            strain_halfwidth     = prm.get_double ("Halfwidth with which Gaussian noise is smoothed out in depth");
            grid_intervals[0]    = prm.get_integer ("Grid intervals for noise X or radius");
            grid_intervals[1]    = prm.get_integer ("Grid intervals for noise Y or longitude");
            
           prm.enter_subsection("Box with lithosphere 2d");
                {
                    x_extent_2d = prm.get_double("X extent in 2d");
                    y_extent_2d_bl = prm.get_double("Y extent in 2d");            
                } 
                prm.leave_subsection();


        // Read in the string of Strain_noise segments
        const std::string temp_all_segments = prm.get("Strain_noise axis line segments");
        // Split the string into segment strings
        const std::vector<std::string> temp_segments = Utilities::split_string_list(temp_all_segments,';');
        // The number of segments, each consisting of a begin and an end point in 3d and one point in 3d
        const unsigned int n_temp_segments = temp_segments.size();
        point_list.resize(n_temp_segments);

        // Default is true, but in case we use a (ellipsoidal) chunk domain, set to false
//         if (dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model()) == NULL)
//           cartesian_domain = false;
         cartesian_domain = true;

        // Loop over the segments to extract the points
        for (unsigned int i_segment = 0; i_segment < n_temp_segments; i_segment++)
          {
            // In 3d a line segment consists of 2 points,
            // in 2d of only 1 (ridge axis orthogonal to x and y).
            // Also, a 3d point has 2 coordinates (x and y),
            // a 2d point only 1 (x).
            //point_list[i_segment].resize(dim-1);

                // Add the point to the list of points for this segment
                // As we're in 2d all segments correspond to 1 point consisting of 1 coordinate
                const double temp_point = Utilities::string_to_double(temp_segments[i_segment]);
                point_list[i_segment][0] = (Point<2>(temp_point, temp_point));
 
          }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Strain_noise,
                                              "strain_noise",
                                              "Specify the strain initial compositional field value based on the distance to a list of line segments "
                                              "and the user-defined Gaussian distribution around these segments, combined with random noise.")
  }
}