 
/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/postprocess/topolayer.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/geometry_model/two_merged_boxes.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <cmath>
#include <limits>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    Topolayer<dim>::execute (TableHandler &statistics)
    {
      //Disallow use of the plugin with sphere geometry model
      AssertThrow(!Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()),
                  ExcMessage("Topolayer postprocessor is not yet implemented "
                             "for the sphere geometry model. "
                             "Consider using a box, spherical shell, or chunk.") );

//       const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

      // Get a quadrature rule that exists only on the corners
#if DEAL_II_VERSION_GTE(9,3,0)
      const QTrapezoid<dim-1> face_corners;
#else
      const QTrapez<dim-1> face_corners;
#endif
//       FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, update_quadrature_points);
        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                         this->get_fe(),
                                         face_corners,
                                         update_values |
                                             update_quadrature_points);
      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output_stats;
      std::ostringstream output_file;

      // Choose stupidly large values for initialization
      double local_max_height = std::numeric_limits<double>::lowest();
      double local_min_height = std::numeric_limits<double>::max();
      
      std::vector<double>compositional_values(face_corners.size()); 

      // loop over all of the surface cells and save the elevation to stored_value
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
        {
          bool composition_tracked_present = false;
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
//             if (cell->face(face_no)->at_boundary())
//               {
//                 if ( cell->face(face_no)->boundary_id() != relevant_boundary)
//                   continue;

                fe_face_values.reinit( cell, face_no);
                fe_face_values[this->introspection().extractors.compositional_fields[composition_tracked]].get_function_values (this->get_solution(),
                          compositional_values);                

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);
                    
                       if (compositional_values[corner] >= 0.2)
                        {
                          composition_tracked_present = true;
                        }                     
                        
                       if(composition_tracked_present)
                        {                                                    
//                             output_file << vertex << ' '<< vertex(1) << std::endl;
                            if ( vertex(1) > local_max_height)
                            local_max_height = vertex(1);
                            if (write_to_file)
                                output_file << vertex << ' ' <<std::endl;                              
//                             if ( vertex(1) < local_min_height)
//                             local_min_height = vertex(1);
                        }
                  }
//               }
        }
        }
      }

      //Calculate min/max topolayer across all processes
      const double max_topolayer = Utilities::MPI::max(local_max_height, this->get_mpi_communicator());
//       const double min_topolayer = Utilities::MPI::min(local_min_height, this->get_mpi_communicator());

//       //Write results to statistics file
//       statistics.add_value ("Minimum topolayer (m)",
//                             min_topolayer);
//       statistics.add_value ("Maximum topolayer (m)",
//                             max_topolayer);
//       const char *columns[] = { "Minimum topolayer (m)",
//                                 "Maximum topolayer (m)"
//                               };
//       for (auto &column : columns)
//         {
//           statistics.set_precision (column, 8);
//           statistics.set_scientific (column, true);
//         }
// 
//       output_stats.precision(4);
//       output_stats << min_topolayer << " m, "
//                    << max_topolayer << " m";

      // Write the solution to file

      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
        }

      // Just return stats if text output is not required at all or not needed at this time
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return std::pair<std::string, std::string> ("Topolayer min/max:",
                                                    output_stats.str());

      std::string filename = this->get_output_directory() +
                             "topolayer." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      const std::vector<std::string> data = Utilities::MPI::gather(this->get_mpi_communicator(), output_file.str());

      // On processor 0, collect all of the data the individual processors sent
      // and concatenate them into one file:
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << ((dim==2)? "x y" : "x y z")
               << " topolayer" << std::endl;

          for (const auto &str : data)
            file << str;
        }

      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
        }

      return std::pair<std::string, std::string> ("Topolayer min/max:",
                                                  output_stats.str());
    }

    template <int dim>
    void
    Topolayer<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Topolayer");
        {
          prm.declare_entry ("Output to file", "false",
                             Patterns::List(Patterns::Bool()),
                             "Whether or not to write topolayer to a text file named named "
                             "'topolayer.NNNNN' in the output directory");

          prm.declare_entry ("Time between text output", "0.",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "text output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
            prm.declare_entry("Composition tracked", "1",
                              Patterns::Integer(),
                              "composition number that is tracked. ");          
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Topolayer<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Topolayer");
        {
          write_to_file        = prm.get_bool ("Output to file");
          output_interval = prm.get_double ("Time between text output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
          composition_tracked= prm.get_integer("Composition tracked");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Topolayer,
                                  "topolayer",
                                  "A postprocessor intended for use with a deforming top surface.  After every step "
                                  "it loops over all the vertices on the top surface and determines the "
                                  "maximum and minimum topolayer relative to a reference datum (initial "
                                  "box height for a box geometry model or initial radius for a "
                                  "sphere/spherical shell geometry model). "
                                  "If 'Topolayer.Output to file' is set to true, also outputs topolayer "
                                  "into text files named `topolayer.NNNNN' in the output directory, "
                                  "where NNNNN is the number of the time step.\n"
                                  "The file format then consists of lines with Euclidean coordinates "
                                  "followed by the corresponding topolayer value."
                                  "Topolayer is printed/written in meters.")
  }
}