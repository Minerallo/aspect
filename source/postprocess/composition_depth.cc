/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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

#include <aspect/postprocess/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/geometry_model/box.h>
#include <aspect/simulator.h>
#include <deal.II/base/quadrature_lib.h>


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <cmath>
#include <limits>

#include <aspect/geometry_model/two_merged_boxes.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the compositional
     * fields, if any.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class CompositionDepth : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some temperature statistics.
         */
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the free surface handling.
         */
        void parse_parameters (ParameterHandler &prm);        

      private:
        /**
         * The compositional field number, min ref level and max ref level
         * for the crust, mantle part of the slab lithosphere and the
         * overriding plate
         */
        unsigned int composition_number; 
        /**
         * Whether or not to produce text files with topography values
         */
        bool write_to_file;
        bool write_to_end_file;
        double ending_depth;

        /**
         * Interval between the generation of text output. This parameter
         * is read from the input file and consequently is not part of the
         * state that needs to be saved and restored.
         */
        double output_interval;

        /**
         * A time (in seconds) at which the last text output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        double last_output_time;
        double analytical_solution_example; 
       /**
         * Extent of the whole model domain in x-, y-, and z-direction (in 3d).
         */
        Point<dim> extents;               
    };

    template <int dim>
    std::pair<std::string,std::string>
    CompositionDepth<dim>::execute (TableHandler &statistics)
    {
      QTrapez<dim-1> face_corners;

      // const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).get_unit_support_points());
        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                         this->get_fe(),
                                         face_corners,
                                         update_values |
                                             update_quadrature_points);
      // FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, quadrature, update_quadrature_points | update_values | update_gradients);

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output_stats;
      // std::ostringstream output_file;


      // Choose stupidly large values for initialization
      double local_min_depth = std::numeric_limits<double>::max();                               

    //   MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
    //   MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
        
        // std::vector<std::vector<double>compositional_values
        // (this->n_compositional_fields(),
        // std::vector<double> (quadrature.size()));

        std::vector<double>compositional_values(face_corners.size()); 
        // double compositional_values =0.0;
        // std::vector<std::vector<double>> temporary_variables(0, std::vector<double>());


      // loop over all cells and find cell with 100 % of the defined composition, then save the elevation to stored_value
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
        {
            bool composition_tracked_present = false;   

          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.compositional_fields[composition_number]].get_function_values (this->get_solution(),
                          compositional_values);
                
                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                      const Point<dim> vertex = fe_face_values.quadrature_point(corner);
                      const double elevation = this->get_geometry_model().height_above_reference_surface(vertex);  

                    // std::cout<<corner<<"  "<<compositional_values[corner]<<std::endl;

                  // for (unsigned int p=0; p<face_corners.size(); ++p)
                  //   {     
                      if (compositional_values[corner] >= 0.95)
                        {
                          composition_tracked_present = true;
                        }                        
                    //  } 

               if(composition_tracked_present)
               {
                    if (elevation < local_min_depth)
                      local_min_depth = elevation;                   
               } 
                }           
          }        
        }
      }

      //Calculate max depth across all processes
      const double min_depth = extents[1] - Utilities::MPI::min(local_min_depth, this->get_mpi_communicator());        
      // std::cout<<"Maximum depth value for composition :"<<local_min_height<<std::endl;

      // statistics.add_value ("Maximum depth value for composition " + this->introspection().name_for_compositional_index(composition_number),
      //                           local_min_height);

      if(write_to_end_file)
      {
        if(min_depth>=ending_depth)
          {
            std::string filename = this->get_output_directory() +"terminate-aspect";            
            std::ofstream file (filename.c_str());
            file << "This is a generated random file generated when the tracked composition reached the defined depth of "<<ending_depth<<" m, it can be later used to end the run using the user request termination criteria " << std::endl;
          }        
      }

      //Write results to statistics file
      statistics.add_value ("Depth (m)",
                            min_depth);
      const char *columns[] = { "Depth (m)"};
      for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
        {
          statistics.set_precision (columns[i], 8);
          statistics.set_scientific (columns[i], true);
        }

      output_stats.precision(4);
      output_stats << min_depth << " m";

      // Just return stats if text output is not required at all or not needed at this time
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return std::pair<std::string, std::string> (" Depth composition seclected (m):",
                                                    output_stats.str());


      return std::pair<std::string, std::string> ("Depth (m):",
                                                  output_stats.str());

      // statistics.add_value ("Depth (m)",
      //                       min_depth);
      // const char columns[] = { "depth (m)"};      

      // // for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
      // //   {
      // //     statistics.set_precision (columns[i], 8);
      // //     statistics.set_scientific (columns[i], true);
      // //   }

      // output_stats.precision(4);
      // output_stats << min_depth<< " m ";


      // // Write the solution to file

      // // if this is the first time we get here, set the last output time
      // // to the current time - output_interval. this makes sure we
      // // always produce data during the first time step
      // if (std::isnan(last_output_time))
      //   {
      //     last_output_time = this->get_time() - output_interval;
      //   }

      // // Just return stats if text output is not required at all or not needed at this time
      // if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
      //                        && (this->get_timestep_number() != 0)))
      //   return std::pair<std::string, std::string> ("Depth",
      //                                               output_stats.str());

      // std::string filename = this->get_output_directory() +
      //                        "Depth" +
      //                        Utilities::int_to_string(this->get_timestep_number(), 5);
      // if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
      //   filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      // const unsigned int max_data_length = Utilities::MPI::max (output_file.str().size()+1,
      //                                                           this->get_mpi_communicator());

      // const unsigned int mpi_tag = 777;

      // // on processor 0, collect all of the data the individual processors send
      // // and concatenate them into one file
      // if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
      //   {
      //     std::ofstream file (filename.c_str());

      //     file << "# Max depth for composition number : " <<composition_number<< std::endl;

      //     // first write out the data we have created locally
      //     file << output_file.str();

      //     std::string tmp;
      //     tmp.resize (max_data_length, '\0');

      //     // then loop through all of the other processors and collect
      //     // data, then write it to the file
      //     for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
      //       {
      //         MPI_Status status;
      //         // get the data. note that MPI says that an MPI_Recv may receive
      //         // less data than the length specified here. since we have already
      //         // determined the maximal message length, we use this feature here
      //         // rather than trying to find out the exact message length with
      //         // a call to MPI_Probe.
      //         const int ierr = MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
      //                                    this->get_mpi_communicator(), &status);
      //         AssertThrowMPI(ierr);

      //         // output the string. note that 'tmp' has length max_data_length,
      //         // but we only wrote a certain piece of it in the MPI_Recv, ended
      //         // by a \0 character. write only this part by outputting it as a
      //         // C string object, rather than as a std::string
      //         file << tmp.c_str();
      //       }
      //   }
      // else
      //   // on other processors, send the data to processor zero. include the \0
      //   // character at the end of the string
      //   {
      //     output_file << "\0";
      //     const int ierr = MPI_Send (&output_file.str()[0], output_file.str().size()+1, MPI_CHAR, 0, mpi_tag,
      //                                this->get_mpi_communicator());
      //     AssertThrowMPI(ierr);
      //   }

      // // if output_interval is positive, then update the last supposed output
      // // time
      // if (output_interval > 0)
      //   {
      //     // We need to find the last time output was supposed to be written.
      //     // this is the last_output_time plus the largest positive multiple
      //     // of output_intervals that passed since then. We need to handle the
      //     // edge case where last_output_time+output_interval==current_time,
      //     // we did an output and std::floor sadly rounds to zero. This is done
      //     // by forcing std::floor to round 1.0-eps to 1.0.
      //     const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
      //     last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
      //   }

      // return std::pair<std::string, std::string> ("Depth:",
      //                                             output_stats.str());
    }      
  


    template <int dim>
    void CompositionDepth<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          // Total box extents
          prm.declare_entry ("X extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in x-direction. Units: $\\si{m}$.");
          prm.declare_entry ("Y extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in y-direction. Units: $\\si{m}$.");            
                  }
        prm.leave_subsection();
      }
      prm.leave_subsection(); 

     prm.enter_subsection("Postprocess");
      {      
        prm.enter_subsection("Composition Depth");
        {        
            prm.declare_entry("Composition tracked number", "1",
                              Patterns::Integer(),
                              "composition number that is tracked. ");
          prm.declare_entry ("Time between text output", "0.",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "text output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise."); 
          prm.declare_entry ("Output to file", "false",
                             Patterns::List(Patterns::Bool()),
                             "Whether or not to write topography to a text file named named "
                             "'topography.NNNNN' in the output directory"); 
          prm.declare_entry ("Output depth dependant ending file", "true",
                             Patterns::List(Patterns::Bool()),
                             "Output a Termin aspect file when a depth that user set has been reached by the choosen composition"
                             ", for example, in a subduction model, when the slab reach a defined depth"
                             ", a fill is outputed and using the request user termination criterion, the simulation can stop");   
          prm.declare_entry ("Ending depth criterion", "300000",
                             Patterns::Double (0.),
                             "The depth choosen to ouput a random file that cam be used as an termine aspect file");                                          
        }
          prm.leave_subsection();        
      }
      prm.leave_subsection();    
    }                           

    template <int dim>
    void
    CompositionDepth<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          // Total box extents
          extents[0]           = prm.get_double ("X extent");
          extents[1]           = prm.get_double ("Y extent");
                  }
        prm.leave_subsection();
      }
      prm.leave_subsection(); 

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition Depth");
        {
          composition_number = prm.get_integer("Composition tracked number");
          write_to_file      = prm.get_bool ("Output to file");
          output_interval    = prm.get_double ("Time between text output");
          write_to_end_file  = prm.get_bool ("Output depth dependant ending file"); 
          ending_depth       = prm.get_double ("Ending depth criterion");                                      
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
    ASPECT_REGISTER_POSTPROCESSOR(CompositionDepth,
                                  "composition depth",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
