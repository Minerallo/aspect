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
    class Trenchposition : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
        
//         const std::set<types::boundary_id> &boundary_id;
    };

    template <int dim>
    std::pair<std::string,std::string>
    Trenchposition<dim>::execute (TableHandler &statistics)
    {
      QTrapez<dim-1> face_corners;
      FEFaceValues<dim> fe_face_values (this->get_mapping(), this->get_fe(), face_corners, update_quadrature_points);


      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output_stats;
      // std::ostringstream output_file;


      // Choose stupidly large values for initialization
      double local_min_z = std::numeric_limits<double>::max();                               
      double local_max_z = std::numeric_limits<double>::min();
      double index_min_z =std::numeric_limits<double>::max(); 
        

        std::vector<double>compositional_values(face_corners.size()); 
        
        const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
//         const bool top_boundary = boundary_ids.find(relevant_boundary) == boundary_ids.begin();

      // loop over all cells and find cell with 100 % of the defined composition, then save the elevation to stored_value
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
        {
            

          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
              if (cell->face(face_no)->at_boundary())
              {
                if (cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;
                
                fe_face_values.reinit(cell, face_no);
                
                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                      const Point<dim> vertex = fe_face_values.quadrature_point(corner);
                      

                    if (vertex(1) < local_min_z)
                    {
                      local_min_z = vertex(1);
                      index_min_z  = vertex(0);
                    }              
                }           
          }        
        }
      }
      }

//       for (unsigned int p=0; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
//            {
//             min_z = std::abs(Utilities::MPI::min(local_min_z
//             if (min_z < local_min_z)
//                     {
//                         index_min_z = p; 
//                     }
//         const double min_z = std::abs(Utilities::MPI::min(local_min_z, this->get_mpi_communicator()));
//            }
        
              const double min_z = std::abs(Utilities::MPI::min(local_min_z, this->get_mpi_communicator()));

//       std::cout<<"min_x :"<<min_x<<" max_x :"<<max_x<<std::endl;    
      
      //Write results to statistics file
      statistics.add_value ("Trench position (m)",
                            min_z);
      const char *columns[] = { "Trench position (m)"};
      for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
        {
          statistics.set_precision (columns[i], 8);
          statistics.set_scientific (columns[i], true);
        }

      output_stats.precision(5);
      output_stats << min_z << " m";

      // Just return stats if text output is not required at all or not needed at this time
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return std::pair<std::string, std::string> (" Trench position (m)",
                                                    output_stats.str());


      return std::pair<std::string, std::string> ("Tench positiom (m)",
                                                  output_stats.str());

    }      
  


    template <int dim>
    void Trenchposition<dim>::declare_parameters(ParameterHandler &prm)
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
        prm.enter_subsection("Trench position");
        {        
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
          prm.declare_entry ("Output trench position dependant ending file", "false",
                             Patterns::List(Patterns::Bool()),
                             "Output a Termin aspect file when a depth that user set has been reached by the choosen composition"
                             ", for example, in a subduction model, when the slab reach a defined depth"
                             ", a fill is outputed and using the request user termination criterion, the simulation can stop");   
          prm.declare_entry ("Ending trench position criterion", "2",
                             Patterns::Double (0.),
                             "The depth choosen to ouput a random file that cam be used as an termine aspect file"); 
        
        }
          prm.leave_subsection();        
      }
      prm.leave_subsection();    
    }                           

    template <int dim>
    void
    Trenchposition<dim>::parse_parameters (ParameterHandler &prm)
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
        prm.enter_subsection("Trench position");
        {
          write_to_file      = prm.get_bool ("Output to file");
          output_interval    = prm.get_double ("Time between text output");                            
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
    ASPECT_REGISTER_POSTPROCESSOR(Trenchposition,
                                  "Trench position",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
