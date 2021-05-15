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
    class ExtentComposition : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
        unsigned int composition_extent; 
        /**
         * Whether or not to produce text files with topography values
         */
        bool write_to_file;
        bool write_to_end_file;
        double ending_depth;
        double composition_depth;

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
    ExtentComposition<dim>::execute (TableHandler &statistics)
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
      double local_min_x = std::numeric_limits<double>::max();                               
      double local_max_x = std::numeric_limits<double>::min();
    //   MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
    //   MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
        
        // std::vector<std::vector<double>compositional_values
        // (this->n_compositional_fields(),
        // std::vector<double> (quadrature.size()));

        std::vector<double>compositional_values(face_corners.size()); 
        // double compositional_values =0.0;
        // std::vector<std::vector<double>> temporary_variables(0, std::vector<double>());

      //std::cout<<"name composition selected :"<<this->introspection().name_for_compositional_index(composition_extent)<<std::endl;
        
        const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
//         const bool top_boundary = boundary_ids.find(relevant_boundary) == boundary_ids.begin();

      // loop over all cells and find cell with 100 % of the defined composition, then save the elevation to stored_value
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
        {
            bool composition_extent_tracked_present = false;   

          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
              if (cell->face(face_no)->at_boundary())
              {
                if (cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;
                
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.compositional_fields[composition_extent]].get_function_values (this->get_solution(),
                          compositional_values);
                
                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                      const Point<dim> vertex = fe_face_values.quadrature_point(corner);
                      

                    // std::cout<<corner<<"  "<<compositional_values[corner]<<std::endl;

                  // for (unsigned int p=0; p<face_corners.size(); ++p)
                  //   {     
                      if (compositional_values[corner] >= 0.95)
                        {
                            // track the composition at the surface
//                             if(vertex(1) > extents[1]-composition_depth)
//                             {
                          composition_extent_tracked_present = true;
//                             }
                        }                        
                    //  } 

               if(composition_extent_tracked_present)
               {
                    if (vertex(0) < local_min_x)
                    {
                      local_min_x = vertex(0);
                    }
                    if (vertex(0) > local_max_x)
                    {
                      local_max_x = vertex(0);
                    }                    
                    
               } 
                }           
          }        
        }
      }
      }

      //Calculate max depth across all processes
      const double min_x = std::abs(Utilities::MPI::min(local_min_x, this->get_mpi_communicator())); 
      const double max_x = std::abs(Utilities::MPI::max(local_max_x, this->get_mpi_communicator()));
      
      double extent_composition = max_x-min_x;
//       std::cout<<"min_x :"<<min_x<<" max_x :"<<max_x<<std::endl;    

      // statistics.add_value ("Maximum depth value for composition " + this->introspection().name_for_compositional_index(composition_extent),
      //                           local_min_height);

      
//       double a_dt = this->get_timestep();
//       if (this->convert_output_to_years())
//       {
//         a_dt = this->get_timestep() / year_in_seconds;
//       }
      
      //Write results to statistics file
      statistics.add_value ("Composition extent (m)",
                            extent_composition);
      const char *columns[] = { "Composition extent (m)"};
      for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
        {
          statistics.set_precision (columns[i], 8);
          statistics.set_scientific (columns[i], true);
        }

      output_stats.precision(5);
      output_stats << extent_composition << " m";

      // Just return stats if text output is not required at all or not needed at this time
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return std::pair<std::string, std::string> (" Extent of the selected composition (m)",
                                                    output_stats.str());


      return std::pair<std::string, std::string> ("Composition extent (m)",
                                                  output_stats.str());

    }      
  


    template <int dim>
    void ExtentComposition<dim>::declare_parameters(ParameterHandler &prm)
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
        prm.enter_subsection("Extent composition");
        {        
            prm.declare_entry("Composition tracked for extent", "1",
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
          prm.declare_entry ("Output extent dependant ending file", "false",
                             Patterns::List(Patterns::Bool()),
                             "Output a Termin aspect file when a depth that user set has been reached by the choosen composition"
                             ", for example, in a subduction model, when the slab reach a defined depth"
                             ", a fill is outputed and using the request user termination criterion, the simulation can stop");   
          prm.declare_entry ("Ending extent criterion", "300000",
                             Patterns::Double (0.),
                             "The depth choosen to ouput a random file that cam be used as an termine aspect file"); 
          prm.declare_entry ("Composition maximum depth", "5000",
                             Patterns::Double (0.),
                             "Maximum depth until which the composition extent sould be tracked");           
        }
          prm.leave_subsection();        
      }
      prm.leave_subsection();    
    }                           

    template <int dim>
    void
    ExtentComposition<dim>::parse_parameters (ParameterHandler &prm)
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
        prm.enter_subsection("Extent composition");
        {
          composition_extent = prm.get_integer("Composition tracked for extent");
          composition_depth = prm.get_double ("Composition maximum depth");
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
    ASPECT_REGISTER_POSTPROCESSOR(ExtentComposition,
                                  "Composition extent",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
