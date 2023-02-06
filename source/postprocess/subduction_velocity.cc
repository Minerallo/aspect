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
    class SubductionVelocity : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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

        double vel_area_x;
        double vel_area_y;         
    };

    template <int dim>
    std::pair<std::string,std::string>
    SubductionVelocity<dim>::execute (TableHandler &statistics)
    {
      QTrapez<dim-1> face_corners;

        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                         this->get_fe(),
                                         face_corners,
                                         update_values |
                                             update_quadrature_points);

      std::ostringstream output_stats;

         std::vector<double> temporary_vel(face_corners.size());
         double mean_vel_proc =0.0;
         double return_value = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
        {
            bool area_tracked_present = false;   

          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
                std::vector<Tensor<1, dim>> vel(face_corners.size());
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel); 

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                  const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                    if(vertex(0) > vel_area_x-5000 && vertex(0) < vel_area_x+5000 && vertex(1) > extents[1]-vel_area_y-5000 && vertex(1) < extents[1]-vel_area_y+5000) 
                    {
                      area_tracked_present = true; 
                      std::cout<<"added"<<std::endl;
                    }
                    if(area_tracked_present)
                    {
                    temporary_vel.push_back(vel[corner][0] * year_in_seconds);                     
                     return_value = return_value + vel[corner][0] * year_in_seconds;
                     mean_vel_proc=return_value/temporary_vel.size();
                    }                                                         

            }
          }        
        }
      }

      //Calculate max velocity across all processes
      const double mean_vel = Utilities::MPI::max(mean_vel_proc, this->get_mpi_communicator());   


      //Write results to statistics file
      statistics.add_value ("Oce_vel(m/yr)",
                            mean_vel);
      const char *columns[] = { "Oce_vel(m/yr)"};
      for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
        {
          statistics.set_precision (columns[i], 8);
          statistics.set_scientific (columns[i], true);
        }

      output_stats.precision(4);
      output_stats << mean_vel << " m/yr";

      // Just return stats if text output is not required at all or not needed at this time
      if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
                             && (this->get_timestep_number() != 0)))
        return std::pair<std::string, std::string> (" Oce_vel(m/yr) :",
                                                    output_stats.str());


      return std::pair<std::string, std::string> ("Oce_vel(m/yr):",
                                                  output_stats.str());
    }      
  


    template <int dim>
    void SubductionVelocity<dim>::declare_parameters(ParameterHandler &prm)
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
        prm.enter_subsection("Subduction Velocity");
        {        
          prm.declare_entry ("Mean velocity area in x", "100000",
                             Patterns::Double (0.),
                             "The position in x choosen to calculate an average of the velocity");     
          prm.declare_entry ("Mean velocity area in y", "100000",
                             Patterns::Double (0.),
                             "The position in depth choosen to calculate an average of the velocity");                                                                    
        }
          prm.leave_subsection();        
      }
      prm.leave_subsection();    
    }                           

    template <int dim>
    void
    SubductionVelocity<dim>::parse_parameters (ParameterHandler &prm)
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
        prm.enter_subsection("Subduction Velocity");
        {
          vel_area_x       = prm.get_double ("Mean velocity area in x");   
          vel_area_y       = prm.get_double ("Mean velocity area in y");                                               
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
    ASPECT_REGISTER_POSTPROCESSOR(SubductionVelocity,
                                  "subduction velocity",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
