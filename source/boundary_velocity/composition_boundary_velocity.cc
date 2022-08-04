// /*
//   Copyright (C) 2011 - 2021 by the authors of the ASPECT code.
// 
//   This file is part of ASPECT.
// 
//   ASPECT is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2, or (at your option)
//   any later version.
// 
//   ASPECT is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
// 
//   You should have received a copy of the GNU General Public License
//   along with ASPECT; see the file LICENSE.  If not see
//   <http://www.gnu.org/licenses/>.
// */
// 
// 
// #include <aspect/utilities.h>
// #include <aspect/global.h>
// #include <deal.II/base/signaling_nan.h>
// 
// #include <deal.II/base/quadrature_lib.h>
// #include <deal.II/fe/fe_values.h>
// #include <aspect/geometry_model/box.h>
// #include <aspect/simulator.h>
// #include <deal.II/base/quadrature_lib.h>
// #include <aspect/simulator_access.h>
// #include <cmath>
// #include <limits>
// #include <aspect/geometry_model/two_merged_boxes.h>
// 
// #include <aspect/boundary_velocity/interface.h>
// #include <deal.II/base/parsed_function.h>
// 
// namespace aspect
// {
//   namespace BoundaryVelocity
//   {
//     using namespace dealii;
// 
//     /**
//      * A class that implements velocity boundary conditions based on a
//      * functional description provided in the input file.
//      *
//      * @ingroup BoundaryVelocities
//      */
//     template <int dim>
//     class Composition_boundary_velocity : public Interface<dim>, public SimulatorAccess<dim>
//     {
//       public:
//         /**
//          * Constructor.
//          */
//         Composition_boundary_velocity ();
// 
//         /**
//          * Return the boundary velocity as a function of position. For the
//          * current class, this function obviously simply returns a zero
//          * tensor.
//          */
//         Tensor<1,dim>
//         boundary_velocity (const types::boundary_id boundary_indicator,
//                            const Point<dim> &position) const override;
// 
//         // avoid -Woverloaded-virtual warning until the deprecated function
//         // is removed from the interface:
//         using Interface<dim>::boundary_velocity;
// 
//         /**
//          * A function that is called at the beginning of each time step to
//          * indicate what the model time is for which the boundary values will
//          * next be evaluated. For the current class, the function passes to
//          * the parsed function what the current time is.
//          */
//         void
//         update () override;
// 
//         /**
//          * Declare the parameters this class takes through input files. The
//          * default implementation of this function does not describe any
//          * parameters. Consequently, derived classes do not have to overload
//          * this function if they do not take any runtime parameters.
//          */
//         static
//         void
//         declare_parameters (ParameterHandler &prm);
// 
//         /**
//          * Read the parameters this class declares from the parameter file.
//          * The default implementation of this function does not read any
//          * parameters. Consequently, derived classes do not have to overload
//          * this function if they do not take any runtime parameters.
//          */
//         void
//         parse_parameters (ParameterHandler &prm) override;
// 
//       private:
//         /**
//          * A function object representing the components of the velocity.
//          */
//         Composition_boundary_velocity::ParsedFunction<dim> boundary_velocity_function;
// 
//         /**
//          * The coordinate representation to evaluate the function. Possible
//          * choices are depth, cartesian and spherical.
//          */
//         Utilities::Coordinates::CoordinateSystem coordinate_system;
// 
//         /**
//          * Whether to specify velocity in x, y, z components, or
//          * r, phi, theta components.
//          */
//         bool use_spherical_unit_vectors;
//     };
// 
// 
//     template <int dim>
//     Composition_boundary_velocity<dim>::Composition_boundary_velocity ()
//       :
//       boundary_velocity_function (dim)
//     {}
// 
// 
// 
//     template <int dim>
//     Tensor<1,dim>
//     Composition_boundary_velocity<dim>::
//     boundary_velocity (const types::boundary_id ,
//                        const Point<dim> &position) const
//     {
//       Tensor<1,dim> velocity;
// 
//       Utilities::NaturalCoordinate<dim> point =
//         this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);
//      
//       QTrapez<dim-1> face_corners;
// 
//       FEFaceValues<dim> fe_face_values(this->get_mapping(),
//                                          this->get_fe(),
//                                          face_corners,
//                                          update_values |
//                                             update_quadrature_points);
//       
//       // Choose stupidly large values for initialization
//       double local_min_x = std::numeric_limits<double>::max();                               
//       double local_max_x = std::numeric_limits<double>::min();
//    
//         std::vector<double>compositional_values(face_corners.size()); 
//         
//     const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");
//     
//      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
//       {
//         if (cell->is_locally_owned())
//         {
//             bool composition_tracked = false;   
// 
//           for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
//           {
//               if (cell->face(face_no)->at_boundary())
//               {
//                 if (cell->face(face_no)->boundary_id() != relevant_boundary)
//                   continue;
//                 
//                 fe_face_values.reinit(cell, face_no);
//                 fe_face_values[this->introspection().extractors.compositional_fields[composition_extent]].get_function_values (this->get_solution(),
//                           compositional_values);
//                 
//                 for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
//                   {
//                       const Point<dim> vertex = fe_face_values.quadrature_point(corner);
//                                
//                       if (compositional_values[corner] >= 0.95)
//                         {
//                            
//                           composition_tracked = true;
//                         
//                         }                        
// 
//                if(composition_extent_tracked_present)
//                {
//                     if (vertex(0) < local_min_x)
//                     {
//                       local_min_x = vertex(0);
//                     }
//                     if (vertex(0) > local_max_x)
//                     {
//                       local_max_x = vertex(0);
//                     }  
//                     
//                 }         
//               }           
//             }        
//           }
//         }
//       }
//       Calculate max depth across all processes
//       const double min_x = std::abs(Utilities::MPI::min(local_min_x, this->get_mpi_communicator())); 
//       const double max_x = std::abs(Utilities::MPI::max(local_max_x, this->get_mpi_communicator()));
//       double extent_composition = max_x-min_x;
//       
//         for (const auto &cell : this->get_dof_handler().active_cell_iterators())
//       {
//         if (cell->is_locally_owned())
//         {
//             bool composition_tracked = false;   
// 
//           for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
//           {
//               if (cell->face(face_no)->at_boundary())
//               {
//                 if (cell->face(face_no)->boundary_id() != relevant_boundary)
//                   continue;
//                 
//                 fe_face_values.reinit(cell, face_no);
//                 fe_face_values[this->introspection().extractors.compositional_fields[composition_extent]].get_function_values (this->get_solution(),
//                           compositional_values);
//                 
//                 for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
//                   {
//                       const Point<dim> vertex = fe_face_values.quadrature_point(corner);
//       
//         for (unsigned int d=0; d<dim; ++d)
//         {
//             if (relevant_boundary && )
//             {
//         velocity[d] = boundary_velocity_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), d);
//             }else{
//        velocity[d] = boundary_velocity_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), d);
// 
//             }
//         }
// 
//       if (use_spherical_unit_vectors)
//         velocity = Utilities::Coordinates::spherical_to_cartesian_vector(velocity, position);
// 
//       // ASPECT always wants things in MKS system. however, as described
//       // in the documentation of this class, we interpret the formulas
//       // given to this plugin as meters per year if the global flag
//       // for using years instead of seconds is given. so if someone
//       // write "5" in their parameter file and sets the flag, then this
//       // means "5 meters/year" and we need to convert it to the ASPECT
//       // time system by dividing by the number of seconds per year
//       // to get MKS units
//       if (this->convert_output_to_years())
//         return velocity / year_in_seconds;
//       else
//         return velocity;
//       
//     }
//       
//     template <int dim>
//     void
//     Composition_boundary_velocity<dim>::update()
//     {
//       // we get time passed as seconds (always) but may want
//       // to reinterpret it in years
//       if (this->convert_output_to_years())
//         boundary_velocity_function.set_time (this->get_time() / year_in_seconds);
//       else
//         boundary_velocity_function.set_time (this->get_time());
//     }
// 
// 
//     template <int dim>
//     void
//     Composition_boundary_velocity<dim>::declare_parameters (ParameterHandler &prm)
//     {
//       prm.enter_subsection("Boundary composition velocity model");
//       {
//         prm.enter_subsection("Function");
//         {
//           prm.declare_entry("Composition tracked", "1",
//                               Patterns::Integer(),
//                               "composition number that is tracked. ");        
//           prm.declare_entry ("Coordinate system", "cartesian",
//                              Patterns::Selection ("cartesian|spherical|depth"),
//                              "A selection that determines the assumed coordinate "
//                              "system for the function variables. Allowed values "
//                              "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
//                              "are interpreted as r,phi or r,phi,theta in 2D/3D "
//                              "respectively with theta being the polar angle. `depth' "
//                              "will create a function, in which only the first "
//                              "parameter is non-zero, which is interpreted to "
//                              "be the depth of the point.");
//           prm.declare_entry ("Use spherical unit vectors", "false",
//                              Patterns::Bool (),
//                              "Specify velocity as $r$, $\\phi$, and $\\theta$ components "
//                              "instead of $x$, $y$, and $z$. Positive velocities point up, east, "
//                              "and north (in 3D) or out and clockwise (in 2D). "
//                              "This setting only makes sense for spherical geometries.");
// 
//           Composition_boundary_velocitys::ParsedFunction<dim>::declare_parameters (prm, dim);
//         }
//         prm.leave_subsection();
//       }
//       prm.leave_subsection();
//     }
// 
// 
// 
//     template <int dim>
//     void
//     Composition_boundary_velocity<dim>::parse_parameters (ParameterHandler &prm)
//     {
//       prm.enter_subsection("Boundary composition velocity model");
//       {
//         prm.enter_subsection("Function");
//         {
//           coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
//           use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
//           if (use_spherical_unit_vectors)
//             AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
//                          ExcMessage ("Spherical unit vectors should not be used "
//                                      "when geometry model is not spherical."));
//         }
//         try
//           {
//             boundary_velocity_function.parse_parameters (prm);
//           }
//         catch (...)
//           {
//             std::cerr << "ERROR: FunctionParser failed to parse\n"
//                       << "\t'Boundary velocity model.Function'\n"
//                       << "with expression\n"
//                       << "\t'" << prm.get("Function expression") << "'"
//                       << "More information about the cause of the parse error \n"
//                       << "is shown below.\n";
//             throw;
//           }
//         prm.leave_subsection();
//       }
//       prm.leave_subsection();
//     }
// 
//   }
// }
// 
// // explicit instantiations
// namespace aspect
// {
//   namespace BoundaryVelocity
//   {
//     ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(Composition_boundary_velocity,
//                                             "composition boundary velocity",
//                                             "Implementation of a model in which the boundary "
//                                             "velocity is given in terms of an explicit formula "
//                                             "that is elaborated in the parameters in section "
//                                             "``Boundary velocity model|Function''. The format of these "
//                                             "functions follows the syntax understood by the "
//                                             "muparser library, see Section~\\ref{sec:muparser-format}."
//                                             "\n\n"
//                                             "The formula you describe in the mentioned "
//                                             "section is a semicolon separated list of velocities "
//                                             "for each of the $d$ components of the velocity vector. "
//                                             "These $d$ formulas are interpreted as having units "
//                                             "m/s, unless the global input parameter ``Use "
//                                             "years in output instead of seconds'' is set, in "
//                                             "which case we interpret the formula expressions "
//                                             "as having units m/year."
//                                             "\n\n"
//                                             "Likewise, since the symbol $t$ indicating time "
//                                             "may appear in the formulas for the prescribed "
//                                             "velocities, it is interpreted as having units "
//                                             "seconds unless the global parameter above has "
//                                             "been set.")
//   }
// }
