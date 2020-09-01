// /*
//   Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

//   This file is part of ASPECT.

//   ASPECT is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2, or (at your option)
//   any later version.

//   ASPECT is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with ASPECT; see the file LICENSE.  If not see
//   <http://www.gnu.org/licenses/>.
// */

// #include <aspect/postprocess/interface.h>

// #include <deal.II/base/quadrature_lib.h>
// #include <deal.II/fe/fe_values.h>
// #include <aspect/geometry_model/box.h>
// #include <aspect/simulator.h>
// #include <deal.II/base/quadrature_lib.h>


// #include <boost/archive/text_oarchive.hpp>
// #include <boost/archive/text_iarchive.hpp>


// #include <aspect/postprocess/interface.h>
// #include <aspect/simulator_access.h>

// #include <cmath>
// #include <limits>

// #include <aspect/geometry_model/two_merged_boxes.h>

// namespace aspect
// {
//   namespace Postprocess
//   {

//     /**
//      * A postprocessor that computes some statistics about the compositional
//      * fields, if any.
//      *
//      * @ingroup Postprocessing
//      */
//     template <int dim>
//     class Trench : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
//     {
//       public:
//         /**
//          * Evaluate the solution for some temperature statistics.
//          */
//         std::pair<std::string,std::string> execute (TableHandler &statistics) override;

//          /**
//          * Parse parameters for the free surface handling.
//          */
//         void parse_parameters (ParameterHandler &prm);    

//       private:

//        /**
//          * Extent of the whole model domain in x-, y-, and z-direction (in 3d).
//          */
//         Point<dim> extents;  
//         // double position_x ;    
//         // double position_y ;                     
//     };

//     template <int dim>
//     std::pair<std::string,std::string>
//     Trench<dim>::execute (TableHandler &statistics)
//     {
//       QTrapez<dim-1> face_corners;

//       const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).get_unit_support_points());
//         FEFaceValues<dim> fe_face_values(this->get_mapping(),
//                                          this->get_fe(),
//                                          face_corners,
//                                          update_values |
//                                              update_quadrature_points);
//       // FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, quadrature, update_quadrature_points | update_values | update_gradients);

//       // have a stream into which we write the data. the text stream is then
//       // later sent to processor 0
//       std::ostringstream output_stats;
//       // std::ostringstream output_file;


//       // Choose stupidly large values for initialization
//       double local_min_depth = std::numeric_limits<double>::max();          

//     //   MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
//     //   MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
        
//         // std::vector<std::vector<double>compositional_values
//         // (this->n_compositional_fields(),
//         // std::vector<double> (quadrature.size()));

//         // std::vector<std::vector<double>> temporary_variables(0, std::vector<double>());



//       // loop over all cells and find cell with 100 % of the defined composition, then save the elevation to stored_value
//       for (const auto &cell : this->get_dof_handler().active_cell_iterators())
//       {
//         if (cell->is_locally_owned())
//         {  

//           for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
//           {
//                 fe_face_values.reinit(cell, face_no);

//                 for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
//                   {
//                   const Point<dim> vertex = fe_face_values.quadrature_point(corner);
//                     const double elevation = this->get_geometry_model().height_above_reference_surface(vertex);  

//                     if (elevation < local_min_depth)
//                       local_min_depth = elevation;                   
//                }

//             }
//           }        
//         }
//       }


//                     //   position_x = vertex(0);      
//                     //   position_y = vertex(1);
//       //Calculate max depth across all processes
//       const double min_depth_trench = extents[1] - Utilities::MPI::min(local_min_depth, this->get_mpi_communicator());                  

//       // std::cout<<"Maximum depth value for composition :"<<local_min_height<<std::endl;

//       // statistics.add_value ("Maximum depth value for composition " + this->introspection().name_for_compositional_index(composition_number),
//       //                           local_min_height);

//       //Write results to statistics file
//       statistics.add_value ("trench depth (m)",
//                             min_depth_trench);
//     //   statistics.add_value ("trench_x (m)",
//     //                         position_x);
//     //   statistics.add_value ("trench_y (m)",
//     //                         position_y);                            
//       const char *columns[] = { "trench depth (m)"
//     //   ,
//     //                             "trench_x (m)",
//     //                             "trench_y (m)"
//                               };

//       for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
//         {
//           statistics.set_precision (columns[i], 8);
//           statistics.set_scientific (columns[i], true);
//         }

//       output_stats.precision(4);
//       output_stats << min_depth_trench << " m";
//     //   output_stats << position_x << " m";    
//     //   output_stats << position_y << " m";           

//       // Just return stats if text output is not required at all or not needed at this time
//     //   if (!write_to_file || ((this->get_time() < last_output_time + output_interval)
//     //                          && (this->get_timestep_number() != 0)))
//     //     return std::pair<std::string, std::string> (" Trench depth and position (m):",
//     //                                                 output_stats.str());


//       return std::pair<std::string, std::string> ("Trench depth and position (m):",
//                                                   output_stats.str());

//     }

//     template <int dim>
//     void
//     Trench<dim>::parse_parameters (ParameterHandler &prm)
//     {
//       prm.enter_subsection("Geometry model");
//       {
//         prm.enter_subsection("Box with lithosphere boundary indicators");
//         {
//           // Total box extents
//           extents[0]           = prm.get_double ("X extent");
//           extents[1]           = prm.get_double ("Y extent");
//                   }prm.leave_subsection();
//       }prm.leave_subsection(); 
//     }  

//   } 
// }    


// // explicit instantiations
// namespace aspect
// {
//   namespace Postprocess
//   {
//     ASPECT_REGISTER_POSTPROCESSOR(Trench,
//                                   "trench",
//                                   "A postprocessor that computes some statistics about "
//                                   "the compositional fields, if present in this simulation. "
//                                   "In particular, it computes maximal and minimal values of "
//                                   "each field, as well as the total mass contained in this "
//                                   "field as defined by the integral "
//                                   "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
//   }
// }
