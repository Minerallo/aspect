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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

// #include <aspect/mesh_refinement/isotherms.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/derivative_approximation.h>
#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>
 
namespace aspect
{
  namespace MeshRefinement
  {


    /**
     * A class that implements a mesh refinement criterion based on the
     * compositional fields (if available).
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class Subduction : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        tag_additional_cells () const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double lithosphere_zone_refined;
        double weakzone_two_upper_depth;
        double weakzone_two_lower_depth;
        double weakzone_one_lower_depth;
        double ridge_two_upper_depth;
        double ridge_two_lower_depth;
        double ridge_one_lower_depth;       
        double crust_zone_refined;
        double oceanic_domain_refined;
        double continental_mantle_refined;
        double abs_plate_velocity;
        /**
         * The compositional field number, min ref level and max ref level
         * for the crust, mantle part of the slab lithosphere and the
         * overriding plate
         */
        std::vector<unsigned int> sediments_refinement;         
        std::vector<unsigned int> upper_crust_refinement;
        std::vector<unsigned int> upper_crust_refinement_two;        
        std::vector<unsigned int> lower_crust_refinement;
        std::vector<unsigned int> lower_crust_refinement_two;        
        std::vector<unsigned int> continental_mantle_refinement;
        std::vector<unsigned int> craton_refinement;
        std::vector<unsigned int> weak_zone_refinement;
        std::vector<unsigned int> weak_zone_refinement_two;
        std::vector<unsigned int> ridge_refinement;
        std::vector<unsigned int> ridge_refinement_two;        
        std::vector<unsigned int> oceanic_crust_refinement; 
        std::vector<unsigned int> hazburgite_mantle_refinement;        
        std::vector<unsigned int> oceanic_mantle_refinement; 
        std::vector<unsigned int> mantle_refinement;
        std::vector<unsigned int> upper_mantle_refinement;
        // std::vector<unsigned int> border_refinement_level;
        /**
         * The absolute minimum refinement level
         */
        unsigned int min_level;
        unsigned int max_level;
        /**
         * The absolute maximum refinement level
         */
    };



    template <int dim>
    void
    Subduction<dim>::tag_additional_cells () const
    {
      if (this->get_dof_handler().n_locally_owned_dofs() == 0)
        return;        

      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).get_unit_support_points());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      QTrapez<dim-1> face_corners;

      // const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).get_unit_support_points());
        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                         this->get_fe(),
                                         face_corners,
                                         update_values |
                                             update_quadrature_points);                               

    //   MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
    //   MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
        
        std::vector<std::vector<double> > prelim_composition_values
        (this->n_compositional_fields(),
        std::vector<double> (quadrature.size()));
        
      // the velocity is given in cm/yr,a the new eastern limit 
      double crust_zone_refined_updated = crust_zone_refined-(this->get_time()/ year_in_seconds)*abs_plate_velocity*0.01;
      std::cout<<crust_zone_refined_updated<<std::endl; 

      for (typename DoFHandler<dim>::active_cell_iterator
           cell = this->get_dof_handler().begin_active();
           cell != this->get_dof_handler().end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
                bool coarsen = false;
                bool refine = false;
                bool clear_refine = false;
                bool clear_coarsen = false;
                bool upper_crust_present = false;
                bool upper_crust_present_two = false;
                bool lower_crust_present = false;
                bool lower_crust_present_two = false;
                bool sediments_present =false;
                bool continental_mantle_present = false;
                bool craton_present = false;
                // bool oceanic_mantle_border_present = false;
                bool ridge_present = false;
                bool ridge_present_two = false;
                bool oceanic_crust_present = false;
                bool oceanic_crust_present_one = false;                
                bool oceanic_mantle_present = false;
                bool hazburgite_mantle_present = false;                
                bool weak_zone_present = false;
                bool weak_zone_present_two = false;                
                // bool UPM_present = false;
                // bool TZ_present = false;
                // bool TZ2_present = false;
                // bool LM_present = false;
                // bool overriding_present = false;
//                 bool in_center_of_compo = false;
//                 bool refine_border_present = false;
                bool mantle_present =false;
                bool upper_mantle_present =false;
                bool box_out=false; 
                bool oceanic_mantle_out=false;
                
              fe_values.reinit(cell);
            //   in.reinit(fe_values, cell, this->introspection(), this->get_solution(), true);
            //   this->get_material_model().evaluate(in, out);                        

                for (unsigned int c = 0; c<this->n_compositional_fields(); c++)
                  {
                    fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values (this->get_solution(),
                        prelim_composition_values[c]);

                // if (refine_border[c] == 1){
                // // if the composition exceeds the threshold, cell is marked for refinement
                // for (unsigned int j=0; j<this->get_fe().base_element(this->introspection().base_elements.compositional_fields).dofs_per_cell; ++j){
                //   if (prelim_composition_values[c][j] > 0.35 && prelim_composition_values[c][j] < 0.65)
                //     {
                //       refine_border_present = true;
                //       break;
                //     }                  
                //     }            
                //   }
                } 

                for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                {
                  fe_face_values.reinit(cell, face_no);

                  for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                    {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);                


                    for (unsigned int p=0; p<quadrature.size(); ++p)
                      {
                        //  if (prelim_composition_values[weak_layer_refinement[0]][p] > 0.01)
                        //   {
                        //     weak_layer_present = true;

                        //   }     
                        // std::cout<<prelim_composition_values[upper_crust_refinement[0]][p]<<std::endl;
                          
                        
                         if (prelim_composition_values[ridge_refinement[0]][p] >= 0.60)
                          {
                            if (vertex(0) < crust_zone_refined_updated+300000){  
                            if(vertex(1) > lithosphere_zone_refined)
                            {
                              if(ridge_two_lower_depth <=vertex(1) && vertex(1)<=ridge_two_upper_depth){  
//                               if (vertex(1) > 1050000){
                              ridge_present_two = true;
                              // if (prelim_composition_values[weak_zone_refinement[0]][p] >= 1.0)
                              // {
                              //   in_center_of_compo = true;
                              // }
                              break;
//                               }
                              }else if (ridge_one_lower_depth <=vertex(1) && vertex(1)<ridge_two_lower_depth){  
                                ridge_present = true;
                                break;
                              }
                              // }
                              }
                        }else{
                            box_out=true;
                        }                              
                            
                          }
                        if (prelim_composition_values[weak_zone_refinement[0]][p] >= 0.60)
                          {
                            if (vertex(0) < crust_zone_refined_updated+300000){
                            if(vertex(1) > lithosphere_zone_refined)
                            {
                              if(weakzone_two_lower_depth <=vertex(1) && vertex(1)<=weakzone_two_upper_depth){  
//                               if (vertex(1) > 1050000){
                              weak_zone_present_two = true;
                              // if (prelim_composition_values[weak_zone_refinement[0]][p] >= 1.0)
                              // {
                              //   in_center_of_compo = true;
                              // }
                              break;
//                               }
                              }else if (weakzone_one_lower_depth <=vertex(1) && vertex(1)<weakzone_two_lower_depth){  
                                weak_zone_present = true;
                                break;
                              }
                              // }
                              }
                         }else{
                            box_out=true;
                        }                             
                            
                          }  

                        
                          
                          
                        if (vertex(1) >= continental_mantle_refined)
                        { 
                            if (prelim_composition_values[upper_crust_refinement[0]][p] > 0.01)
                            {
                             if (vertex(0) < crust_zone_refined_updated){  
                               upper_crust_present = true; 
                             }else{
                                 box_out=true;
                            }
                            }
                          }else if(vertex(1) < continental_mantle_refined)
                          {
                           if (prelim_composition_values[upper_crust_refinement[0]][p] > 0.5)
                           {
                            if (vertex(0) < crust_zone_refined_updated){   
                            upper_crust_present = true;
                                
                            }else{
                                box_out=true;
                            }
                           }
                        }

                        
                      
                        
                        if (vertex(1) >= continental_mantle_refined)
                        { 
                            if (prelim_composition_values[lower_crust_refinement[0]][p] > 0.01)
                            {
                               if (vertex(0) < crust_zone_refined_updated){
                               lower_crust_present = true; 
                               }else{
                            box_out = true;
                            }                            
                            }
                          }else if(vertex(1) < continental_mantle_refined)
                          {
                           if (prelim_composition_values[lower_crust_refinement[0]][p] > 0.5)
                           {
                               if (vertex(0) < crust_zone_refined_updated){
                            lower_crust_present = true; 
                               }else{
                                   box_out=true;
                            }
                        }
                        }
                      
                        if (prelim_composition_values[sediments_refinement[0]][p] > 0.01)
                          {
                            if(vertex(1) > lithosphere_zone_refined)
                            {
                                if (vertex(0) < crust_zone_refined_updated){
                            sediments_present = true;}
                            else{
                                box_out=true;
                            }
                            // break;
                            }
                          }
                        if (prelim_composition_values[oceanic_crust_refinement[0]][p] > 0.1)
                          {
                            if(vertex(1) > lithosphere_zone_refined)
                            { 
                              if (vertex(0) < crust_zone_refined_updated+300000){  
                              if (vertex(0) > oceanic_domain_refined){                            
                              oceanic_crust_present = true;
                              break;                           
                              }else{
                                oceanic_crust_present_one = true;                            
                              }
                              }else{
                                  box_out=true;
                            }
                            }
                          }
                          

                        
                        if (prelim_composition_values[oceanic_mantle_refinement[0]][p] > 0.1)
                        {
                            if(vertex(1) > lithosphere_zone_refined)
                            {
                            if (vertex(0) < crust_zone_refined_updated+600000){    
                            oceanic_mantle_present = true;
                            break;
                            }else{
                            oceanic_mantle_out=true;     
                            }
//                             if (prelim_composition_values[oceanic_mantle_refinement[0]][p] >= 1.0)
//                             {
//                               in_center_of_compo = true;
//                             }
                            // if (prelim_composition_values[oceanic_mantle_refinement[0]][p] < 0.1)
                            // // if(prelim_composition_values[oceanic_mantle_refinement[0]][p] > 0.35 && prelim_composition_values[oceanic_mantle_refinement[0]][p] < 0.65)
                            // {
                            //   refine_border_present = true;
                            // }
                            }   
                        }
                        if (prelim_composition_values[hazburgite_mantle_refinement[0]][p] > 0.1)
                        {
                            if(vertex(1) > lithosphere_zone_refined)
                            {
                            if (vertex(0) < crust_zone_refined_updated+600000){    
                            hazburgite_mantle_present = true;
                            break; 
                            }else{
                            oceanic_mantle_out=true;   
                            }
//                             if (prelim_composition_values[hazburgite_mantle_refinement[0]][p] >= 1.0)
//                             {
//                               in_center_of_compo = true;
//                             }
                            // if (prelim_composition_values[oceanic_mantle_refinement[0]][p] < 0.1)
                            // // if(prelim_composition_values[oceanic_mantle_refinement[0]][p] > 0.35 && prelim_composition_values[oceanic_mantle_refinement[0]][p] < 0.65)
                            // {
                            //   refine_border_present = true;
                            // }
                            }   
                        } 
                        

                        
                        
                        if (prelim_composition_values[continental_mantle_refinement[0]][p] > 0.5)
                          {
                            if(vertex(1) > lithosphere_zone_refined)
                            {                                   
                              if (vertex(1) > continental_mantle_refined){                            
                              continental_mantle_present = true;
                              // if (prelim_composition_values[continental_mantle_refinement[0]][p] >= 1.0)
                              // {
                              //   in_center_of_compo = true;
                              // }
                              }
                            }
                          }
                        if (prelim_composition_values[craton_refinement[0]][p] > 0.1)
                          {
                            if(vertex(1) > lithosphere_zone_refined)
                            {                                   
                            craton_present = true;
//                             if (prelim_composition_values[craton_refinement[0]][p] >= 1.0)
//                             {
//                               in_center_of_compo = true;
//                             }
//                             }
                          }
                          }
                        if (prelim_composition_values[upper_mantle_refinement[0]][p] > 0.80)
                          {
                            if(vertex(1) > lithosphere_zone_refined)
                            {                                   
                            upper_mantle_present = true;
                            }
                          }                          
                        if (prelim_composition_values[0][p] > 0.85)
                          {
                            mantle_present = true;
                          }                          
                          
                                                                                                                                                          
                      }
                  }
                }
        
               //Only continue if at least one is true

                    int refinement_level = min_level;
                    
                    if(oceanic_mantle_out){
                        refinement_level =1;
                    }
                    else if(box_out){
                        refinement_level = 2;
                    }
                    else if (upper_crust_present)
                      {
                        // std::cout<<"UPC"<<std::endl; 
                        refinement_level =  upper_crust_refinement[1];

                      }
                    else if (upper_crust_present_two)
                      {
                        // std::cout<<"UPC"<<std::endl; 
                        refinement_level =  upper_crust_refinement_two[1];

                      }                      
                    else if (lower_crust_present)
                      {
                        // std::cout<<"UPL"<<std::endl; 
                        refinement_level =  lower_crust_refinement[1];
                
                      }  
                    else if (lower_crust_present_two)
                      {
                        // std::cout<<"UPL"<<std::endl; 
                       refinement_level =  lower_crust_refinement_two[1];
                  
                      }                        
                    else if (sediments_present)
                      {
                        refinement_level =  sediments_refinement[1];
               
                      }
                    // else if (weak_layer_present)
                    //   {
                    //     minimum_refinement_level = weak_layer_refinement[1];
                    //     maximum_refinement_level = weak_layer_refinement[2];
                    //   }
                    else if (oceanic_crust_present)
                      {
                        refinement_level = oceanic_crust_refinement[1];
                     
                      }                        
                    else if (oceanic_crust_present_one)
                      {
                        refinement_level =  oceanic_crust_refinement[1]-1;
                
                      }                         
                    else if (oceanic_mantle_present)
                      {
                        // std::cout<<"OM"<<std::endl; 
                        refinement_level =  oceanic_mantle_refinement[1];
              
                        // if (in_center_of_compo)
                        // minimum_refinement_level = oceanic_mantle_refinement[1];
                        // maximum_refinement_level = oceanic_mantle_refinement[1];
                        // std::cout<<oceanic_mantle_refinement[0]+1<<std::endl;
                        // std::cout<<refine_border[oceanic_mantle_refinement[0]+1]<<std::endl;
                        // if(refine_border[oceanic_mantle_refinement[0]+1]==1){
                      }
                    else if (hazburgite_mantle_present)
                      {
                        // std::cout<<"OM"<<std::endl; 
                        refinement_level =  hazburgite_mantle_refinement[1];
              
                        // if (in_center_of_compo)
                        // minimum_refinement_level = oceanic_mantle_refinement[1];
                        // maximum_refinement_level = oceanic_mantle_refinement[1];
                        // std::cout<<oceanic_mantle_refinement[0]+1<<std::endl;
                        // std::cout<<refine_border[oceanic_mantle_refinement[0]+1]<<std::endl;
                        // if(refine_border[oceanic_mantle_refinement[0]+1]==1){
                      }                      
                    // else if(oceanic_mantle_present == true  && refine_border_present == true )
                    //       {
                    //         minimum_refinement_level = border_refinement_level;
                    //         maximum_refinement_level = border_refinement_level;                            
                    //       }                          
                        // }
                      // }
                    else if (weak_zone_present)
                      {
                        refinement_level =  weak_zone_refinement[1];
             
                      }    
                    else if (weak_zone_present_two)
                      {
                        refinement_level =  weak_zone_refinement_two[1];
                
                      }
                    else if (ridge_present)
                      {
                        refinement_level =  ridge_refinement[1];
             
                      }    
                    else if (ridge_present_two)
                      {
                       refinement_level =  ridge_refinement_two[1];
                  
                      }                       
                    else if (continental_mantle_present)
                      {
                        // std::cout<<"CM"<<std::endl; 
                        refinement_level =  continental_mantle_refinement[1];
            
                      }
                     else if (craton_present)
                      {
                        refinement_level =  craton_refinement[1];
              
                      }          
                      // This mantle present will be helpfull to not refine tiny amount of composition that get lost in the model
                    else if (upper_mantle_present)
                      {
                        refinement_level =  mantle_refinement[1];
                 
                      }                                                                                                                                                                       
                    else
                      {
                        refinement_level =  mantle_refinement[0];
                
                      }

                    
                    // if (refine_border_present){
                    //   maximum_refinement_level = border_refinement_level;
                    // }   

                    const int cell_level = cell->level();
                    if (cell_level >= refinement_level)
                      {
                        clear_refine = false;
                      }
                    if (cell_level >  refinement_level)
                      {
                        coarsen = true;
                      }
                    if (cell_level <= refinement_level)
                      {
                        clear_coarsen = true;
                      }
                    if (cell_level < refinement_level)
                      {
                        refine = true;
                      }   
              // if both coarsen and refine are true, give preference to refinement
              if (coarsen == true && refine == true)
                {
                  coarsen = false;
                  clear_refine = false;
                }

              // Perform the actual placement of the coarsening and refinement flags
              // We want to make sure that the refiment never goes below the minimum
              // or above the maximum, so we first check/set the coarsen/refine flag,
              // and then check/set the clear coarsen/refine flag.
              if (coarsen == true)
                {
                  cell->set_coarsen_flag ();
                }
              if (clear_coarsen == true)
                {
                  cell->clear_coarsen_flag ();
                }
              if (refine == true)
                {
                  cell->set_refine_flag ();
                }
              if (clear_refine == true)
                {
                  cell->clear_refine_flag ();
                }
              }
        }
    }

    template <int dim>
    void
    Subduction<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Subduction");
        {
          prm.declare_entry("Sediments refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");            
          prm.declare_entry("Upper Crust refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Upper Crust 2 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");                            
          prm.declare_entry("Lower Crust refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Lower Crust 2 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");                            
          prm.declare_entry("Continental Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");  
          prm.declare_entry("Craton refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");                                                           
          // prm.declare_entry("Weak Layer refinement","",
          //                   Patterns::List (Patterns::Integer(0)),
          //                   "The compositional field number of the crust, its minimum refinement level and "
          //                   "its maximum refinement level.");
          prm.declare_entry("Oceanic Crust refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level."); 
          prm.declare_entry("Oceanic Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Hazburgite Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");          
          prm.declare_entry("Weak Zone refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");  
          prm.declare_entry("Weak Zone 2 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Ridge refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");  
          prm.declare_entry("Ridge 2 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");  
          prm.declare_entry("Upper Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");          
          prm.declare_entry("Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Refine border","",
                            Patterns::List (Patterns::Double()),
                            "List of 1 and 0, indicating if the compositions border need or not to be refined.");                              
          prm.declare_entry("Refine border level","",
                            Patterns::List (Patterns::Integer(0)),
                            "The refinement level for border should be only 1 value for now.");
         prm.declare_entry ("Apply refinements if on top of", "250000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine weak zone 2 if under", "1090000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine weak zone 2 if on top of", "1050000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine weak zone 1 if on top of", "980000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$"); 
         prm.declare_entry ("Refine ridge 2 if under", "1090000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine ridge 2 if on top of", "1050000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine ridge 1 if on top of", "980000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");         
         prm.declare_entry ("Refine continental domain if x inferior at", "1670000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$"); 
         prm.declare_entry ("Refine oceanic domain if x superior at", "350000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine continental mantle if on top of", "800000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");  
         prm.declare_entry ("Overriding plate velocity", "2", Patterns::Double(0),
                             "The velocity is given in cm/yr. Units: $1 / s$");        
         
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Subduction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        min_level = prm.get_integer("Minimum refinement level");
        max_level = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
        prm.enter_subsection("Subduction");
        {
        lithosphere_zone_refined = prm.get_double("Apply refinements if on top of"); 
        weakzone_two_upper_depth = prm.get_double("Refine weak zone 2 if under"); 
        weakzone_two_lower_depth = prm.get_double("Refine weak zone 2 if on top of"); 
        weakzone_one_lower_depth = prm.get_double("Refine weak zone 1 if on top of");
        ridge_two_upper_depth = prm.get_double("Refine ridge 2 if under"); 
        ridge_two_lower_depth = prm.get_double("Refine ridge 2 if on top of"); 
        ridge_one_lower_depth = prm.get_double("Refine ridge 1 if on top of");        
        crust_zone_refined = prm.get_double("Refine continental domain if x inferior at");
        oceanic_domain_refined = prm.get_double("Refine oceanic domain if x superior at");
        continental_mantle_refined = prm.get_double("Refine continental mantle if on top of");
        abs_plate_velocity = prm.get_double("Overriding plate velocity");

        
//           const std::vector<double> refine_border
//             = Utilities::string_to_double(
//                 Utilities::split_string_list(prm.get("Refine border")));
// 
//           AssertThrow (refine_border.size() == this->n_compositional_fields(),
//                        ExcMessage ("The number of refine border given here must be "
//                                    "equal to the number of compositionnal field."));
//          
//          border_refinement_level = prm.get_integer("Refine border level");

          // const std::vector<int> border_refinement_level
          //   = Utilities::string_to_int(
          //       Utilities::split_string_list(prm.get("Refine border level")));
          
          // border_refinement_level = std::vector<unsigned int> (border_refinement_level.begin(),border_refinement_level.end());

          // AssertThrow (border_refinement_level.size() == 1,
          //              ExcMessage ("Should be one level value for all fields for now"));


          const std::vector<int> sediments
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Sediments refinement")));

          sediments_refinement = std::vector<unsigned int> (sediments.begin(),sediments.end());

          AssertThrow (sediments_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (sediments_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (sediments_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));




          const std::vector<int> upper_crust
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Upper Crust refinement")));

          upper_crust_refinement = std::vector<unsigned int> (upper_crust.begin(),upper_crust.end());

          AssertThrow (upper_crust_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (upper_crust_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (upper_crust_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));


          const std::vector<int> upper_crust_two
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Upper Crust 2 refinement")));

          upper_crust_refinement_two = std::vector<unsigned int> (upper_crust_two.begin(),upper_crust_two.end());

          AssertThrow (upper_crust_refinement_two.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (upper_crust_refinement_two[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (upper_crust_refinement_two[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

                                  

          const std::vector<int> lower_crust
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Lower Crust refinement")));

          lower_crust_refinement = std::vector<unsigned int> (lower_crust.begin(),lower_crust.end());

          AssertThrow (lower_crust_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (lower_crust_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (lower_crust_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));



          const std::vector<int> lower_crust_two
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Lower Crust 2 refinement")));

          lower_crust_refinement_two = std::vector<unsigned int> (lower_crust_two.begin(),lower_crust_two.end());

          AssertThrow (lower_crust_refinement_two.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (lower_crust_refinement_two[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (lower_crust_refinement_two[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));


          const std::vector<int> continental_mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Continental Mantle refinement")));

          continental_mantle_refinement = std::vector<unsigned int> (continental_mantle.begin(),continental_mantle.end());

          AssertThrow (continental_mantle_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (continental_mantle_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (continental_mantle_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));




          const std::vector<int> craton
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Craton refinement")));

          craton_refinement = std::vector<unsigned int> (craton.begin(),craton.end());

          AssertThrow (craton_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (craton_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (craton_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));




          // const std::vector<int> weak_layer
          //   = Utilities::string_to_int(
          //       Utilities::split_string_list(prm.get("Weak Layer refinement")));

          // weak_layer_refinement = std::vector<unsigned int> (weak_layer.begin(),weak_layer.end());

          // AssertThrow (weak_layer_refinement.size() == 3,
          //              ExcMessage ("The number of refinement data given here must be "
          //                          "equal to 3 (field number + min level + max level). "));

          // AssertThrow (weak_layer_refinement[0] < this->n_compositional_fields(),
          //              ExcMessage ("The number of compositional field to refine (starting "
          //                          "from 0) should be smaller than the number of fields. "));

          // AssertThrow (weak_layer_refinement[1] >= min_level,
          //              ExcMessage ("The minimum refinement for the crust cannot be "
          //                          "smaller than the minimum level of the whole model. "));

          // AssertThrow (weak_layer_refinement[2] <= max_level,
          //              ExcMessage ("The maximum refinement for the crust cannot be "
          //                          "greater than the maximum level of the whole model. "));


          const std::vector<int> oceanic_crust
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Oceanic Crust refinement")));

          oceanic_crust_refinement = std::vector<unsigned int> (oceanic_crust.begin(),oceanic_crust.end());

          AssertThrow (oceanic_crust_refinement.size() ==2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (oceanic_crust_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (oceanic_crust_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));




          const std::vector<int> oceanic_mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Oceanic Mantle refinement")));

          oceanic_mantle_refinement = std::vector<unsigned int> (oceanic_mantle.begin(),oceanic_mantle.end());

          AssertThrow (oceanic_mantle_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (oceanic_mantle_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (oceanic_mantle_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          
          const std::vector<int> hazburgite_mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Hazburgite Mantle refinement")));

          hazburgite_mantle_refinement = std::vector<unsigned int> (hazburgite_mantle.begin(),hazburgite_mantle.end());

          AssertThrow (hazburgite_mantle_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (hazburgite_mantle_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (hazburgite_mantle_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));
         

          const std::vector<int> ridge
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Ridge refinement")));

          ridge_refinement = std::vector<unsigned int> (ridge.begin(),ridge.end());

          AssertThrow (ridge_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (ridge_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (ridge_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));


          

          const std::vector<int> ridge_two
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Ridge 2 refinement")));

          ridge_refinement_two = std::vector<unsigned int> (ridge_two.begin(),ridge_two.end());

          AssertThrow (ridge_refinement_two.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (ridge_refinement_two[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (ridge_refinement_two[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          // AssertThrow (weak_zone_refinement_two[2] <= max_level,
          //              ExcMessage ("The maximum refinement for the crust cannot be "
          //                          "greater than the maximum level of the whole model. "));
          
          const std::vector<int> weak_zone
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Weak Zone refinement")));

          weak_zone_refinement = std::vector<unsigned int> (weak_zone.begin(),weak_zone.end());

          AssertThrow (weak_zone_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (weak_zone_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (weak_zone_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          // AssertThrow (weak_zone_refinement[2] <= max_level,
          //              ExcMessage ("The maximum refinement for the crust cannot be "
          //                          "greater than the maximum level of the whole model. "));
          

          const std::vector<int> weak_zone_two
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Weak Zone 2 refinement")));

          weak_zone_refinement_two = std::vector<unsigned int> (weak_zone_two.begin(),weak_zone_two.end());

          AssertThrow (weak_zone_refinement_two.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (weak_zone_refinement_two[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (weak_zone_refinement_two[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          // AssertThrow (weak_zone_refinement_two[2] <= max_level,
          //              ExcMessage ("The maximum refinement for the crust cannot be "
          //                          "greater than the maximum level of the whole model. "));
          

          const std::vector<int> upper_mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Upper Mantle refinement")));

          upper_mantle_refinement = std::vector <unsigned int> (upper_mantle.begin(),
                                                          upper_mantle.end());

          AssertThrow (upper_mantle_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 2 (min level + max level). "));

          AssertThrow (weak_zone_refinement_two[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (weak_zone_refinement_two[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));
          
          const std::vector<int> mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Mantle refinement")));

          mantle_refinement = std::vector <unsigned int> (mantle.begin(),
                                                          mantle.end());

          AssertThrow (mantle_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 2 (min level + max level). "));

          AssertThrow (mantle_refinement[0] >= min_level,
                       ExcMessage ("The minimum refinement for the mantle cannot be "
                                   "smaller than the minimum level of the whole model. "));


        //   AssertThrow (crust_refinement[0] != slab_mantle_refinement[0] && \
        //                crust_refinement[0] != overriding_refinement[0]  && \
        //                slab_mantle_refinement[0] != overriding_refinement[0], 
        //                ExcMessage ("Defined refinement fields the same. "));

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
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Subduction,
                                              "subduction",
                                              "subduction refinement")
  }
}
