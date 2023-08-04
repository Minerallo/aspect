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
        
         double Z_zone_refined; 
         double X_two_zone_refined;
         double X_one_zone_refined;
         double Y_one_zone_refined;
         double Y_two_zone_refined;
        
        
        /**
         * The compositional field number, min ref level and max ref level
         * for the crust, compo_else part of the slab lithosphere and the
         * overriding plate
         */
        
        std::vector<unsigned int> compo_one_refinement;         
        std::vector<unsigned int> compo_two_refinement;
        std::vector<unsigned int> compo_three_refinement;
        std::vector<unsigned int> compo_four_refinement;
        std::vector<unsigned int> compo_five_refinement;        
        std::vector<unsigned int> compo_else_refinement;       
        std::vector<unsigned int> interface_refinement;
   

        // std::vector<unsigned int> border_refinement_level;
        /**
         * The absolute minimum refinement level
         */
        unsigned int min_level;
        /**
         * The absolute maximum refinement level
         */
        unsigned int max_level;

//         std::vector<double> refine_border;
//         unsigned int border_refinement_level;
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
                bool compo_one_present = false;
                bool compo_two_present = false;
                bool compo_three_present =false;
                bool compo_four_present = false;
                bool compo_five_present = false;                
                bool compo_else_present =false;
                bool interface_present = false;


              fe_values.reinit(cell);                     

                for (unsigned int c = 0; c<this->n_compositional_fields(); c++)
                  {
                    fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values (this->get_solution(),
                        prelim_composition_values[c]);

                } 

                for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                {
                  fe_face_values.reinit(cell, face_no);

                  for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                    {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);                


                    for (unsigned int p=0; p<quadrature.size(); ++p)
                      {
                           if(vertex(0)>X_one_zone_refined && vertex(0)<X_two_zone_refined && vertex(1)>Y_one_zone_refined && vertex(1)<Y_two_zone_refined && vertex(2) > Z_zone_refined)
                            {
                         if (prelim_composition_values[interface_refinement[0]][p] >= 0.005)
                          {    
                                interface_present = true;
                                break;
                              
                          }                         
                        if (prelim_composition_values[compo_one_refinement[0]][p] > 0.01)
                          {                         
    
                              compo_one_present = true;
   
                          }
                        if (prelim_composition_values[compo_two_refinement[0]][p] > 0.01)
                          {
                                    
                            compo_two_present = true;

                            
                          }
                        if (prelim_composition_values[compo_three_refinement[0]][p] > 0.01)
                          {
                                   
                            compo_three_present = true;
                            
                          }
                        if (prelim_composition_values[compo_four_refinement[0]][p] > 0.01)
                          {
                                        
                            compo_four_present = true;
                                                                                                    
                            }
                        if (prelim_composition_values[compo_five_refinement[0]][p] > 0.01)
                          {
                                        
                            compo_five_present = true;
                                                                                                    
                            }                            
                            }
                            else 
                          {
                            compo_else_present = true;
                          }
                      
                            
                                                                                                                                                          
                      }
                  }
                }
        
               //Only continue if at least one is true

                    int refinement_level = min_level;              

                    if (compo_one_present)
                      {
                        refinement_level = compo_one_refinement[1];
                      }                    
                    else if (compo_two_present)
                      {
                        refinement_level = compo_two_refinement[1];

                      }                       
                    else if (compo_three_present)
                      {
                        refinement_level = compo_three_refinement[1];

                      }
                     else if (compo_four_present)
                      {
                        refinement_level = compo_four_refinement[1];
                      }  
                     else if (compo_five_present)
                      {
                        refinement_level = compo_five_refinement[1];
                      }                       
                    else if (interface_present)
                      {
                        refinement_level = interface_refinement[1];
                      }                                                                                                                                                                                                     
                    else if (compo_else_present)
                      {
                        refinement_level = compo_else_refinement[0];
                      }
                    

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
           prm.declare_entry("Interface refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");                 
          prm.declare_entry("Compo 1 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");  
            prm.declare_entry("Compo 2 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");           
          prm.declare_entry("Compo 3 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level."); 
          prm.declare_entry("Compo 4 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level."); 
          prm.declare_entry("Compo 5 refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level."); 
          prm.declare_entry("Compo else refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");          
         prm.declare_entry ("Apply refinements if on top of", "250000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine domain if x inferior at", "1670000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$"); 
         prm.declare_entry ("Refine domain if x superior at", "350000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
         prm.declare_entry ("Refine domain if y inferior at", "1670000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$"); 
         prm.declare_entry ("Refine domain if y superior at", "350000", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");       
         
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
        Z_zone_refined = prm.get_double("Apply refinements if on top of");       
        X_two_zone_refined = prm.get_double("Refine domain if x inferior at");
        X_one_zone_refined = prm.get_double("Refine domain if x superior at");
        Y_one_zone_refined = prm.get_double("Refine domain if y superior at");
        Y_two_zone_refined = prm.get_double("Refine domain if y inferior at");
    
        


          const std::vector<int> compo_one
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Compo 1 refinement")));

          compo_one_refinement = std::vector<unsigned int> (compo_one.begin(),compo_one.end());

          AssertThrow (compo_one_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (compo_one_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (compo_one_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));


                    const std::vector<int> compo_two
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Compo 2 refinement")));

          compo_two_refinement = std::vector<unsigned int> (compo_two.begin(),compo_two.end());

          AssertThrow (compo_two_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (compo_two_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (compo_two_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));
          
          const std::vector<int> compo_three
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Compo 3 refinement")));

          compo_three_refinement = std::vector<unsigned int> (compo_three.begin(),compo_three.end());

          AssertThrow (compo_three_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (compo_three_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (compo_three_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));  
          
          const std::vector<int> compo_four
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Compo 4 refinement")));

          compo_four_refinement = std::vector<unsigned int> (compo_four.begin(),compo_four.end());

          AssertThrow (compo_four_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (compo_four_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (compo_four_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. ")); 
          
           const std::vector<int> compo_five
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Compo 5 refinement")));

          compo_five_refinement = std::vector<unsigned int> (compo_four.begin(),compo_four.end());

          AssertThrow (compo_five_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (compo_five_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (compo_five_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));          
          
          
           const std::vector<int> compo_else
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Compo else refinement")));

          compo_else_refinement = std::vector<unsigned int> (compo_else.begin(),compo_else.end());

          AssertThrow (compo_else_refinement.size() == 1,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (compo_else_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));
      
          
            
          const std::vector<int> interface
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Interface refinement")));

          interface_refinement = std::vector<unsigned int> (interface.begin(),interface.end());

          AssertThrow (interface_refinement.size() == 2,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (interface_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (interface_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));




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