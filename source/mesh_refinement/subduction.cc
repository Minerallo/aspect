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
        /**
         * The compositional field number, min ref level and max ref level
         * for the crust, mantle part of the slab lithosphere and the
         * overriding plate
         */
        std::vector<unsigned int> weak_zone_refinement;         
        std::vector<unsigned int> oceanic_crust_refinement;         
        std::vector<unsigned int> oceanic_mantle_refinement; 
        std::vector<unsigned int> mantle_refinement;

        /**
         * The absolute minimum refinement level
         */
        unsigned int min_level;
        /**
         * The absolute maximum refinement level
         */
        unsigned int max_level;
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
                bool oceanic_crust_present = false;
                bool oceanic_crust_present_one = false;                
                bool oceanic_mantle_present = false;                
                bool weak_zone_present = false;
                bool weak_zone_present_two = false;                
                bool mantle_present =false;


              fe_values.reinit(cell);
            //   in.reinit(fe_values, cell, this->introspection(), this->get_solution(), true);
            //   this->get_material_model().evaluate(in, out);                        

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
                        if (prelim_composition_values[weak_zone_refinement[0]][p] >= 0.005)
                          {
                            if(vertex(1) > 250000)
                            {
                              if (vertex(1) > 980000){
                              weak_zone_present = true;
                              break;
                              }
                            }
                          }                                      
                        if (prelim_composition_values[oceanic_crust_refinement[0]][p] > 0.1)
                          {
                            if(vertex(1) > 250000)
                            {                                   
                              if (vertex(0) > 350000){                            
                              oceanic_crust_present = true;
                              break;                           
                              }else{
                                oceanic_crust_present_one = true;                            
                              }
                            }
                          }                                                                             
                        if (prelim_composition_values[oceanic_mantle_refinement[0]][p] > 0.1)
                        {
                            if(vertex(1) > 250000)
                            {                                 
                            oceanic_mantle_present = true;
                            break; 
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

                    int maximum_refinement_level = max_level;
                    int minimum_refinement_level = min_level;               


                    if (oceanic_crust_present)
                      {
                        minimum_refinement_level = oceanic_crust_refinement[1];
                        maximum_refinement_level = oceanic_crust_refinement[2];
                      }   
                    else if (oceanic_crust_present_one)
                      {
                        minimum_refinement_level = oceanic_crust_refinement[1]-1;
                        maximum_refinement_level = oceanic_crust_refinement[2]-1;
                      }                         
                    else if (oceanic_mantle_present)
                      {
                        minimum_refinement_level = oceanic_mantle_refinement[1];
                        maximum_refinement_level = oceanic_mantle_refinement[2];
                      }
                    else if (weak_zone_present)
                      {
                        minimum_refinement_level = weak_zone_refinement[1];
                        maximum_refinement_level = weak_zone_refinement[2];
                      }    
                                          
                    else if (mantle_present)
                      {
                        minimum_refinement_level = mantle_refinement[0];
                        maximum_refinement_level = mantle_refinement[1];
                      }                                                                                                                                                                       
                    else
                      {
                        minimum_refinement_level = mantle_refinement[0];
                        maximum_refinement_level = mantle_refinement[1];
                      }



                    const int cell_level = cell->level();
                    if (cell_level >= maximum_refinement_level)
                      {
                        clear_refine = false;
                      }
                    if (cell_level >  maximum_refinement_level)
                      {
                        coarsen = true;
                      }
                    if (cell_level <= minimum_refinement_level)
                      {
                        clear_coarsen = true;
                      }
                    if (cell_level < minimum_refinement_level)
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
          prm.declare_entry("Oceanic Crust refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level."); 
          prm.declare_entry("Oceanic Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");
          prm.declare_entry("Weak Zone refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");                                                                                                                                                  
          prm.declare_entry("Mantle refinement","",
                            Patterns::List (Patterns::Integer(0)),
                            "The compositional field number of the crust, its minimum refinement level and "
                            "its maximum refinement level.");                          
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

          const std::vector<int> oceanic_crust
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Oceanic Crust refinement")));

          oceanic_crust_refinement = std::vector<unsigned int> (oceanic_crust.begin(),oceanic_crust.end());

          AssertThrow (oceanic_crust_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (oceanic_crust_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (oceanic_crust_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (oceanic_crust_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the crust cannot be "
                                   "greater than the maximum level of the whole model. "));


          const std::vector<int> oceanic_mantle
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Oceanic Mantle refinement")));

          oceanic_mantle_refinement = std::vector<unsigned int> (oceanic_mantle.begin(),oceanic_mantle.end());

          AssertThrow (oceanic_mantle_refinement.size() == 3,
                       ExcMessage ("The number of refinement data given here must be "
                                   "equal to 3 (field number + min level + max level). "));

          AssertThrow (oceanic_mantle_refinement[0] < this->n_compositional_fields(),
                       ExcMessage ("The number of compositional field to refine (starting "
                                   "from 0) should be smaller than the number of fields. "));

          AssertThrow (oceanic_mantle_refinement[1] >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (oceanic_mantle_refinement[2] <= max_level,
                       ExcMessage ("The maximum refinement for the crust cannot be "
                                   "greater than the maximum level of the whole model. "));


          const std::vector<int> weak_zone
            = Utilities::string_to_int(
                Utilities::split_string_list(prm.get("Weak Zone refinement")));

          weak_zone_refinement = std::vector<unsigned int> (weak_zone.begin(),weak_zone.end());

          AssertThrow (weak_zone_refinement.size() == 3,
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

          AssertThrow (mantle_refinement[1] <= max_level,
                       ExcMessage ("The maximum refinement for the mantle cannot be "
                                   "greater than the maximum level of the whole model. "));

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
