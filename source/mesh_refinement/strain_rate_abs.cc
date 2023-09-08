/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


// #include <aspect/mesh_refinement/strain_rate.h>

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
     * A class that implements a mesh refinement criterion based on
     * the strain rate field.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class StrainRateAbs : public Interface<dim>,
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
          bool use_superior;
          double strain_rate_cut;              
          double strain_rate_level;
          double lim_inf;
          double lim_sup;
        /**
         * The absolute minimum refinement level
         */
        unsigned int min_level;
        /**
         * The absolute maximum refinement level
         */
        unsigned int max_level;  
        double start_strain_refinement;        

      };





    template <int dim>
    void
    StrainRateAbs<dim>::tag_additional_cells() const
    {

      // double a_dt = this->get_timestep();
      // double b_dt = this->get_old_timestep(); 
      // //I need it to get a better approximation of the real time but this doesn not really make sense
      // double time =this->get_time();
      // if (this->convert_output_to_years())
      // {
      //   a_dt = this->get_timestep() / year_in_seconds;
      //   b_dt = this->get_old_timestep() / year_in_seconds;
      //   time = this->get_time()/year_in_seconds;
      // }
      // double full_time = time+a_dt+b_dt/2;//not right but I don't find the right time for some reason. time + a_dt should be enough.
      // // std::cout<<"refinement, model time : "<<a_dt<<"  "<<this->get_time()<<"  "<<this->get_timestep()<<"  "<<this->get_old_timestep()<<std::endl;

      if (this->get_time()/ year_in_seconds > start_strain_refinement){
      // std::cout<<a_dt<<std::endl; 
        // std::cout<<"refinement, model b time : "<<full_time<<std::endl; 
        // std::cout<<"finally working"<<std::endl;
      const QMidpoint<dim> quadrature;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      std::vector<SymmetricTensor<2,dim> > strain_rates (quadrature.size());

            QTrapez<dim-1> face_corners;

      // const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).get_unit_support_points());
        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                         this->get_fe(),
                                         face_corners,
                                         update_values |
                                             update_quadrature_points);  


      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
                bool coarsen = false;
                bool refine = false;
                bool clear_refine = false;
                bool clear_coarsen = false;
                bool smaller_strain_rate =  false; 
                // bool bigger_strain_rate = false ; 

            // const unsigned int idx = cell->active_cell_index();
            fe_values.reinit(cell);

            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                strain_rates);

               for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                {
                  fe_face_values.reinit(cell, face_no);

                  for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                    {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);                

                    for (unsigned int p=0; p<quadrature.size(); ++p)
                      {

                        if (lim_inf <=vertex(1) && vertex(1)<=lim_sup){                        
                          if(use_superior){
                            if (strain_rate_cut<= strain_rates[0].norm()){
                                smaller_strain_rate = true;
                            }                  
                          }else {
                            if (strain_rate_cut>= strain_rates[0].norm()){
                                smaller_strain_rate = true;
                            }                         
                          }
                        } 
                      }
                    }
                }               

            // else{
            //     bigger_strain_rate = true;
            // }


            // int maximum_refinement_level = max_level;
            int minimum_refinement_level = min_level;

            if (smaller_strain_rate){
                minimum_refinement_level = strain_rate_level;
            }            

            const int cell_level = cell->level();
            if (cell_level >= minimum_refinement_level)
                {
                clear_refine = false;
                }
            if (cell_level >  minimum_refinement_level)
                {
                coarsen = false;
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
      }else{
        std::cout<<"not yet"<<std::endl;
      }
    }

    template <int dim>
    void
    StrainRateAbs<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Strain rate cut");
        {
          prm.declare_entry ("Start refine with strain rate", "1000", Patterns::Double(0),
                             "One may want to start to refine later to avoid to over refine due to the isostatic reequilibirum. Units: year");           
          prm.declare_entry("Refine strain superior", "true",
                            Patterns::Bool(),
                            "Flag to refine superior or lower strain");          
          prm.declare_entry ("Strain rate threshold", "3e-13", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");          
          prm.declare_entry ("Strain rate minimum refinement level", "4", Patterns::Integer(0),
                             "The level of refinement of the mesh, cannot be more than the maximum level of refinement of the model");  
          prm.declare_entry ("Interval min", "1050", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");
          prm.declare_entry ("Interval max", "1095", Patterns::Double(0),
                             "the strain rate at which the mesh should start to be refined. Units: $1 / s$");              
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    StrainRateAbs<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        min_level = prm.get_integer("Minimum refinement level");
        max_level = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
        prm.enter_subsection("Strain rate cut");
        {
          start_strain_refinement = prm.get_double("Start refine with strain rate");
          // start_strain_refinement *= year_in_seconds;
          use_superior = prm.get_bool("Refine strain superior");          
          strain_rate_cut = prm.get_double("Strain rate threshold");              
          strain_rate_level = prm.get_integer("Strain rate minimum refinement level");
          lim_inf = prm.get_double("Interval min");  
          lim_sup = prm.get_double("Interval max");  
          AssertThrow (strain_rate_level >= min_level,
                       ExcMessage ("The minimum refinement for the crust cannot be "
                                   "smaller than the minimum level of the whole model. "));

          AssertThrow (strain_rate_level <= max_level,
                       ExcMessage ("The maximum refinement for the crust cannot be "
                                   "greater than the maximum level of the whole model. "));

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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(StrainRateAbs,
                                              "strain_rate_cut",
                                              "A usefull mesh criterion based on a threshold,This will help to refine some specific strained area ")
  }
}