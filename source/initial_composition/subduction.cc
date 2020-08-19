/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.
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


//#include "2Dsubd_compo.h"
#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/boundary_temperature/box.h>

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{

  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the porosity field
     * by computing the equilibrium melt fraction for the given initial
     * condition and reference pressure profile. Note that this plugin only
     * works if there is a compositional field called 'porosity', and the
     * used material model implements the 'MeltFractionModel' interface.
     * All compositional fields except porosity are not changed by this plugin.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class  Subduction : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position,
                                    const unsigned int n_comp) const override;
        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        double hdepth;
        double Xres;
        // double Yres=1100000;
        double border_left;
        double border_right;
        double xtrench;
        double xsedim;
        double xcraton;
        double xforearc;
        double dforearc;// how thicker is the continental lithospheric mantle at the forearc
        double dsedim;
        double lsedim;
        double dweak;
        double dcrust1;
        double dcrust2;
        double dlwcrust1;
        double dlwcrust2;
        double dmantcont;
        double dcraton;
        double dmantoce;
        double dgabbro;
        double dsticky;
        double TZ1;
        double TZ2; 
        double TZ3; 

        // double dlithoforearc = dcrust1+dlwcrust1+dmantcont+dsticky+dforearc;
        double dlithocont;
        double dlithocraton;
        double hcrust3;
        double hlwcrust3;

        double angleforearc;
        double angle4;
        // double Ts=293;
        // double tlitho1=1601;
        // double tlitho2=1613;
        // double tlitho=1573;
        // double Tp=1610;

        double weakbottom;
        double segment;  
        double angle;

        //Find the radius by rearranging Arclength = (2* numbers::PI*r)/(360/angle)
        double radius;
        double radiusweak;
        double radiusgabbro;
        double radiusmantoce;
        double rpyweak;
        double rpygabbro;
        double rpymantoce;

        //set it so the bounds of the circle are at the surface and the
        //center is aligned with the trench location.

        double rpx = xtrench;

        //Determine the y location of the end of the segment
        double segment_y;
        double tip;

                   
    };
  
  // namespace InitialComposition
  // {


    // template <int dim>
    // void
    //   Subduction<dim>::initialize ()
    //  {}


    template <int dim>
    double
    Subduction<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      // UPM,TZ1,TZ2,LM,Weak_Layer,Oceanic_crust,Oceanic_mantle,Sediments,Upper_Crust,Lower_Crust,Continental_mantle,Craton,plastic_strain
      const unsigned int UPM_index = this->introspection().compositional_index_for_name("UPM");
      const unsigned int TZ1_index = this->introspection().compositional_index_for_name("TZ");
      const unsigned int TZ2_index = this->introspection().compositional_index_for_name("TZ2");
      const unsigned int LM_index = this->introspection().compositional_index_for_name("LM");
      const unsigned int WKL_index = this->introspection().compositional_index_for_name("Weak_Layer");
      const unsigned int OC_index = this->introspection().compositional_index_for_name("Oceanic_crust");
      const unsigned int OM_index = this->introspection().compositional_index_for_name("Oceanic_mantle");
      const unsigned int Sed_index = this->introspection().compositional_index_for_name("Sediments");
      const unsigned int UPC_index = this->introspection().compositional_index_for_name("Upper_Crust");
      const unsigned int LC_index = this->introspection().compositional_index_for_name("Lower_Crust");
      const unsigned int CM_index = this->introspection().compositional_index_for_name("Continental_mantle");
      const unsigned int Craton_index = this->introspection().compositional_index_for_name("Craton");
      // const unsigned int SA_index = this->introspection().compositional_index_for_name("Sticky_Air");                


      // //Sticky Air
      // if(position[1]>hdepth-dsticky)
      //         {
      //           return (n_comp == SA_index ) ? 1 : 0;
      //         }
      //Weak layer
      // else if (position[0]>=border_left && position[0]<=xtrench && position[1]>hdepth-dsticky-dweak)
      if (position[0]>=border_left && position[0]<=xtrench && position[1]>hdepth-dsticky-dweak)
              {
                return (n_comp == WKL_index ) ? 1 : 0;
              }
      //Weak layer
      else if (sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpyweak),2)) <= radius+weakbottom && sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpyweak),2)) > radiusweak &&  \
                position[0] > xtrench && position[1] >= segment_y && position[1]>=hdepth-dgabbro-dmantoce-dweak-tip)
              {
                return (n_comp == WKL_index ) ? 1 : 0;  
              }
      // Oceanic_crust
      else if (position[0]>=border_left && position[0]<=xtrench && position[1]<=hdepth-dsticky-dweak && position[1]>=hdepth-dgabbro-dweak-dsticky)
              {
                return (n_comp == OC_index ) ? 1 : 0;
              }
      // Oceanic_crust
      else if (sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpygabbro),2)) <= radius && sqrt( pow((position[0]-rpx),2) + pow((position[1]-rpygabbro),2)) > radiusgabbro &&  \
                position[0] > xtrench && position[1] >= segment_y && position[1]>=hdepth-dgabbro-dmantoce-dweak-tip)
              {
                return (n_comp == OC_index ) ? 1 : 0;           
              }
      // lithosphere oceanic
      else if (position[0]>=border_left && position[0]<=xtrench && position[1]<hdepth-dsticky-dweak-dgabbro && position[1]>=hdepth-dgabbro-dmantoce-dweak-dsticky)
              {
                return (n_comp == OM_index ) ? 1 : 0;        
              }
      else if (sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpymantoce),2)) <= radius && sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpymantoce),2)) >= radiusmantoce &&  \
                position[0] > xtrench && position[1] >= segment_y && position[1]>=hdepth-dgabbro-dmantoce-dweak-tip)
              {
                return (n_comp == OM_index ) ? 1 : 0;            
              }        
      //Continental sediments
      else if(position[0]>=xsedim && position[0]<=xsedim+lsedim && position[1]>=hdepth-dsedim-dsticky)
              {
                return (n_comp == Sed_index ) ? 1 : 0;           
              } 
      //Upper_crust
      else if (sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpyweak),2)) > radius+weakbottom && position[0]> xtrench && position[1]>=hdepth-dcrust1-dsticky && position[0]<xsedim || position[0]>=xsedim && position[0]<=xsedim+lsedim && position[1]<hdepth-dsedim-dsticky && position[1]>=hdepth-dsedim-dcrust2-dsticky || position[0]>xsedim+lsedim && position[0]<=border_right && position[1]>=hdepth-hcrust3-dsticky)
              {
                return (n_comp == UPC_index ) ? 1 : 0;        
              }  
      //Lower_crust
      else if(sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpyweak),2)) > radius+weakbottom && position[0]> xtrench && position[1]<hdepth-dcrust1-dsticky && position[1]>=hdepth-dcrust1-dlwcrust1-dsticky && position[0]<xsedim || position[0]>=xsedim && position[0]<=xsedim+lsedim && position[1]<hdepth-dsedim-dcrust2-dsticky && position[1]>=hdepth-dsedim-dcrust2-dlwcrust2-dsticky || position[0]>xsedim+lsedim && position[0]<=border_right && position[1]<hdepth-hcrust3-dsticky && position[1]>=hdepth-hlwcrust3-dsticky )
              {
                return (n_comp == LC_index ) ? 1 : 0;           
              }  
      //Litho mantle continenatl (Forearc)
      else if(sqrt( pow((position[0]-rpx),2) + pow((position[1]-rpyweak),2)) > radius+weakbottom && position[0]> xtrench && position[1]<hdepth-dcrust1-dlwcrust1-dsticky && position[1]>=(hdepth-dsticky-dcrust1-dlwcrust1-dmantcont+tan(angleforearc*( numbers::PI/180))*(-xforearc+position[0])) && position[0]<=xforearc)
              {
                return (n_comp == CM_index ) ? 1 : 0;            
              }  
      else if(position[0]>=xforearc && position[0]<xcraton && position[1]<hdepth-dcrust1-dlwcrust1-dsticky && position[1]>=hdepth-dcrust1-dlwcrust1-dmantcont-dsticky)  
              {
                return (n_comp == CM_index ) ? 1 : 0;           
              }  
      //Craton
      else if(position[0]>=xcraton && position[0]<=xcraton+((dlithocraton-dlithocont)/tan(angle4* numbers::PI/180)) && position[1]>=(hdepth-dsticky-dcrust1-dlwcrust1-dmantcont-tan(angle4*( numbers::PI/180))*(-xcraton+position[0])) && position[1]<hdepth-dcrust2-dlwcrust2-dsticky || position[0]>xcraton+((dlithocraton-dlithocont)/tan(angle4* numbers::PI/180)) && position[1]<hdepth-dsedim-dcrust2-dlwcrust2-dsticky && position[1]>=hdepth-dsedim-dcrust2-dlwcrust2-dcraton-dsticky && position[0]<=border_right)
              {
                return (n_comp == Craton_index ) ? 1 : 0;           
              }  
      // Asthenosphere Transition zone
      else if(position[1]<=hdepth-TZ1 && position[1]>hdepth-TZ2)
              {
                return (n_comp == TZ1_index ) ? 1 : 0;           
              }  
      else if(position[1]<=hdepth-TZ2 && position[1]>hdepth-TZ3)
              {
                return (n_comp == TZ2_index ) ? 1 : 0;             
              }  
      else if(position[1]<=hdepth-TZ3)
              {
                return (n_comp == LM_index ) ? 1 : 0;              
              }  
      else
              {
                return (n_comp == UPM_index ) ? 1 : 0;
              }
        
    }
    

    template <int dim>
    void Subduction<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Subduction");
        {
          prm.declare_entry("X extent", "2208000",
                            Patterns::Double(),
                            "Set a X extent in meters.");  
          prm.declare_entry("Y extent", "1100000",
                            Patterns::Double(),
                            "Set a Y extent (depth) in meters.");                              
          // const double Yres=1100000;
          prm.declare_entry("Left border", "0",
                            Patterns::Double(),
                            "Set a border or not where there is no lithosphere");  
          prm.declare_entry("Right border", "0",
                            Patterns::Double(),
                            "Set a border or not where there is no lithosphere");                             
          //const double border_right=Xres-0;
          prm.declare_entry("Trench position", "600000",
                            Patterns::Double(),
                            "Set the position of the trench in meters");
          prm.declare_entry("Sediments position", "1800000",
                            Patterns::Double(),
                            "Set the position of the foreland sediments in meters");                            
          prm.declare_entry("Craton position", "1800000",
                            Patterns::Double(),
                            "Set the position of the Craton in meters"); 
          prm.declare_entry("Forearc extension", "600000",
                            Patterns::Double(),
                            "Set the extension of the forearc from the trench in meters");                            
          //const double xforearc=xtrench+600000;
          prm.declare_entry("Forearc additional thickness", "20000",
                            Patterns::Double(),
                            "Set the additional thickness of the lithospheric mantle of the forearc in meters, this is known as lithospheric indentor");          
          prm.declare_entry("Sediments thickness", "8000",
                            Patterns::Double(),
                            "Set thickness of the foreland sediments in meters"); 
          prm.declare_entry("Sediments extension", "300000",
                            Patterns::Double(),
                            "Set extension of the foreland sediments in meters");  
          prm.declare_entry("Weak layer thickness", "5000",
                            Patterns::Double(),
                            "Set thickness of the weak layer in meters"); 
          prm.declare_entry("Orogenic Upper crust", "33000",
                            Patterns::Double(),
                            "Set thickness of the orogenic upper crust in meters");
          prm.declare_entry("Foreland Upper crust thickness", "26000",
                            Patterns::Double(),
                            "Set thickness of the foreland domain upper crust in meters");
          prm.declare_entry("Orogenic Lower crust thickness", "12000",
                            Patterns::Double(),
                            "Set thickness of the foreland domain lower crust in meters");                            
          prm.declare_entry("Foreland Lower crust thickness", "12000",
                            Patterns::Double(),
                            "Set thickness of the foreland domain lower crust in meters");
          prm.declare_entry("Orogenic mantle lithosphere thickness", "55000",
                            Patterns::Double(),
                            "Set thickness of the Orogenic mantle lithosphere in meters");
          prm.declare_entry("Cratonic mantle lithosphere thickness", "100000",
                            Patterns::Double(),
                            "Set thickness of the Cratonic mantle lithosphere in meters");
          prm.declare_entry("Oceanic mantle lithosphere thickness", "73000",
                            Patterns::Double(),
                            "Set thickness of the Oceanic mantle lithosphere in meters");
          prm.declare_entry("Oceanic crust thickness", "7000",
                            Patterns::Double(),
                            "Set thickness of the Oceanic crust in meters");
          prm.declare_entry("Sticky air thickness", "0",
                            Patterns::Double(),
                            "Set thickness of the Sticky air in meters");
          prm.declare_entry("Transition zone depth 1", "410000",
                            Patterns::Double(),
                            "Set depth of the first transition zone in meters");      
          prm.declare_entry("Transition zone depth 2", "520000",
                            Patterns::Double(),
                            "Set depth of the first transition zone in meters");
          prm.declare_entry("Transition zone depth 3", "670000",
                            Patterns::Double(),
                            "Set depth of the first transition zone in meters");
          prm.declare_entry("Forearc thickening angle", "5.71",
                            Patterns::Double(),
                            "Set the thickening angle of the orogenic lithosphere of the forearc in degree");
          prm.declare_entry("Craton thickening angle", "35",
                            Patterns::Double(),
                            "Set the thickening angle of the cratonic mantle lithosphere in degree");                                
          prm.declare_entry("Additional interface thickness", "7000",
                            Patterns::Double(),
                            "Set an additional thickness for the weak layer at the interface between the plates in meters");
          prm.declare_entry("Radius of curvature", "650000",
                            Patterns::Double(),
                            "Set the radius of curvature of the subducting oceanic plate in meters, this is similar to the Arc length");
          prm.declare_entry("Angle of curvature", "35",
                            Patterns::Double(),
                            "Set the angle of curvature of subducting oceanic plate in degree");
          prm.declare_entry("Slab tip depth", "50000",
                            Patterns::Double(),
                            "Set the initial depth of the oceanic plate, this value has to be consistent with the radius of curvature as well as the angle.");                                                                                                                                    
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();      
    }


    template <int dim>
    void Subduction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {      
        prm.enter_subsection("Subduction");
        {
           Xres = prm.get_double("X extent");
           hdepth = prm.get_double("Y extent");
           border_left = prm.get_double("Left border");
           border_right = Xres - prm.get_double("Left border");
           xtrench = prm.get_double("Trench position");
           xsedim = prm.get_double("Sediments position");
           xcraton = prm.get_double("Craton position");
           xforearc = xtrench +  prm.get_double("Forearc extension");
           dforearc = prm.get_double("Forearc additional thickness");
           dsedim = prm.get_double("Sediments thickness");
           lsedim = prm.get_double("Sediments extension");
           dweak = prm.get_double("Weak layer thickness");
           dcrust1 = prm.get_double("Orogenic Upper crust");
           dcrust2 = prm.get_double("Foreland Upper crust thickness");
           dlwcrust1 = prm.get_double("Orogenic Lower crust thickness");
           dlwcrust2 = prm.get_double("Foreland Lower crust thickness");
           dmantcont = prm.get_double("Orogenic mantle lithosphere thickness");
           dcraton = prm.get_double("Cratonic mantle lithosphere thickness");
           dmantoce = prm.get_double("Oceanic mantle lithosphere thickness");
           dgabbro = prm.get_double("Oceanic crust thickness");
           dsticky = prm.get_double("Sticky air thickness");
           TZ1 = prm.get_double("Transition zone depth 1");
           TZ2 = prm.get_double("Transition zone depth 2");
           TZ3 = prm.get_double("Transition zone depth 3");
          //  dlithoforearc = dcrust1+dlwcrust1+dmantcont+dsticky+dforearc;
           angleforearc = prm.get_double("Forearc thickening angle");
           angle4 = prm.get_double("Craton thickening angle");
           weakbottom = prm.get_double("Additional interface thickness");
           segment = prm.get_double("Radius of curvature");
           angle = prm.get_double("Angle of curvature");
           tip = prm.get_double("Slab tip depth");

           dlithocont = dcrust1+dlwcrust1+dmantcont+dsticky;
           dlithocraton = dsedim+dcrust2+dlwcrust2+dcraton+dsticky;
           hcrust3=dsedim+dcrust2;
           hlwcrust3=dsedim+dcrust2+dlwcrust2;

          //Find the radius by rearranging Arclength = (2* numbers::PI*r)/(360/angle)
           radius = (segment*(360/angle))/(2* numbers::PI)-dsticky;
           radiusweak=radius -dweak;
           radiusgabbro = radius - dgabbro;
           radiusmantoce = radius - dmantoce;
           rpyweak=hdepth-dsticky - radius;
           rpygabbro = hdepth-dsticky - radius-dweak;
           rpymantoce = hdepth-dsticky - radius-dweak-dgabbro;

          //set it so the bounds of the circle are at the surface and the
          //center is aligned with the trench location.

           rpx = xtrench;

          //Determine the y location of the end of the segment
           segment_y = std::sin(90-angle)*(radiusmantoce)+rpymantoce;
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();         
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Subduction,
                                              "subduction",
                                              "Build you subduction model in 2D")
  }
}


