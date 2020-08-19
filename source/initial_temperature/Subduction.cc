// /*
//   Copyright (C) 2011 - 2020 by the authors of the ASPECT code.
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
//   along with ASPECT; see the file doc/COPYING.  If not see
//   <http://www.gnu.org/licenses/>.
// */
// #include <aspect/postprocess/interface.h>
// #include <aspect/geometry_model/box.h>
// #include <aspect/boundary_temperature/box.h>
// #include <aspect/utilities.h>
// //#include <aspect/global.h>
// #include <aspect/geometry_model/interface.h>
// #include <aspect/initial_temperature/interface.h>
// #include <aspect/simulator_access.h>
// #include <deal.II/numerics/vector_tools.h>

// #include <deal.II/base/parsed_function.h>

// #include <iostream>

// #include <algorithm>
// #include <iomanip>
// #include <vector>


// namespace aspect
// {
//   namespace InitialTemperature
//   {

//         using namespace dealii;

//         /**
//         * A class that implements temperature initial conditions based on a
//         * functional description provided in the input file.
//         *
//         * @ingroup InitialTemperatures
//         */
//         template <int dim>
//         class Subduction : public Interface<dim>, public SimulatorAccess<dim>
//         {
//         public:

//             virtual void initialize ();
            
//             /**
//             * Constructor.
//             */
//             //Function ();

//             /**
//             * Return the initial temperature as a function of position.
//             */
//             double initial_temperature (const Point<dim> &position) const override;

//             /**
//             * Declare the parameters this class takes through input files. The
//             * default implementation of this function does not describe any
//             * parameters. Consequently, derived classes do not have to overload
//             * this function if they do not take any runtime parameters.
//             */
//             static
//             void
//             declare_parameters (ParameterHandler &prm);

//             /**
//             * Read the parameters this class declares from the parameter file.
//             * The default implementation of this function does not read any
//             * parameters. Consequently, derived classes do not have to overload
//             * this function if they do not take any runtime parameters.
//             */
//             void
//             parse_parameters (ParameterHandler &prm) override;

//             // double interpolate(std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate );
//             double interpolate(const double x,bool interpa,bool interpb, bool interpc,bool interpd,bool interpe);
//             // const std::vector<double> xData;
//             // const std::vector<double> yData;
//             // const char outer3;
//             // const char outer4;
//             // const char outer5;            
//             bool extrapolate;
//             bool interpa;
//             bool interpb;
//             bool interpc;
//             bool interpd;
//             bool interpe;
//             const std::vector<double> xData;
//             const std::vector<double> yData;
//         private:

//             // char outer3;
//             // char outer4;
//             // char outer5;    

//             double Ts;
//             double tlitho1;
//             double tlitho2;
//             double tlitho;
//             //double Tp;

//             double hdepth;
//             double Xres;
//             // double hdepth=1100000;
//             double border_left;
//             double border_right;
//             double xtrench;
//             double xsedim;
//             double xcraton;
//             double xforearc;
//             double dforearc;// how thicker is the continental lithospheric mantle at the forearc
//             double dsedim;
//             double lsedim;
//             double dweak;
//             double dcrust1;
//             double dcrust2;
//             double dlwcrust1;
//             double dlwcrust2;
//             double dmantcont;
//             double dcraton;
//             double dmantoce;
//             double dgabbro;
//             double dsticky;
//             double TZ1;
//             double TZ2; 
//             double TZ3; 

//             // double dlithoforearc = dcrust1+dlwcrust1+dmantcont+dsticky+dforearc;
//             double dlithocont;
//             double dlithocraton;
//             double hcrust3;
//             double hlwcrust3;

//             double angleforearc;
//             double angle4;
//             // double Ts=293;
//             // double tlitho1=1601;
//             // double tlitho2=1613;
//             // double tlitho=1573;
//             // double Tp=1610;

//             double weakbottom;
//             double segment;  
//             double angle;


//             //Find the radius by rearranging Arclength = (2* numbers::PI*r)/(360/angle)
//             double radius;
//             double radiusweak;
//             double radiusgabbro;
//             double radiusmantoce;
//             double rpyweak;
//             double rpygabbro;
//             double rpymantoce;

//             // double htopslab;
//             // double hbottomslab;  
//             // double hslabtip;

//             const std::vector<double> htopslab;
//             const std::vector<double> hbottomslab;  
//             const std::vector<double> hslabtip;  
//             const std::vector<double> vec_X;          

//             //set it so the bounds of the circle are at the surface and the
//             //center is aligned with the trench location.

//             double rpx = xtrench;

//             //Determine the y location of the end of the segment
//             double segment_y;
//             double tip;
//             // dealii::Table<2,double>Slab;
//             // std::vector<double> t;
//             const std::vector<double> circle_top_slab_x;
//             const std::vector<double> circle_top_slab_y;
//             std::vector<double> interp_top_slab_y;
//             const std::vector<double> circle_bottom_slab_x;
//             const std::vector<double> circle_bottom_slab_y;
//             std::vector<double> interp_bottom_slab_y;   
//             // std::vector<double> vec;      
//         };


//         // // Returns interpolated value at x from parallel arrays ( xData, yData )
//         // //   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//         // //   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
//         // template <int dim>
//         // double
//         // interpolate( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate )
//         // {
//         // int size = xData.size();

//         // int i = 0;                                                                  // find left end of interval for interpolation
//         // if ( x >= xData[size - 2] )                                                 // special case: beyond right end
//         // {
//         //     i = size - 2;
//         // }
//         // else
//         // {
//         //     while ( x > xData[i+1] ) i++;
//         // }
//         // double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
//         // if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
//         // {
//         //     if ( x < xL ) yR = yL;
//         //     if ( x > xR ) yL = yR;
//         // }

//         // double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

//         // return yL + dydx * ( x - xL );                                              // linear interpolation
//         // }
//         template <int dim>        
//         double 
//         interpolate(const double x,bool interpa,bool interpb, bool interpc,bool interpd,bool interpe);
//         {
//             bool extrapolate = false;
//             if (interpa){
//                 int size = circle_top_slab_x.size();
//                 xData=circle_top_slab_x; 
//                 yData=circle_top_slab_y;
//             }
//             else if(interpb){
//                 int size = circle_bottom_slab_x.size();
//                 xData=circle_bottom_slab_x; 
//                 yData=circle_bottom_slab_y;
//             }
//             else if(interpc){
//                 int size = vec_X.size();
//                 xData=vec_X; 
//                 yData=htopslab;
//                 x=position_x;
//             }
//             else if(interpd){
//                 int size = vec_X.size();
//                 xData=vec_X;
//                 yData=hbottomslab;
//                 x=position_x;               
//             }
//             else if(interpe){
//                 int size = vec_X.size();
//                 xData=vec_X;
//                 yData=hslabtip;
//                 x=position_x; 
//             }
//             // if (input=="1"){
//             //     int size = circle_top_slab_x.size();
//             //     xData=circle_top_slab_x; 
//             //     yData=circle_top_slab_y;
//             // }
//             // else if(input=="2"){
//             //     int size = circle_bottom_slab_x.size();
//             //     xData=circle_bottom_slab_x; 
//             //     yData=circle_bottom_slab_y;
//             // }
//             // else if(input=="3"){
//             //     int size = vec_X.size();
//             //     xData=vec_X; 
//             //     yData=htopslab;
//             //     x=position_x;
//             // }
//             // else if(input=="4"){
//             //     int size = vec_X.size();
//             //     xData=vec_X;
//             //     yData=hbottomslab;
//             //     x=position_x;               
//             // }
//             // else if(input=="5"){
//             //     int size = vec_X.size();
//             //     xData=vec_X;
//             //     yData=hslabtip;
//             //     x=position_x; 
//             // }

//         int i = 0;                                                                  // find left end of interval for interpolation
//         if ( x >= xData[size - 2] )                                                 // special case: beyond right end
//         {
//             i = size - 2;
//         }
//         else
//         {
//             while ( x > xData[i+1] ) i++;
//         }
//         double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
//         if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
//         {
//             if ( x < xL ) yR = yL;
//             if ( x > xR ) yL = yR;
//         }

//         double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

//         return yL + dydx * ( x - xL );           
//         }

//     template <int dim>
//     void
//     Subduction<dim>::initialize()
//     {

    
//         // std::cout << "t: ";
//         // for (auto iv: t) {
//         //     std::cout << iv << " ";
//         // }
//         // std::cout << "\n";


//         // Option 2 : use circle functions
//         const int adeg=90; const int bdeg=0; //degree for the orientation of the section of the trigonometric circle
//         const double rpytop=hdepth-dsticky-radius+weakbottom;

//         // std::vector<double> t(100);

//         std::vector<int> t;

//         // set some values:
//         for (int i=0; i<adeg; ++i){
//             double vec =90-i;
//             t.push_back(vec);
//         }    
   
//         std::vector<double>circle_top_slab_x(adeg);
//         std::vector<double>circle_top_slab_y(adeg);
//         std::vector<double>circle_bottom_slab_x(adeg);
//         std::vector<double>circle_bottom_slab_y(adeg);

//         // Set up some points for interpolation in xVals
//         const int NPTS = Xres;
//         std::vector<double> xVals, interp_bottom_slab_y;
//         for ( int i = xtrench; i < NPTS; i++ ) xVals.push_back( (double)i );
        

//         for (unsigned int i=0; i < t.size(); ++i)
//         {
//         circle_top_slab_x[i] = radius*std::cos(t[i]) + rpx;
//         circle_top_slab_y[i] = radius*std::sin(t[i]) + rpytop;
//         circle_bottom_slab_x[i] = radiusmantoce*std::cos(t[i]) + rpx;
//         circle_bottom_slab_y[i] = radiusmantoce*std::sin(t[i]) + rpymantoce;        
//         }

//         // Interpolate
//         for ( double x : xVals )
//         {
//             // double y_top = interpolate( circle_top_slab_x, circle_top_slab_y, x, false );
//             // double y_bottom = interpolate( circle_bottom_slab_x, circle_bottom_slab_y, x, false );
//             double y_top = interpolate(x, true,false,false,false,false);
//             double y_bottom = interpolate(x, false,true,false,false,false);            
//             interp_bottom_slab_y.push_back( y_bottom );
//             interp_top_slab_y.push_back( y_top );
//         }

//         std::vector<double>htopslab(Xres);
//         std::vector<double>hbottomslab(Xres);
//         std::vector<double>hslabtip(Xres);

//         for (unsigned int j=0; j < Xres; ++j){
//             if (j<= xtrench){
//                 htopslab[j] = hdepth-dsticky;
//                 hbottomslab[j] = hdepth-dsticky-dgabbro-dweak-dmantoce;
//             }
//             else if (j> xtrench){
//                 if(interp_bottom_slab_y[j-xtrench]>hdepth-dsticky-dgabbro-dmantoce-dweak-tip){
//                     interp_bottom_slab_y[j-xtrench]=interp_bottom_slab_y[j-xtrench];
//                 }
//                 else{
//                     interp_bottom_slab_y[j-xtrench]=hdepth-dsticky-dgabbro-dmantoce-dweak-tip;
//                 }
//                 if (interp_top_slab_y[j-xtrench]>hdepth){
//                     htopslab[j]= hdepth-dsticky;
//                     hbottomslab[j]= interp_bottom_slab_y[j-xtrench];
//                 }
//                 else{
//                     htopslab[j]= interp_top_slab_y[j-xtrench];
//                     hbottomslab[j]= interp_bottom_slab_y[j-xtrench];
//                 }
//             }
//             if(htopslab[j]>hbottomslab[j]){
//                 hslabtip[j]=htopslab[j]-hbottomslab[j];
//             }
//             else{
//                 hslabtip[j]=1;
//             }            
//         }

//             std::vector<double> vec_X (Xres);

//             // set some values:
//             for (int i=0; i<Xres; ++i){
//                 vec_X.push_back(i);
//             }           
//             // char outer3 = "3";
//             // char outer4 = "4";
//             // char outer5 = "5";
//     }
                          
//         // template <int dim>
//         // double
//         // interpolate_point( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate )
//         // {

//         // }

//         template <int dim>
//         double
//         Subduction<dim>::
//         initial_temperature (const Point<dim> &position) const
//         {            
//             // double Top = interpolate( vec_X, htopslab, position[0], false );
//             // double Bottom = interpolate( vec_X, hbottomslab, position[0], false );
//             // double Tip = interpolate( vec_X, hslabtip, position[0], false );
//             const double position_x=position[0];
//             double Top = interpolate( position_x,false,false,true,false,false);
//             double Bottom = interpolate( position_x,false,false,false,true,false);
//             double Tip = interpolate(position_x,false,false,false,false,true);

//             if(position[1]<=Top && position[1]>=Bottom){
//                 return Ts+(tlitho-Ts)*(Top-position[1])/(Tip); 
//             }       
//             else if(sqrt(pow((position[0]-rpx),2) + pow((position[1]-rpyweak),2)) > radius+weakbottom && position[0]> xtrench && position[1]>=(hdepth-dsticky-dcrust1-dlwcrust1-dmantcont+std::tan(angleforearc*(numbers::PI/180))*(-xforearc+position[0])) && position[0]<=xforearc)
//             {
//                 return Ts +(tlitho1-Ts)*(hdepth-dsticky- position[1])/((dcrust1+dlwcrust1+dmantcont)-std::tan(angleforearc*(numbers::PI/180))*(-xforearc+position[0]));
//             }
//             else if(position[0]>=xforearc && position[0]<xcraton && position[1]>=hdepth-dcrust1-dlwcrust1-dmantcont-dsticky){
//                 return Ts+(tlitho1-Ts)*(hdepth-dsticky- position[1])/(dcrust1+dlwcrust1+dmantcont);
//             }  
//             else if(position[0]>=xcraton && position[0]<=xcraton+((dlithocraton-dlithocont)/std::tan(angle4*numbers::PI/180)) && position[1]>=(hdepth-dsticky-dcrust1-dlwcrust1-dmantcont-std::tan(angle4*(numbers::PI/180))*(-xcraton+position[0]))){
//                 return Ts +(tlitho2-Ts)*(hdepth-dsticky- position[1])/((dcrust1+dlwcrust1+dmantcont)+std::tan(angle4*(numbers::PI/180))*(-xcraton+position[0]));
//             }
//             else if(position[0]>xcraton+((dlithocraton-dlithocont)/std::tan(angle4*numbers::PI/180)) && position[1]>=hdepth-dsedim-dcrust2-dlwcrust2-dcraton-dsticky){
//                 return Ts +(tlitho2-Ts)*(hdepth-dsticky- position[1])/(dsedim+dcrust2+dlwcrust2+dcraton);
//             }
//             else{
//                 return 2400; // this temperature is a threshold that is replaced later by an adiabtic profil.
//             }
//         }

//         template <int dim>
//         void
//         Subduction<dim>::declare_parameters (ParameterHandler &prm)
//         {
//         prm.enter_subsection ("Initial temperature model");
//         {
//             prm.enter_subsection("Subduction");
//             {
//             prm.declare_entry("Surface temperature", "293",
//                                 Patterns::Double(),
//                                 "Surface temperature in Kelvin, use 293 instead of 273K has for advantage to avoid any zero temperature over time and to not refine the cold Sticky Air if there is based on the isotherms (temperature)");  
//             prm.declare_entry("Orogenic LAB temperature", "1601",
//                                 Patterns::Double(),
//                                 "Orogenic Lithosphere-Astehnosphere boundary temperature");
//             prm.declare_entry("Cratonic LAB temperature", "1613",
//                                 Patterns::Double(),
//                                 "Cratonic Lithosphere-Astehnosphere boundary temperature");                                         
//             }
//             prm.leave_subsection();
//         }
//         prm.leave_subsection();        
//         }

//         template <int dim>
//         void
//         Subduction<dim>::parse_parameters (ParameterHandler &prm)
//         {
//         prm.enter_subsection ("Initial temperature model");
//         {
//             prm.enter_subsection("Subduction");
//             {
//                 Ts = prm.get_double("Surface temperature");
//                 tlitho1 = prm.get_double("Orogenic LAB temperature");
//                 tlitho2 = prm.get_double("Cratonic LAB temperature");
//                 //Tp = prm.get_double("LAB temperature");
//             }
//             prm.leave_subsection();
//         }
//         prm.leave_subsection();

//                 prm.enter_subsection("Initial composition model");
//         {      
//             prm.enter_subsection("Subduction");
//             {
//             Xres = prm.get_double("X extent");
//             hdepth = prm.get_double("Y extent");
//             border_left = prm.get_double("Left border");
//             border_right = Xres - prm.get_double("Left border");
//             xtrench = prm.get_double("Trench position");
//             xsedim = prm.get_double("Sediments position");
//             xcraton = prm.get_double("Craton position");
//             xforearc = xtrench +  prm.get_double("Forearc extension");
//             dforearc = prm.get_double("Forearc additional thickness");
//             dsedim = prm.get_double("Sediments thickness");
//             lsedim = prm.get_double("Sediments extension");
//             dweak = prm.get_double("Weak layer thickness");
//             dcrust1 = prm.get_double("Orogenic Upper crust");
//             dcrust2 = prm.get_double("Foreland Upper crust thickness");
//             dlwcrust1 = prm.get_double("Orogenic Lower crust thickness");
//             dlwcrust2 = prm.get_double("Foreland Lower crust thickness");
//             dmantcont = prm.get_double("Orogenic mantle lithosphere thickness");
//             dcraton = prm.get_double("Cratonic mantle lithosphere thickness");
//             dmantoce = prm.get_double("Oceanic mantle lithosphere thickness");
//             dgabbro = prm.get_double("Oceanic crust thickness");
//             dsticky = prm.get_double("Sticky air thickness");
//             TZ1 = prm.get_double("Transition zone depth 1");
//             TZ2 = prm.get_double("Transition zone depth 2");
//             TZ3 = prm.get_double("Transition zone depth 3");
//             //  dlithoforearc = dcrust1+dlwcrust1+dmantcont+dsticky+dforearc;
//             angleforearc = prm.get_double("Forearc thickening angle");
//             angle4 = prm.get_double("Craton thickening angle");
//             weakbottom = prm.get_double("Additional interface thickness");
//             segment = prm.get_double("Radius of curvature");
//             angle = prm.get_double("Angle of curvature");
//             tip = prm.get_double("Slab tip depth");
            
//             dlithocont = dcrust1+dlwcrust1+dmantcont+dsticky;
//             dlithocraton = dsedim+dcrust2+dlwcrust2+dcraton+dsticky;
//             hcrust3=dsedim+dcrust2;
//             hlwcrust3=dsedim+dcrust2+dlwcrust2;

//             dsticky=0;
//             //Find the radius by rearranging Arclength = (2* numbers::PI*r)/(360/angle)
//             radius = (segment*(360/angle))/(2* numbers::PI)-dsticky;
//             radiusweak=radius - dweak;
//             radiusgabbro = radius - dgabbro;
//             radiusmantoce = radius - dmantoce;
//             rpyweak=hdepth-dsticky - radius;
//             rpygabbro = hdepth-dsticky - radius-dweak;
//             rpymantoce = hdepth-dsticky - radius-dweak-dgabbro;

//             //set it so the bounds of the circle are at the surface and the
//             //center is aligned with the trench location.

//             rpx = xtrench;

//             //Determine the y location of the end of the segment
//             segment_y = std::sin(90-angle)*(radiusmantoce)+rpymantoce;
            
//             }
//             prm.leave_subsection();
//         }
//       }
//     }
// }

// // explicit instantiations
// namespace aspect
// {
//   namespace InitialTemperature
//   {
//     ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Subduction,
//                                               "subduction",
//                                               "subduction temperature file")
//   }
// }


    

