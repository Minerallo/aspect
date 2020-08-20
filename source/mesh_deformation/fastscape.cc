/*
  Copyright (C) 2019 by the authors of the ASPECT code.
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

#include <aspect/mesh_deformation/fastscape.h>
#include <aspect/geometry_model/box.h>
#include <deal.II/numerics/vector_tools.h>
#include <ctime>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshDeformation
  {

    template <int dim>
    void
    FastScape<dim>::initialize()
    {
      const GeometryModel::Box<dim> *geometry = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());

      AssertThrow(geometry != nullptr,
                  ExcMessage("Fastscape can only be run with a box model"));

      //Initialize parameters for restarting fastscape
      restart = this->get_parameters().resume_computation;
      restart_step = 0;

      //second is for maximum coordiantes, first for minimum.
      if (use_boxlitho_2d)
      {
        grid_extent[0].first = 0;
        grid_extent[1].first = 0;
        grid_extent[0].second = x_extent_2d;
        grid_extent[1].second = y_extent_2d;
        x_repetitions = x_repetitions_2d;
        y_repetitions = y_repetitions_2d;
      }
      else
      {
        for (unsigned int i = 0; i < dim; ++i)
        {
          grid_extent[i].first = 0;
          grid_extent[i].second = geometry->get_extents()[i];
        }
        //TODO: There has to be a better type to use to get this.
//         const std::pair<int, int> repetitions = geometry->get_repetitions();
//         x_repetitions = repetitions.first;
//         y_repetitions = repetitions.second;
      }

      //Set nx and dx, as these will be the same regardless of dimension.
      nx = 3 + std::pow(2, surface_resolution + additional_refinement) * x_repetitions;
      dx = grid_extent[0].second / (nx - 3);
      x_extent = grid_extent[0].second + 2 * dx; //geometry->get_extents()[0]

      //sub intervals are 1 less than points.
      table_intervals[0] = nx - 3;
      //TODO: it'd be best to not have to use dim-1 intervals at all.
      table_intervals[dim - 1] = 1;

      if (dim == 2)
      {
        dy = dx;
        y_extent = round(grid_extent[1].second / dy) * dy + 2 * dy; //x_extent; //grid_extent[0].second*3+2*dy;
        ny = 1 + y_extent / dy;
      }

      if (dim == 3)
      {
        ny = 3 + std::pow(2, surface_resolution + additional_refinement) * y_repetitions;
        dy = grid_extent[1].second / (ny - 3);
        table_intervals[1] = ny - 3;
        y_extent = geometry->get_extents()[1] + 2 * dy;
      }

      //Determine array size to send to fastscape
      array_size = nx * ny - 1;
      std::cout << "nx: " << nx << "ny: " << ny << "y_extent: " << y_extent << "x_extent: " << x_extent << "x repetition: "<<x_repetitions << " dy " << dy << "  " << std::endl;
    }


    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                             ConstraintMatrix &mesh_velocity_constraints,
                                                             const std::set<types::boundary_id> &boundary_ids) const
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Fastscape plugin");
      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const bool top_boundary = boundary_ids.find(relevant_boundary) == boundary_ids.begin();
      int current_timestep = this->get_timestep_number();

      double a_dt = this->get_timestep();
      if (this->convert_output_to_years())
      {
        a_dt = this->get_timestep() / year_in_seconds;
      }

      //We only want to run fastscape if there was a change in time, and if we're at the top boundary.
      if (a_dt > 0 && top_boundary)
      {
        /*
           * Initialize a vector of temporary variables to hold: z component, index, Vx, Vy, and Vz.
           */
        std::vector<std::vector<double>> temporary_variables(dim + 3, std::vector<double>());
        std::vector<double> V(array_size);
        double precision = 0.4;

        // Get a quadrature rule that exists only on the corners, and increase the refinement if specified.
        const QIterated<dim - 1> face_corners(QTrapez<1>(),
                                              pow(2, additional_refinement + resolution_difference));

        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                         this->get_fe(),
                                         face_corners,
                                         update_values |
                                             update_quadrature_points);

        typename DoFHandler<dim>::active_cell_iterator
            cell = this->get_dof_handler().begin_active(),
            endc = this->get_dof_handler().end();

        for (; cell != endc; ++cell)
          if (cell->is_locally_owned() && cell->at_boundary())
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
              if (cell->face(face_no)->at_boundary())
              {
                if (cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                std::vector<Tensor<1, dim>> vel(face_corners.size());
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                {
                  const Point<dim> vertex = fe_face_values.quadrature_point(corner);
                  //This tells us which x point we're at
                  double indx = 2 + vertex(0) / dx;

                  // If our x or y index isn't close to a whole number, then it's likely an artifact
                  // from using an over-resolved quadrature rule, in that case ignore it.
                  //std::cout<<indx<<"  "<<round(indx)<<"  "<<abs(indx - round(indx))<<std::endl;
                  if (abs(indx - round(indx)) >= precision)
                    continue;

                  /*
                         * If we're in 2D, we want to take the values and apply them to every row of X points.
                         * TODO: I need to test that 2D still works after changing how we store these variables.
                         */
                  if (dim == 2)
                  {
                    for (int ys = 0; ys < ny; ys++)
                    {
                      /*
                                 * Fastscape indexes from 1 to n, starting at X and Y = 0, and increases
                                 * across the X row. At the end of the row, it jumps back to X = 0
                                 * and up to the next X row in increasing Y direction. We track
                                 * this to correctly place the variables later on.
                                 */
                      double index = round(indx) + nx * ys;

                      temporary_variables[0].push_back(vertex(dim - 1));
                      temporary_variables[1].push_back(index - 1);

                      //std::cout<<index-1<<"  "<<vertex(dim-1)<<"  ";

                      for (unsigned int i = 0; i < dim; ++i)
                        temporary_variables[i + 2].push_back(vel[corner][i] * year_in_seconds);

                      if (diffuse_area && use_intervals)
                      {
                        if (vertex(0) > x_min_diffusion && vertex(0) < x_max_diffusion)
                        {
                          temporary_variables[4].push_back(1);
                        }
                        else
                        {
                          temporary_variables[4].push_back(0);
                        }
                      }
                      else{
                        temporary_variables[4].push_back(0);
                      }
                      // else if (diffuse_area && minimum_topography_diffusion){
                      //     //double min_h = 99999999;
                      //     mutable double idx_h = 0;    
                      //   if (current_timestep == 0)                    
                      //     temporary_variables[4].push_back(0);
                      //   else if (current_timestep != 0){
                      //     if (vertex(0) > ((idx_h)-2 * dx) - 50000 && vertex(0) < ((idx_h)-2 * dx) + 50000)
                      //     {
                      //       temporary_variables[4].push_back(1);
                      //       // std::cout << "Yes" << std::endl;
                      //     }
                      //     else{
                      //       temporary_variables[4].push_back(0);
                      //       // std::cout << "No" << std::endl;
                      //     }
                      //   }
                      //   // std::cout << "New temporary variable " << std::endl;
                      // }
                    }
                    //std::cout<<std::endl;
                  }

                  if (dim == 3)
                  {
                    double indy = 2 + vertex(1) / dy;
                    if (indy - floor(indy) >= precision)
                      continue;

                    double index = (indy - 1) * nx + indx;

                    temporary_variables[1].push_back(index - 1);
                    temporary_variables[0].push_back(vertex(dim - 1)); //z component

                    for (unsigned int i = 0; i < dim; ++i)
                      temporary_variables[i + 2].push_back(vel[corner][i] * year_in_seconds);
                  }
                }
              }

        //Run fastscape on single processor.
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Fastscape 1 proc");
          /*
               * Initialize the variables that will be sent to fastscape.
               * These have to be doubles of array_size, which C++ doesn't like,
               * so they're initialized this way.
               */
          std::unique_ptr<double[]> h(new double[array_size]);
          std::unique_ptr<double[]> vx(new double[array_size]);
          std::unique_ptr<double[]> vy(new double[array_size]);
          std::unique_ptr<double[]> vz(new double[array_size]);
          std::unique_ptr<double[]> kf(new double[array_size]);
          std::unique_ptr<double[]> kd(new double[array_size]);
          std::unique_ptr<double[]> slopep(new double[array_size]);
          std::unique_ptr<double[]> test(new double[array_size]);
          int istep = 0;
          int steps = nstep;
          std::srand(fs_seed);
          std::vector<double> h_old(array_size);

          //Create variables for output directory and restart file
          std::string filename;
          filename = this->get_output_directory();
          const char *c = filename.c_str();
          int length = filename.length();
          const std::string restart_filename = filename + "fastscape_h_restart.txt";
          const std::string restart_step_filename = filename + "fastscape_steps_restart.txt";

          //Initialize all FastScape variables.
          for (int i = 0; i <= array_size; i++)
          {
            h[i] = 0;
            vx[i] = 0;
            vz[i] = 0;
            vy[i] = 0;
            slopep[i] = 0;
            kf[i] = kff;
            kd[i] = kdd;
          }

          //Get info from first processor.
          //std::cout<<"nx: "<<nx<<"ny: "<<ny<<"array size: "<<array_size<<"  "<<temporary_variables[1].size()<<std::endl;
          for (unsigned int i = 0; i < temporary_variables[1].size(); i++)
          {
            h[temporary_variables[1][i]] = temporary_variables[0][i];
            vx[temporary_variables[1][i]] = temporary_variables[2][i];
            vz[temporary_variables[1][i]] = temporary_variables[dim + 1][i];
            //std::cout<<temporary_variables[1][i]<<"  "<<temporary_variables[0][i]<<std::endl;
            test[temporary_variables[1][i]] = temporary_variables[4][i];

            if (dim == 2)
              vy[temporary_variables[1][i]] = 0;

            if (dim == 3)
              vy[temporary_variables[1][i]] = temporary_variables[3][i];
          }

              for (unsigned int p = 1; p < Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
          {

            //First, find out the size of the array a processor wants to send.
            MPI_Status status;
            MPI_Probe(p, 42, this->get_mpi_communicator(), &status);
            int incoming_size = 0;
            MPI_Get_count(&status, MPI_DOUBLE, &incoming_size);

            //Resize the array so it fits whatever the processor sends.
            for (unsigned int i = 0; i < temporary_variables.size(); ++i)
            {
              temporary_variables[i].resize(incoming_size);
            }

            for (unsigned int i = 0; i < temporary_variables.size(); i++)
              MPI_Recv(&temporary_variables[i][0], incoming_size, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

            //Now, place the numbers into the correct place based off the index.
            for (unsigned int i = 0; i < temporary_variables[1].size(); i++)
            {
              h[temporary_variables[1][i]] = temporary_variables[0][i];
              vx[temporary_variables[1][i]] = temporary_variables[2][i];
              vz[temporary_variables[1][i]] = temporary_variables[dim + 1][i];
              test[temporary_variables[1][i]] = temporary_variables[4][i];
              // std::cout << "Yes2" << std::endl;

              if (dim == 2)
                vy[temporary_variables[1][i]] = 0;
              if (dim == 3)
                vy[temporary_variables[1][i]] = temporary_variables[3][i];
            }
          }

            //Initializaion from the beginning
            if (diffuse_area && new_diffusion_step == 0)
            {
                for (int i = 0; i <= array_size; i++)
                {
                  if (test[i] == 1)
                    kd[i] = new_diffusion;
                  else
                    kd[i] = kdd;
                }
            }
            else if (diffuse_area && new_diffusion_step != 0)
            {
              if (current_timestep >= new_diffusion_step){
                for (int i = 0; i <= array_size; i++)
                {
                  if (test[i] == 1)
                    kd[i] = new_diffusion;
                  else
                    kd[i] = kdd;
                }
              }
            }

              if (current_timestep == 1 || restart)
              {

                this->get_pcout() << "   Initializing FastScape... " << (1 + surface_resolution + additional_refinement) << " levels, cell size: " << dx << " m." << std::endl;

                //If we are restarting from a checkpoint, load h values for fastscape so we don't lose resolution.
                if (restart)
                {
                  this->get_pcout() << "      Loading FastScape restart file... " << std::endl;
                  restart = false;

                  //Load in h values.
                  std::ifstream in;
                  in.open(restart_filename.c_str());
                  if (in)
                  {
                    int line = 0;

                    while (line <= array_size)
                    {
                      in >> h[line];
                      line++;
                    }

                    in.close();
                  }
                  else if (!in)
                    AssertThrow(false, ExcMessage("Cannot open file to restart FastScape."));

                  /*
                       * Now load the fastscape istep at time of restart.
                       * Reinitializing fastscape always resets this to 0, so here
                       * we keep it in a separate variable to keep track for visualization files.
                       */
                  std::ifstream in_step;
                  in_step.open(restart_step_filename.c_str());
                  if (in)
                  {

                    in_step >> restart_step;
                    in_step.close();
                  }
                  else if (!in_step)
                    AssertThrow(false, ExcMessage("Cannot open file to restart fastscape."));
                }

                //Initialize fastscape with grid and extent.
                fastscape_init_();
                fastscape_set_nx_ny_(&nx, &ny);
                fastscape_setup_();
                fastscape_set_xl_yl_(&x_extent, &y_extent);

                //set boundary conditions
                fastscape_set_bc_(&bc);

                //Initialize topography
                fastscape_init_h_(h.get());

                //Set erosional parameters. May have to move this if sed values are updated over time.
                fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);

                if (use_marine)
                  fastscape_set_marine_parameters_(&sl, &p1, &p2, &z1, &z2, &r, &l, &kds1, &kds2);

                //if (use_strat)
                //  folder_output_(&length, &restart_step, c);
              }
              else
              {
                //If it isn't the first timestep we just want to know current h values in fastscape,
                //and if we aren't advecting the surface then we want to continually send aspect heights.
                if (use_velocity)
                  fastscape_copy_h_(h.get());
              }

          /*
               * Copy the slopes at each point, this will be used to set an H
               * at the ghost nodes if a boundary mass flux is given.
               */
          fastscape_copy_slope_(slopep.get());

          //Now we set the ghost nodes at the left and right boundaries.
          for (int j = 0; j < ny; j++)
          {
            double index_left = nx * j + 1;
            double index_right = nx * (j + 1);
            double slope = 0;

            //Generally, they will always be set to the same values as the
            //inner nodes next to them.
            vz[index_right - 1] = vz[index_right - 2];
            vz[index_left - 1] = vz[index_left];

            vy[index_right - 1] = vy[index_right - 2];
            vy[index_left - 1] = vy[index_left];

            vx[index_right - 1] = vx[index_right - 2];
            vx[index_left - 1] = vx[index_left];

            if (current_timestep == 1 || left_flux == 0)
            {
              //If its the first timestep add in initial slope. If we have no flux,
              //set the ghost node to the node next to it.
              slope = left_flux / kdd;
              h[index_left - 1] = h[index_left] + slope * 2 * dx;
            }
            else
            {
              //If we have flux through boundary, we need to update the height to keep the correct slope.
              //Because the corner nodes always show a slope of zero, this will update them according to
              //the nodes next to them.
              if (j == 0)
                slope = left_flux / kdd - std::tan(slopep[index_left + nx] * boost::math::double_constants::pi / 180);
              else if (j == (ny - 1))
                slope = left_flux / kdd - std::tan(slopep[index_left - nx] * boost::math::double_constants::pi / 180);
              else
                slope = left_flux / kdd - std::tan(slopep[index_left] * boost::math::double_constants::pi / 180);

              h[index_left - 1] = h[index_left - 1] + slope * 2 * dx;
            }

            if (current_timestep == 1 || right_flux == 0)
            {
              slope = right_flux / kdd;
              h[index_right - 1] = h[index_right - 2] + slope * 2 * dx;
            }
            else
            {
              if (j == 0)
                slope = right_flux / kdd - std::tan(slopep[index_right + nx - 2] * boost::math::double_constants::pi / 180);
              else if (j == (ny - 1))
                slope = right_flux / kdd - std::tan(slopep[index_right - nx - 2] * boost::math::double_constants::pi / 180);
              else
                slope = right_flux / kdd - std::tan(slopep[index_right - 2] * boost::math::double_constants::pi / 180);

              h[index_right - 1] = h[index_right - 1] + slope * 2 * dx;
            }

            //If we set the boundaries as periodic, then reset any values to the
            //nodes on the opposite side.
            if (left == 0 && right == 0)
            {
              double side = index_left;
              int jj = 0;

              //If they aren't going the same direction, don't set the ghost nodes.
              if (vx[index_right - 2] > 0 && vx[index_left] >= 0)
              {
                side = index_right;
                jj = 2;
              }
              else if (vx[index_right - 2] <= 0 && vx[index_left] < 0)
              {
                side = index_left;
                jj = 0;
              }
              else
                continue;

              vz[index_right - 1] = vz[side - jj];
              vz[index_left - 1] = vz[side - jj];

              vy[index_right - 1] = vy[side - jj];
              vy[index_left - 1] = vy[side - jj];

              vx[index_right - 1] = vx[side - jj];
              vx[index_left - 1] = vx[side - jj];

              h[index_right - 1] = h[side - jj];
              h[index_left - 1] = h[side - jj];
            }
          }

          //Now do the same for the top and bottom ghost nodes.
          for (int j = 0; j < nx; j++)
          {
            double index_bot = j + 1;
            double index_top = nx * (ny - 1) + j + 1;
            double slope = 0;

            vz[index_bot - 1] = vz[index_bot + nx - 1];
            vz[index_top - 1] = vz[index_top - nx - 1];

            vy[index_bot - 1] = vy[index_bot + nx - 1];
            vy[index_top - 1] = vy[index_top - nx - 1];

            vx[index_bot - 1] = vx[index_bot + nx - 1];
            vx[index_top - 1] = vx[index_top - nx - 1];

            if (current_timestep == 1 || top_flux == 0)
            {
              slope = top_flux / kdd;
              h[index_top - 1] = h[index_top - nx - 1] + slope * 2 * dx;
            }
            else
            {
              if (j == 0)
                slope = top_flux / kdd - std::tan(slopep[index_top - nx] * boost::math::double_constants::pi / 180);
              else if (j == (nx - 1))
                slope = top_flux / kdd - std::tan(slopep[index_top - nx - 2] * boost::math::double_constants::pi / 180);
              else
                slope = top_flux / kdd - std::tan(slopep[index_top - nx - 1] * boost::math::double_constants::pi / 180);

              h[index_top - 1] = h[index_top - 1] + slope * 2 * dx;
            }

            if (current_timestep == 1 || bottom_flux == 0)
            {
              slope = bottom_flux / kdd;
              h[index_bot - 1] = h[index_bot + nx - 1] + slope * 2 * dx;
            }
            else
            {
              if (j == 0)
                slope = bottom_flux / kdd - std::tan(slopep[index_bot + nx] * boost::math::double_constants::pi / 180);
              else if (j == (nx - 1))
                slope = bottom_flux / kdd - std::tan(slopep[index_bot + nx - 2] * boost::math::double_constants::pi / 180);
              else
                slope = bottom_flux / kdd - std::tan(slopep[index_bot + nx - 1] * boost::math::double_constants::pi / 180);

              h[index_bot - 1] = h[index_bot - 1] + slope * 2 * dx;
            }

            if (bottom == 0 && top == 0)
            {
              double side = index_bot;
              int jj = nx;

              if (vy[index_bot + nx - 1] > 0 && vy[index_top - nx - 1] >= 0)
              {
                side = index_top;
                jj = -nx;
              }
              else if (vy[index_bot + nx - 1] <= 0 && vy[index_top - nx - 1] < 0)
              {
                side = index_bot;
                jj = nx;
              }
              else
                continue;

              vz[index_bot - 1] = vz[side + jj - 1];
              vz[index_top - 1] = vz[side + jj - 1];

              vy[index_bot - 1] = vy[side + jj - 1];
              vy[index_top - 1] = vy[side + jj - 1];

              vx[index_bot - 1] = vx[side + jj - 1];
              vx[index_top - 1] = vx[side + jj - 1];

              h[index_bot - 1] = h[side + jj - 1];
              h[index_top - 1] = h[side + jj - 1];
            }
          }

          /*
               * Keep initial h values so we can calculate velocity later.
               * In the first timestep, h will be given from other processors.
               * In other timesteps, we copy h directly from fastscape.
               */
          for (int i = 0; i <= array_size; i++)
          {
            //Initialize random topography noise first time fastscape is called.
            if (current_timestep == 1 && use_velocity)
            {
              double h_seed = (std::rand() % 2000) / 100;
              h[i] = h[i] + h_seed;
            }
            //std::cout<<i<<"  "<<h[i]<<std::endl;

            h_old[i] = h[i];
          }

          //Get current fastscape timestep.
          fastscape_get_step_(&istep);

          //Write a file to store h & step in case of restart.
          //TODO: there's probably a faster way to write these.
          if ((this->get_parameters().checkpoint_time_secs == 0) &&
              (this->get_parameters().checkpoint_steps > 0) &&
              (current_timestep % this->get_parameters().checkpoint_steps == 0))
          {
            std::ofstream out_h(restart_filename.c_str());
            std::ofstream out_step(restart_step_filename.c_str());

            out_step << (istep + restart_step) << std::endl;

            for (int i = 0; i <= array_size; i++)
              out_h << h[i] << std::endl;
          }

          //Find a fastscape timestep that is below our maximum timestep.
          double f_dt = a_dt / steps;
          while (f_dt > max_timestep)
          {
            steps = steps * 2;
            f_dt = a_dt / steps;
          }

          //Set time step
          fastscape_set_dt_(&f_dt);

          //Set velocity components and h.
          if (use_velocity)
          {
            fastscape_set_u_(vz.get());
            fastscape_set_v_(vx.get(), vy.get());
          }
   
          //Change diffusion parameter in time
          if (diffuse_area && minimum_topography_diffusion)
          {      
            this->get_pcout() << "   Tracking the minimum topography for specific diffusion... " << (1 + surface_resolution + additional_refinement) << " levels, cell size: " << dx << " m." << std::endl;
          // Find the minimum elevation
          // There might be a better way to track the minimum of the elevation using min_element
          double min_h = 99999999;
          double idx_h = 0;
          for (int i = 0; i <= nx; i++)
          {
            if (h[nx * 5 + i] < min_h)
            {
              min_h = h[nx * 5 + i];
              idx_h = i;
              // std::cout << "Elevations: " << h[nx * 5 + i] <<"  "<< idx_h << " / " << nx << std::endl;
            }
          }
          std::cout << " Minimum elevation (inverted) : " << min_h << " index : " << idx_h<< " Position in km : " << ((idx_h-2) * dx) /1000 << std::endl;
          // Since the variable are not carried from one timestep to another I redefine the diffusion grid here
          //one index is equal to one element or about 1dx
   
            if (diffusion_sinusoid){

               if(this->get_time()/ year_in_seconds>=time_to_diffuse)     
                  {    
                  std::cout << " It is time to use sinuosid function for Diffusion " <<std::endl;                   
                  for (int i = 0; i <= array_size-nx; i++){
                    double current_ny = i/nx;
                    double current_cut = nx*trunc(current_ny)+idx_h-round(dx_diffusion/2);                        
                    if (i>=current_cut && i <=nx*trunc(current_ny)+idx_h+round(dx_diffusion/2)){
                      kd[i]=new_diffusion*(1-(cos((2*boost::math::double_constants::pi/dx_diffusion)*(i-current_cut))+1)/2);
                    }
                    else{
                      kd[i]=kdd;
                    }
                  }
                }else{
                  std::cout << " Not the time to use sinuosid function for Diffusion " <<std::endl; 
                  for (int i = 0; i <= array_size-nx; i++){ 
                    kd[i]=kdd;
                  }                 
                }  

              // for (int i = 0; i <= array_size-nx; i++){
              //   double current_ny = i/nx;
              //   double current_cut = nx*trunc(current_ny)+idx_h-round(dx_diffusion/2);
              //   if (i>=current_cut && i <=nx*trunc(current_ny)+idx_h+round(dx_diffusion/2)){
              //     kd[i]=new_diffusion*(1-(cos((2*boost::math::double_constants::pi/dx_diffusion)*(i-current_cut))+1)/2);
              //   }
              //   else{
              //     kd[i]=kdd;
              //   }
              // }  
              // std::cout << " Use sinuosid function for Diffusion " <<std::endl; 
              fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);               
            }
            else{
              for (int i = 0; i <= array_size-nx; i++){
                double current_ny = i/nx;
                //std::cout << "current_ny :" << current_ny <<std::endl;
                if (i>=nx*trunc(current_ny)+idx_h-round(dx_diffusion/2) && i <=nx*trunc(current_ny)+idx_h+round(dx_diffusion/2)){
                  kd[i]=new_diffusion;
                  //std::cout << "new_diffusion :"<<std::endl; 
                }
                else{
                  kd[i]=kdd;
                  //std::cout << "old_diffusion :"<<std::endl; 
                }
              }  
              fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);               
            } 
          }
          else if (diffuse_area && new_diffusion_step != 0 && use_intervals )
          {
            if (current_timestep == new_diffusion_step)
            {
              this->get_pcout() << "   Initializing FastScape with new diffusion... " << (1 + surface_resolution + additional_refinement) << " levels, cell size: " << dx << " m." << std::endl;
              for (int i = 0; i <= array_size; i++)
              {
                if (test[i] == 1)
                  kd[i] = new_diffusion;
                else
                  kd[i] = kdd;
              }
            }
            fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);
          }

          fastscape_set_h_(h.get());

            //Initialize first time step, and update steps.
            int visualization_step = istep + restart_step;
            steps = istep + steps;

            this->get_pcout() << "   Calling FastScape... " << (steps - istep) << " timesteps of " << f_dt << " years." << std::endl;
            {
              auto t_start = std::chrono::high_resolution_clock::now();

              /*
                 * If we use stratigraphy it'll handle visualization and not the normal function.
                 * TODO: The frequency in this needs to be the same as the total timesteps fastscape will
                 * run for, need to figure out how to work this in better.
                 */
              if (use_strat && current_timestep == 1)
                fastscape_strati_(&nstepp, &nreflectorp, &steps, &vexp);
              else if (!use_strat)
                fastscape_named_vtk_(h.get(), &vexp, &visualization_step, c, &length);

              do
              {
                //execute step, this increases timestep counter
                fastscape_execute_step_();

                //get value of time step counter
                fastscape_get_step_(&istep);

                //outputs new h values
                fastscape_copy_h_(h.get());
              } while (istep < steps);

              auto t_end = std::chrono::high_resolution_clock::now();
              double r_time = std::chrono::duration<double>(t_end - t_start).count();
              keep_time += r_time;
              this->get_pcout() << "      FastScape runtime... " << round(r_time * 1000) / 1000 << "s" << std::endl;
              }

              //If we've reached the end time, destroy fastscape.
              if (this->get_time() + a_dt >= end_time)
              {
                this->get_pcout() << "   Destroying FastScape..." << std::endl;

                visualization_step = visualization_step + 1;

                if (!use_strat)
                  fastscape_named_vtk_(h.get(), &vexp, &visualization_step, c, &length);

                fastscape_destroy_();
              }

              //Find out our velocities from the change in height.
              for (int i = 0; i <= array_size; i++)
              {
                V[i] = (h[i] - h_old[i]) / a_dt;
              }

              MPI_Bcast(&V[0], array_size + 1, MPI_DOUBLE, 0, this->get_mpi_communicator());
            }
        else
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Fastscape 1 proc");

          for (unsigned int i = 0; i < temporary_variables.size(); i++)
            MPI_Ssend(&temporary_variables[i][0], temporary_variables[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          MPI_Bcast(&V[0], array_size + 1, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }

        TableIndices<dim> size_idx;
        for (unsigned int d = 0; d < dim; ++d)
        {
          size_idx[d] = table_intervals[d] + 1;
        }

        Table<dim, double> data_table;
        data_table.TableBase<dim, double>::reinit(size_idx);
        TableIndices<dim> idx;

        //this variable gives us how many slices near the boundaries to ignore,
        //this helps avoid boundary conditions effecting the topography.
        if (dim == 2)
        {
          std::vector<double> V2(nx);

          for (int i = 1; i < (nx - 1); i++)
          {
            if (slice)
            {
              int index = i + nx * (round((ny - 1) / 2));
              V2[i - 1] = V[index];
            }
            else
            {
              for (int ys = 1; ys < (ny - 1); ys++)
              {
                int index = i + nx * ys;
                V2[i - 1] += V[index];
              }
              V2[i - 1] = V2[i - 1] / (ny - 2);
            }
          }

          for (unsigned int i = 0; i < data_table.size()[1]; ++i)
          {
            idx[1] = i;

            for (unsigned int j = 0; j < (data_table.size()[0]); ++j)
            {
              idx[0] = j;

              //Convert from m/yr to m/s if needed
              if (i == 1)
              {
                if (this->convert_output_to_years())
                  data_table(idx) = V2[j] / year_in_seconds;
                else
                  data_table(idx) = V2[j];
              }
              else
                data_table(idx) = 0;
            }
          }
        }

        if (dim == 3)
        {
          //Indexes through z, y, and then x.
          for (unsigned int k = 0; k < data_table.size()[2]; ++k)
          {
            idx[2] = k;

            for (unsigned int i = 0; i < data_table.size()[1]; ++i)
            {
              idx[1] = i;

              for (unsigned int j = 0; j < data_table.size()[0]; ++j)
              {
                idx[0] = j;

                //We only care about the surface
                if (k == 1)
                {
                  if (this->convert_output_to_years())
                    data_table(idx) = V[(nx + 1) + nx * i + j] / year_in_seconds;
                  else
                    data_table(idx) = V[(nx + 1) + nx * i + j];
                }
                else
                  data_table(idx) = 0;
              }
            }
          }
        }

        Functions::InterpolatedUniformGridData<dim> *velocities;
        velocities = new Functions::InterpolatedUniformGridData<dim>(grid_extent,
                                                                     table_intervals,
                                                                     data_table);

        auto lambda = [&](const Point<dim> &p) -> double {
          return velocities->value(p);
        };

        VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
            lambda,
            dim - 1,
            dim);

        VectorTools::interpolate_boundary_values(mesh_deformation_dof_handler,
                                                 *boundary_ids.begin(),
                                                 vector_function_object,
                                                 mesh_velocity_constraints);
      }
    }

    //TODO: Give better explanations of variables and cite the fastscape documentation.
    template <int dim>
    void FastScape<dim>::declare_parameters(ParameterHandler &prm)
    {

      prm.declare_entry("End time",
                        /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max() /
                                              year_in_seconds) = */
                        "5.69e+300",
                        Patterns::Double(),
                        "The end time of the simulation. The default value is a number "
                        "so that when converted from years to seconds it is approximately "
                        "equal to the largest number representable in floating point "
                        "arithmetic. For all practical purposes, this equals infinity. "
                        "Units: Years if the "
                        "'Use years in output instead of seconds' parameter is set; "
                        "seconds otherwise.");

      prm.enter_subsection("Mesh deformation");
      {
        prm.enter_subsection("Fastscape");
        {
          prm.declare_entry("Number of steps", "10",
                            Patterns::Integer(),
                            "Number of steps per ASPECT timestep");
          prm.declare_entry("Maximum timestep", "10e3",
                            Patterns::Double(0),
                            "Maximum timestep for fastscape.");
          prm.declare_entry("Vertical exaggeration", "3",
                            Patterns::Double(),
                            "Vertical exaggeration for fastscape's VTK file.");
          prm.declare_entry("Additional fastscape refinement", "0",
                            Patterns::Integer(),
                            "Refinement level expected at surface to determine"
                            "proper nx and ny values");
          prm.declare_entry("Use center slice for 2d", "false",
                            Patterns::Bool(),
                            "If this is set to true, then a 2D model will only consider the "
                            "center slice fastscape gives. If set to false, then aspect will"
                            "average the inner third of what fastscape calculates.");
          prm.declare_entry("Fastscape seed", "1000",
                            Patterns::Integer(),
                            "Seed used for adding an initial 0-1 m noise to fastscape topography.");
          prm.declare_entry("Surface resolution", "1",
                            Patterns::Integer(),
                            "This should be set to the highest ASPECT resolution level you expect at the surface.");
          prm.declare_entry("Resolution difference", "0",
                            Patterns::Integer(),
                            "The difference between the lowest and highest resolution level at surface. So if three resolution "
                            "levels are expected, this would be set to 2.");
          prm.declare_entry("Use marine parameters", "false",
                            Patterns::Bool(),
                            "Flag to use marine parameters");
          prm.declare_entry("Use stratigraphy", "false",
                            Patterns::Bool(),
                            "Flag to use stratigraphy");
          prm.declare_entry("Total steps", "100000",
                            Patterns::Integer(),
                            "Total number of steps you expect in the FastScape model, only used if stratigraphy is turned on.");
          prm.declare_entry("Number of horizons", "1",
                            Patterns::Integer(),
                            "Number of horizons to track and visualize in FastScape..");
          prm.declare_entry("Use velocities", "true",
                            Patterns::Bool(),
                            "Flag on whether to have FastScape deal with velocities");
          prm.declare_entry("Minimum timestep to run fastscape", "0",
                            Patterns::Double(),
                            "Set a time if you don't want FastScape to run every timestep.");
          prm.declare_entry("Use box with lithosphere 2d", "true",
                            Patterns::Bool(),
                            "Flag on to autorize the user to enter its own model dimension X extent and XY repetitions");
          prm.declare_entry("Diffuse specific area", "false",
                            Patterns::Bool(),
                            "Flag on whether to diffuse a specific area or not, only work in 2D");



          prm.enter_subsection("Box with lithosphere 2d");
          {
          prm.declare_entry("X extent in 2d", "2208000",
                            Patterns::Double(),
                            "Set a X extent in meters.");  
          prm.declare_entry("Y extent in 2d", "100000",
                            Patterns::Double(),
                            "Y extent when used with 2D");                                                
          prm.declare_entry("X repetitions in 2d", "69",
                            Patterns::Double(),
                            "Set a X repetitions.");
          prm.declare_entry("Y repetitions in 2d", "35",
                            Patterns::Double(),
                            "Set a Y repetitions.");
          }    
          prm.leave_subsection();                              


          prm.enter_subsection("Specific area");
          {
            prm.declare_entry("Diffuse minimum topography", "true",
                              Patterns::Bool(),
                              "Flag on whether to diffuse depending on an interval given");
            prm.declare_entry("Use intervals", "false",
                              Patterns::Bool(),
                              "Flag on whether to diffuse depending on an interval given");
            prm.declare_entry("Diffusion interval x1", "90000",
                            Patterns::Double(),
                            "Set a the minimum X value to diffuse a length interval.");
            prm.declare_entry("Diffusion interval x2", "140000",
                            Patterns::Double(),
                            "Set a the maximum X value to diffuse a length interval.");
            prm.declare_entry("Specific bedrock diffusivity", "5000",
                              Patterns::Double(),
                              "Diffusivity of a specific bedrock area.");
            prm.declare_entry("Start diffusing at timestep", "20",
                              Patterns::Double(),
                              "timestep when the specific area should start diffusing."); 
            prm.declare_entry("Number of dx to diffuse", "10",
                            Patterns::Double(),
                            "Set a the number of dx to diffuse near the minimum topography."); 
            prm.declare_entry("Use sinusoid diffusion", "true",
                              Patterns::Bool(),
                              "sinuosid function of diffusion. The minimum topography is diffuse at the specific diffusion and decrease to reach normal Kd away ");  
            prm.declare_entry("Start diffusing", "100e3",
                              Patterns::Double(),
                              "timestep when the specific area should start diffusing.");                                                                                                                                    
          }
          prm.leave_subsection();

          

          prm.enter_subsection("Boundary conditions");
          {
            prm.declare_entry("Bottom", "1",
                              Patterns::Integer(0, 1),
                              "Bottom boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry("Right", "1",
                              Patterns::Integer(0, 1),
                              "Right boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry("Top", "1",
                              Patterns::Integer(0, 1),
                              "Top boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry("Left", "1",
                              Patterns::Integer(0, 1),
                              "Left boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry("Left mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through left boundary (m^2/yr)");
            prm.declare_entry("Right mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through right boundary (m^2/yr)");
            prm.declare_entry("Top mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through top boundary (m^2/yr)");
            prm.declare_entry("Bottom mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through bottom boundary (m^2/yr)");
          }
          prm.leave_subsection();

          prm.enter_subsection("Erosional parameters");
          {
            prm.declare_entry("Drainage area exponent", "0.4",
                              Patterns::Double(),
                              "Exponent for drainage area.");
            prm.declare_entry("Slope exponent", "1",
                              Patterns::Double(),
                              "The  slope  exponent  for  SPL (n).  Generally  m/n  should  equal  approximately 0.4");
            prm.declare_entry("Multi-direction slope exponent", "1",
                              Patterns::Double(),
                              "Exponent to determine the distribution from the SPL to neighbor nodes, with"
                              "10 being steepest decent and 1 being more varied.");
            prm.declare_entry("Bedrock deposition coefficient", "-1",
                              Patterns::Double(),
                              "Deposition coefficient for bedrock");
            prm.declare_entry("Sediment deposition coefficient", "-1",
                              Patterns::Double(),
                              "Deposition coefficient for sediment");
            prm.declare_entry("Bedrock river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for bedrock");
            prm.declare_entry("Sediment river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for sediment ");
            prm.declare_entry("Bedrock diffusivity", "1",
                              Patterns::Double(),
                              "Diffusivity of bedrock.");
            prm.declare_entry("Sediment diffusivity", "-1",
                              Patterns::Double(),
                              "Diffusivity of sediment.");
          }
          prm.leave_subsection();

          prm.enter_subsection("Marine parameters");
          {
            prm.declare_entry("Sea level", "0",
                              Patterns::Double(),
                              "Sea level in meters. ");
            prm.declare_entry("Sand porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of sand. ");
            prm.declare_entry("Shale porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of shale. ");
            prm.declare_entry("Sand e-folding depth", "1e3",
                              Patterns::Double(),
                              "e-folding depth for the exponential of the sand porosity law.");
            prm.declare_entry("Shale e-folding depth", "1e3",
                              Patterns::Double(),
                              "e-folding depth for the exponential of the shale porosity law.");
            prm.declare_entry("Sand-shale ratio", "0.5",
                              Patterns::Double(),
                              "ratio of sand to shale for material leaving continent.");
            prm.declare_entry("Depth averaging thickness", "1e2",
                              Patterns::Double(),
                              "Depth averaging for the sand-shale equation in meters.");
            prm.declare_entry("Sand transport coefficient", "5e2",
                              Patterns::Double(),
                              "Transport coefficient for sand in m^2/yr ");
            prm.declare_entry("Shale transport coefficient", "2.5e2",
                              Patterns::Double(),
                              "Transport coefficient for shale in m^/2yr");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Mesh refinement");
      {
        prm.declare_entry("Initial global refinement", "2",
                          Patterns::Integer(0),
                          "The number of global refinement steps performed on "
                          "the initial coarse mesh, before the problem is first "
                          "solved there.");
      }
      prm.leave_subsection();
    }

    template <int dim>
    void FastScape<dim>::parse_parameters(ParameterHandler &prm)
    {
      end_time = prm.get_double("End time");
      if (prm.get_bool("Use years in output instead of seconds") == true)
        end_time *= year_in_seconds;
      prm.enter_subsection("Mesh deformation");
      {
        prm.enter_subsection("Fastscape");
        {
          nstep = prm.get_integer("Number of steps");
          max_timestep = prm.get_double("Maximum timestep");
          vexp = prm.get_double("Vertical exaggeration");
          additional_refinement = prm.get_integer("Additional fastscape refinement");
          slice = prm.get_bool("Use center slice for 2d");
          fs_seed = prm.get_integer("Fastscape seed");
          surface_resolution = prm.get_integer("Surface resolution");
          resolution_difference = prm.get_integer("Resolution difference");
          use_marine = prm.get_bool("Use marine parameters");
          use_strat = prm.get_bool("Use stratigraphy");
          nstepp = prm.get_integer("Total steps");
          nreflectorp = prm.get_integer("Number of horizons");
          use_velocity = prm.get_bool("Use velocities");
          minimum_time = prm.get_double("Minimum timestep to run fastscape");
          use_boxlitho_2d = prm.get_bool("Use box with lithosphere 2d");
          diffuse_area = prm.get_bool("Diffuse specific area");


          prm.enter_subsection("Box with lithosphere 2d");
          {
            x_extent_2d = prm.get_double("X extent in 2d");
            y_extent_2d = prm.get_double("Y extent in 2d");            
            x_repetitions_2d = prm.get_double("X repetitions in 2d");
            y_repetitions_2d = prm.get_double("Y repetitions in 2d");
          }
          prm.leave_subsection();

          prm.enter_subsection("Specific area");
          {
            minimum_topography_diffusion = prm.get_bool("Diffuse minimum topography");
            use_intervals = prm.get_bool("Use intervals");
            x_min_diffusion = prm.get_double("Diffusion interval x1");
            x_max_diffusion = prm.get_double("Diffusion interval x2");
            new_diffusion = prm.get_double("Specific bedrock diffusivity");
            new_diffusion_step = prm.get_double("Start diffusing at timestep");
            time_to_diffuse = prm.get_double("Start diffusing");            
            dx_diffusion = prm.get_double("Number of dx to diffuse");
            diffusion_sinusoid = prm.get_bool("Use sinusoid diffusion");
          }
          prm.leave_subsection();

          prm.enter_subsection("Boundary conditions");
          {
            bottom = prm.get_integer("Bottom");
            right = prm.get_integer("Right");
            top = prm.get_integer("Top");
            left = prm.get_integer("Left");
            left_flux = prm.get_double("Left mass flux");
            right_flux = prm.get_double("Right mass flux");
            top_flux = prm.get_double("Top mass flux");
            bottom_flux = prm.get_double("Bottom mass flux");

            bc = bottom * 1000 + right * 100 + top * 10 + left;

            if ((left_flux != 0 && top_flux != 0) || (left_flux != 0 && bottom_flux != 0) ||
                (right_flux != 0 && bottom_flux != 0) || (right_flux != 0 && top_flux != 0))
              AssertThrow(false, ExcMessage("Currently the plugin does not support mass flux through adjacent boundaries."));
          }
          prm.leave_subsection();

          prm.enter_subsection("Erosional parameters");
          {
            m = prm.get_double("Drainage area exponent");
            n = prm.get_double("Slope exponent");
            kfsed = prm.get_double("Sediment river incision rate");
            kff = prm.get_double("Bedrock river incision rate");
            kdsed = prm.get_double("Sediment diffusivity");
            kdd = prm.get_double("Bedrock diffusivity");
            g = prm.get_double("Bedrock deposition coefficient");
            gsed = prm.get_double("Sediment deposition coefficient");
            p = prm.get_double("Multi-direction slope exponent");
          }
          prm.leave_subsection();

          prm.enter_subsection("Marine parameters");
          {
            sl = prm.get_double("Sea level");
            p1 = prm.get_double("Sand porosity");
            p2 = prm.get_double("Shale porosity");
            z1 = prm.get_double("Sand e-folding depth");
            z2 = prm.get_double("Shale e-folding depth");
            r = prm.get_double("Sand-shale ratio");
            l = prm.get_double("Depth averaging thickness");
            kds1 = prm.get_double("Sand transport coefficient");
            kds2 = prm.get_double("Shale transport coefficient");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Mesh refinement");
      {
        initial_global_refinement = prm.get_integer("Initial global refinement");
      }
      prm.leave_subsection();
    }
  } // namespace MeshDeformation
} // namespace aspect

// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FastScape,
                                           "fastscape",
                                           "A plugin which uses the program FastScape to add surface processes "
                                           "such as diffusion, sedimentation, and the Stream Power Law to ASPECT. "
                                           "FastScape is initialized with ASPECT height and velocities, and then "
                                           "continues to run in the background, updating with new ASPECT velocities "
                                           "when called. Once FastScape has run for the given amount of timesteps, it "
                                           "compares the initial and final heights to send a velocity back to the mesh "
                                           "deformation plugin.")
  }
} // namespace aspect
