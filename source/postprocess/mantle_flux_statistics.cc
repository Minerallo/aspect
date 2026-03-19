#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/global.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <queue>
#include <set>
#include <limits>
#include <algorithm>
#include <cmath>

namespace aspect
{
namespace Postprocess
{

template <int dim>
double spherical_distance(const Point<dim> &p1,
                          const Point<dim> &p2)
{
  const double R1 = p1.norm();
  const double R2 = p2.norm();
  const double R  = 0.5 * (R1 + R2);

  double dot = 0.0;
  for (unsigned int d=0; d<dim; ++d)
    dot += p1[d] * p2[d];

  const double cos_angle = std::max(-1.0,
                                    std::min(1.0, dot / (R1 * R2)));

  return R * std::acos(cos_angle);
}



template <int dim>
class MantleFluxStatistics :
  public Interface<dim>,
  public SimulatorAccess<dim>
{
public:
  std::pair<std::string,std::string>
  execute(TableHandler &statistics) override;

  static void declare_parameters(ParameterHandler &prm);
  void parse_parameters(ParameterHandler &prm) override;

private:
  std::vector<double> depths;
  double depth_tolerance;
  double slab_temperature_threshold;
  double plume_temperature_threshold;
  double minimum_slab_length;
  double cluster_radius_factor;

  unsigned int minimum_cluster_size;

  bool write_to_file;
  double output_interval;
  double last_output_time = std::numeric_limits<double>::quiet_NaN();
};



template <int dim>
void MantleFluxStatistics<dim>::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Postprocess");
  prm.enter_subsection("Mantle flux statistics");
  {
    prm.declare_entry("Depths", "440e3",
                      Patterns::List(Patterns::Double()),
                      "Depths (m) where fluxes are evaluated.");

    prm.declare_entry("Depth tolerance", "20000",
                      Patterns::Double(),
                      "Half thickness of shell used for integration.");

    prm.declare_entry("Slab temperature threshold", "-200",
                      Patterns::Double(),
                      "Nonadiabatic temperature threshold for slab detection.");

    prm.declare_entry("Plume temperature threshold", "200",
                      Patterns::Double(),
                      "Nonadiabatic temperature threshold for plume detection.");

    prm.declare_entry("Minimum slab length", "0",
                      Patterns::Double(),
                      "Minimum spherical arc length (m) for a connected cold "
                      "downwelling cluster to be counted as a slab. "
                      "Useful to filter out small drip-like downwellings.");

    prm.declare_entry("Cluster radius factor","4.0",
                      Patterns::Double(),
                      "Cluster linking distance in units of depth tolerance.");

    prm.declare_entry("Minimum cluster size","5",
                      Patterns::Integer(1),
                      "Minimum number of points required for a cluster.");

    prm.declare_entry("Output to file", "false",
                      Patterns::Bool(),
                      "Whether or not to write detected plume and slab cluster "
                      "points to a text file named 'mantle_flux_clusters.NNNNN' "
                      "in the output directory.");

    prm.declare_entry("Time between text output", "0.",
                      Patterns::Double(0.),
                      "The time interval between each generation of text output "
                      "files. A value of zero indicates that output should be "
                      "generated in each time step. Units: years if "
                      "'Use years in output instead of seconds' is set; "
                      "seconds otherwise.");
  }
  prm.leave_subsection();
  prm.leave_subsection();
}



template <int dim>
void MantleFluxStatistics<dim>::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Postprocess");
  prm.enter_subsection("Mantle flux statistics");
  {
    depths =
      Utilities::string_to_double(
        Utilities::split_string_list(prm.get("Depths")));

    depth_tolerance = prm.get_double("Depth tolerance");
    slab_temperature_threshold = prm.get_double("Slab temperature threshold");
    plume_temperature_threshold = prm.get_double("Plume temperature threshold");
    minimum_slab_length = prm.get_double("Minimum slab length");
    cluster_radius_factor = prm.get_double("Cluster radius factor");
    minimum_cluster_size = prm.get_integer("Minimum cluster size");

    write_to_file = prm.get_bool("Output to file");
    output_interval = prm.get_double("Time between text output");
    if (this->convert_output_to_years())
      output_interval *= year_in_seconds;
  }
  prm.leave_subsection();
  prm.leave_subsection();
}



template <int dim>
std::pair<std::string,std::string>
MantleFluxStatistics<dim>::execute(TableHandler &statistics)
{
  const MPI_Comm comm = this->get_mpi_communicator();
  const unsigned int rank = Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(comm);

  const double m3s_to_km3yr = year_in_seconds / 1e9;

  const Quadrature<dim> &quadrature =
    this->introspection().quadratures.velocities;

  FEValues<dim> fe_values(this->get_mapping(),
                          this->get_fe(),
                          quadrature,
                          update_values |
                          update_quadrature_points |
                          update_JxW_values);

  const unsigned int n_q = quadrature.size();

  std::vector<Tensor<1,dim>> velocity_values(n_q);
  std::vector<double> temperature_values(n_q);

  std::ostringstream screen_text;
  screen_text.precision(4);

  std::ostringstream output_file;
  if (rank == 0)
    output_file << "# x y z depth_km type cluster_id\n";

  unsigned int global_cluster_id = 0;

  for (const double depth : depths)
  {
    const double shell_thickness = 2.0 * depth_tolerance;

    double local_slab_flux = 0.0;
    double local_plume_flux = 0.0;
    double local_slab_bflux = 0.0;
    double local_plume_bflux = 0.0;

    std::vector<Point<dim>> local_plume_points;
    std::vector<Point<dim>> local_slab_points;

    for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        fe_values[this->introspection().extractors.velocities]
          .get_function_values(this->get_solution(), velocity_values);

        fe_values[this->introspection().extractors.temperature]
          .get_function_values(this->get_solution(), temperature_values);

        const auto &points = fe_values.get_quadrature_points();

        bool cell_has_plume = false;
        bool cell_has_slab  = false;

        for (unsigned int q=0; q<n_q; ++q)
        {
          const Point<dim> &x = points[q];

          const double current_depth =
            this->get_geometry_model().depth(x);

          if (std::abs(current_depth - depth) > depth_tolerance)
            continue;

          const Tensor<1,dim> r_hat = x / x.norm();
          const double vr = velocity_values[q] * r_hat;

          const double Tadi =
            this->get_adiabatic_conditions().temperature(x);

          const double Tprime = temperature_values[q] - Tadi;
          const double dA = fe_values.JxW(q) / shell_thickness;

          if (Tprime < slab_temperature_threshold && vr < 0.0)
          {
            local_slab_flux += (-vr) * dA;
            local_slab_bflux += ((-Tprime) * (-vr)) * dA;
            cell_has_slab = true;
          }

          if (Tprime > plume_temperature_threshold && vr > 0.0)
          {
            local_plume_flux += vr * dA;
            local_plume_bflux += (Tprime * vr) * dA;
            cell_has_plume = true;
          }
        }

        if (cell_has_plume)
          local_plume_points.push_back(cell->center());

        if (cell_has_slab)
          local_slab_points.push_back(cell->center());
      }

    const double slab_flux =
      Utilities::MPI::sum(local_slab_flux, comm);

    const double plume_flux =
      Utilities::MPI::sum(local_plume_flux, comm);

    const double slab_bflux =
      Utilities::MPI::sum(local_slab_bflux, comm);

    const double plume_bflux =
      Utilities::MPI::sum(local_plume_bflux, comm);

    if (rank != 0)
    {
      Utilities::MPI::isend(local_plume_points, comm, 0, 100);
      Utilities::MPI::isend(local_slab_points,  comm, 0, 101);
    }

    unsigned int plume_count = 0;
    unsigned int slab_count  = 0;
    double total_slab_length = 0.0;
    double average_slab_length = 0.0;
    double average_plume_radius = 0.0;

    if (rank == 0)
    {
      std::vector<Point<dim>> all_plumes = local_plume_points;
      std::vector<Point<dim>> all_slabs  = local_slab_points;

      for (unsigned int p=1; p<n_procs; ++p)
      {
        auto recv_plumes =
          Utilities::MPI::irecv<std::vector<Point<dim>>>(comm, p, 100).get();

        auto recv_slabs =
          Utilities::MPI::irecv<std::vector<Point<dim>>>(comm, p, 101).get();

        all_plumes.insert(all_plumes.end(),
                          recv_plumes.begin(), recv_plumes.end());

        all_slabs.insert(all_slabs.end(),
                         recv_slabs.begin(), recv_slabs.end());
      }

      const double cluster_radius = cluster_radius_factor * depth_tolerance;

      // Plume clustering
      {
        std::vector<bool> visited(all_plumes.size(), false);
        double total_plume_radius = 0.0;

        for (unsigned int i=0; i<all_plumes.size(); ++i)
        {
          if (visited[i])
            continue;

          std::queue<unsigned int> q;
          q.push(i);

          std::vector<Point<dim>> cluster_points;

          while (!q.empty())
          {
            const unsigned int j = q.front();
            q.pop();

            if (visited[j])
              continue;

            visited[j] = true;
            cluster_points.push_back(all_plumes[j]);

            for (unsigned int k=0; k<all_plumes.size(); ++k)
              if (!visited[k] &&
                  spherical_distance(all_plumes[j], all_plumes[k]) < cluster_radius)
                q.push(k);
          }

          if (cluster_points.size() < minimum_cluster_size)
            continue;

          const unsigned int cluster_id = global_cluster_id++;
          ++plume_count;

          for (const auto &p : cluster_points)
            output_file << p[0] << ' '
                        << p[1] << ' '
                        << p[2] << ' '
                        << depth/1000.0 << ' '
                        << "plume" << ' '
                        << cluster_id << '\n';

          Point<dim> center;
          for (const auto &p : cluster_points)
            center += p;

          center /= static_cast<double>(cluster_points.size());

          const double mean_radius = cluster_points[0].norm();
          if (center.norm() > 0.0)
            center *= (mean_radius / center.norm());

          double plume_radius = 0.0;
          for (const auto &p : cluster_points)
            plume_radius += spherical_distance(center, p);

          plume_radius /= static_cast<double>(cluster_points.size());

          total_plume_radius += plume_radius;
        }

        if (plume_count > 0)
          average_plume_radius = total_plume_radius / static_cast<double>(plume_count);
      }

      // Slab clustering
      {
        std::vector<bool> visited(all_slabs.size(), false);
        double total_valid_slab_length = 0.0;

        for (unsigned int i=0; i<all_slabs.size(); ++i)
        {
          if (visited[i])
            continue;

          std::queue<unsigned int> q;
          q.push(i);

          std::vector<Point<dim>> cluster_points;

          while (!q.empty())
          {
            const unsigned int j = q.front();
            q.pop();

            if (visited[j])
              continue;

            visited[j] = true;
            cluster_points.push_back(all_slabs[j]);

            for (unsigned int k=0; k<all_slabs.size(); ++k)
              if (!visited[k] &&
                  spherical_distance(all_slabs[j], all_slabs[k]) < cluster_radius)
                q.push(k);
          }

          if (cluster_points.size() < minimum_cluster_size)
            continue;

          const unsigned int N = cluster_points.size();
          std::vector<std::vector<unsigned int>> neighbors(N);

          for (unsigned int a=0; a<N; ++a)
            for (unsigned int b=a+1; b<N; ++b)
              if (spherical_distance(cluster_points[a],cluster_points[b]) < cluster_radius)
              {
                neighbors[a].push_back(b);
                neighbors[b].push_back(a);
              }

          auto bfs_farthest = [&](unsigned int start)
          {
            std::vector<bool> visited_bfs(N,false);
            std::queue<unsigned int> q;

            std::vector<double> dist(N,0.0);

            q.push(start);
            visited_bfs[start] = true;

            unsigned int farthest = start;

            while(!q.empty())
            {
              unsigned int u = q.front();
              q.pop();

              for (auto v : neighbors[u])
                if (!visited_bfs[v])
                {
                  visited_bfs[v] = true;

                  dist[v] = dist[u] +
                    spherical_distance(cluster_points[u],cluster_points[v]);

                  q.push(v);

                  if (dist[v] > dist[farthest])
                    farthest = v;
                }
            }

            return std::make_pair(farthest,dist[farthest]);
          };

          auto p1 = bfs_farthest(0);
          auto p2 = bfs_farthest(p1.first);

          double slab_length = p2.second;

          if (slab_length >= minimum_slab_length)
          {
            const unsigned int cluster_id = global_cluster_id++;
            ++slab_count;
            total_valid_slab_length += slab_length;

            for (const auto &p : cluster_points)
              output_file << p[0] << ' '
                          << p[1] << ' '
                          << p[2] << ' '
                          << depth/1000.0 << ' '
                          << "slab" << ' '
                          << cluster_id << '\n';
          }
        }

        total_slab_length = total_valid_slab_length;

        if (slab_count > 0)
          average_slab_length = total_valid_slab_length / static_cast<double>(slab_count);
      }
    }

    plume_count = Utilities::MPI::broadcast(comm, plume_count, 0);
    slab_count  = Utilities::MPI::broadcast(comm, slab_count, 0);
    total_slab_length = Utilities::MPI::broadcast(comm, total_slab_length, 0);
    average_slab_length = Utilities::MPI::broadcast(comm, average_slab_length, 0);
    average_plume_radius = Utilities::MPI::broadcast(comm, average_plume_radius, 0);

    const std::string depth_label =
      Utilities::to_string(depth/1e3) + " km";

    const std::string plume_flux_name =
      "Plume flux depth " + depth_label + " (km^3/yr)";
    statistics.add_value(plume_flux_name, plume_flux * m3s_to_km3yr);
    statistics.set_precision(plume_flux_name, 8);
    statistics.set_scientific(plume_flux_name, true);

    const std::string slab_flux_name =
      "Slab flux depth " + depth_label + " (km^3/yr)";
    statistics.add_value(slab_flux_name, slab_flux * m3s_to_km3yr);
    statistics.set_precision(slab_flux_name, 8);
    statistics.set_scientific(slab_flux_name, true);

    const std::string plume_bflux_name =
      "Plume buoyancy flux depth " + depth_label + " (K km^3/yr)";
    statistics.add_value(plume_bflux_name, plume_bflux * m3s_to_km3yr);
    statistics.set_precision(plume_bflux_name, 8);
    statistics.set_scientific(plume_bflux_name, true);

    const std::string slab_bflux_name =
      "Slab buoyancy flux depth " + depth_label + " (K km^3/yr)";
    statistics.add_value(slab_bflux_name, slab_bflux * m3s_to_km3yr);
    statistics.set_precision(slab_bflux_name, 8);
    statistics.set_scientific(slab_bflux_name, true);

    const std::string plume_count_name =
      "Number of plumes depth " + depth_label;
    statistics.add_value(plume_count_name, plume_count);

    const std::string slab_count_name =
      "Number of slabs depth " + depth_label;
    statistics.add_value(slab_count_name, slab_count);

    const std::string total_slab_length_name =
      "Total slab length depth " + depth_label + " (km)";
    statistics.add_value(total_slab_length_name, total_slab_length / 1000.0);
    statistics.set_precision(total_slab_length_name, 8);
    statistics.set_scientific(total_slab_length_name, true);

    const std::string average_slab_length_name =
      "Average slab length depth " + depth_label + " (km)";
    statistics.add_value(average_slab_length_name, average_slab_length / 1000.0);
    statistics.set_precision(average_slab_length_name, 8);
    statistics.set_scientific(average_slab_length_name, true);

    const std::string average_plume_radius_name =
      "Average plume radius depth " + depth_label + " (km)";
    statistics.add_value(average_plume_radius_name, average_plume_radius / 1000.0);
    statistics.set_precision(average_plume_radius_name, 8);
    statistics.set_scientific(average_plume_radius_name, true);

    screen_text
      << "Depth " << depth_label
      << " plume_flux=" << plume_flux * m3s_to_km3yr << " km^3/yr"
      << " slab_flux=" << slab_flux * m3s_to_km3yr << " km^3/yr"
      << " plume_bflux=" << plume_bflux * m3s_to_km3yr << " K km^3/yr"
      << " slab_bflux=" << slab_bflux * m3s_to_km3yr << " K km^3/yr"
      << " plumes=" << plume_count
      << " slabs=" << slab_count
      << " total_slab_length=" << total_slab_length / 1000.0 << " km"
      << " avg_slab_length=" << average_slab_length / 1000.0 << " km"
      << " avg_plume_radius=" << average_plume_radius / 1000.0 << " km"
      << " ";
  }

  if (std::isnan(last_output_time))
    last_output_time = this->get_time() - output_interval;

  if (write_to_file &&
      !((this->get_time() < last_output_time + output_interval)
         && (this->get_timestep_number() != 0)))
  {
    std::string filename = this->get_output_directory() +
                           "mantle_flux_clusters." +
                           Utilities::int_to_string(this->get_timestep_number(), 5);

    if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
      filename.append("." + Utilities::int_to_string(this->get_nonlinear_iteration(), 4));

    Utilities::collect_and_write_file_content(filename,
                                              output_file.str(),
                                              this->get_mpi_communicator());

    if (output_interval > 0)
    {
      const double magic = 1.0 + 2.0*std::numeric_limits<double>::epsilon();
      last_output_time =
        last_output_time +
        std::floor((this->get_time()-last_output_time)/output_interval*magic)
        * output_interval/magic;
    }
  }

  return std::pair<std::string,std::string>(
    "Mantle flux statistics:",
    screen_text.str());
}



ASPECT_REGISTER_POSTPROCESSOR(
  MantleFluxStatistics,
  "mantle flux statistics",
  "Computes slab/plume fluxes, buoyancy fluxes, counts plume and slab "
  "clusters using MPI gather, and reports total slab length, average slab "
  "length, and average plume radius using spherical arc distances. "
  "If 'Output to file' is set to true, this postprocessor also writes "
  "plume and slab cluster points into text files named "
  "'mantle_flux_clusters.NNNNN' in the output directory.")
}
}