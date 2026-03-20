#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>
#include <algorithm>

namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      /**
       * Compute Euclidean distance between two points.
       *
       * Used when estimating the length of an edge on the surface mesh
       * for the continent perimeter calculation.
       */
      template <int dim>
      double point_distance(const Point<dim> &a, const Point<dim> &b)
      {
        return a.distance(b);
      }


      /**
       * Convert a point into a string key after quantizing each coordinate.
       *
       * Why?
       * ----
       * When we reconstruct the surface mesh connectivity on rank 0,
       * the same geometric vertex may appear several times from different cells
       * and even different MPI ranks. We need a robust way to identify that
       * those vertices are the same.
       *
       * The idea is:
       *   - divide each coordinate by a tiny tolerance
       *   - round it
       *   - store the rounded integers as a string
       *
       * This is not mathematically perfect, but it is very practical and
       * commonly sufficient for mesh connectivity bookkeeping.
       */
      template <int dim>
      std::string point_key(const Point<dim> &p,
                            const double tol = 1e-10)
      {
        std::ostringstream out;
        out << std::fixed << std::setprecision(0);

        for (unsigned int d = 0; d < dim; ++d)
          {
            const double q = std::round(p[d] / tol);
            out << q;
            if (d + 1 < dim)
              out << "_";
          }

        return out.str();
      }


      /**
       * Convert Cartesian coordinates to longitude/latitude in degrees.
       *
       * Only meaningful for dim=3 spherical geometry.
       * The compile-time branch `if constexpr (dim == 3)` must be used
       * before calling this function.
       */
      template <int dim>
      std::pair<double,double> cartesian_to_lon_lat_deg(const Point<dim> &p)
      {
        static_assert(dim == 3, "lon/lat conversion only implemented for dim=3.");

        const double r = p.norm();
        if (r <= 0.0)
          return std::make_pair(0.0, 0.0);

        const double lon = std::atan2(p[1], p[0]) * 180.0 / numbers::PI;
        const double lat = std::asin(p[2] / r) * 180.0 / numbers::PI;
        return std::make_pair(lon, lat);
      }


      /**
       * Simple disjoint-set / union-find structure.
       *
       * Purpose:
       * --------
       * We classify each top surface face as continent or not.
       * Then we want to group neighboring continent faces into connected blocks.
       *
       * DSU lets us merge connected faces efficiently.
       */
      struct DSU
      {
        std::vector<unsigned int> parent;
        std::vector<unsigned int> rank;

        explicit DSU(const unsigned int n)
          : parent(n), rank(n, 0)
        {
          for (unsigned int i = 0; i < n; ++i)
            parent[i] = i;
        }

        unsigned int find(unsigned int x)
        {
          if (parent[x] != x)
            parent[x] = find(parent[x]);
          return parent[x];
        }

        void unite(unsigned int a, unsigned int b)
        {
          a = find(a);
          b = find(b);

          if (a == b)
            return;

          if (rank[a] < rank[b])
            std::swap(a, b);

          parent[b] = a;

          if (rank[a] == rank[b])
            ++rank[a];
        }
      };


      /**
       * Summary information for one top-boundary face.
       *
       * Instead of storing all quadrature-point values forever,
       * we integrate them and keep one compact face record.
       */
      template <int dim>
      struct FaceRecord
      {
        Point<dim> center;                  // face center
        std::vector<Point<dim>> vertices;   // face vertices for connectivity/perimeter

        double area = 0.0;                  // face area
        double continent_value = 0.0;       // area-averaged sum of selected continent fields
        double speed = 0.0;                 // area-averaged speed magnitude on the face

        bool is_continent = false;          // thresholded continent mask
        int block_id = -1;                  // assigned later after connectivity analysis
      };


      /**
       * Summary information for one connected continent block.
       */
      template <int dim>
      struct BlockSummary
      {
        unsigned int id = numbers::invalid_unsigned_int;
        double area = 0.0;
        double speed_area_integral = 0.0;
        unsigned int n_faces = 0;
        Tensor<1,dim> centroid_numerator;
      };
    }



    /**
     * Postprocessor computing continent statistics from user-selected
     * compositional fields on the top surface.
     *
     * It outputs:
     *   - continental area
     *   - continent perimeter
     *   - number of continental blocks
     *   - largest continental block area
     *   - fragmentation index
     *   - average continental drift speed
     *
     * It can also write debug files:
     *   - one per top-surface face
     *   - one per connected continent block
     */
    template <int dim>
    class ContinentStatistics :
      public Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string>
        execute(TableHandler &statistics) override;

        void
        parse_parameters(ParameterHandler &prm) override;

        static
        void
        declare_parameters(ParameterHandler &prm);

      private:
        /**
         * Names of user-selected compositional fields representing "continent".
         *
         * Example:
         *   continent
         *   upper_crust, lower_crust
         *   upper_crust, lower_crust, continental_lithosphere
         */
        std::vector<std::string> continent_field_names;

        /**
         * Internal indices corresponding to the names above.
         */
        std::vector<unsigned int> continent_field_indices;

        /**
         * Threshold on the summed continent composition.
         *
         * A face is continental if:
         *   sum(selected fields) > continent_threshold
         */
        double continent_threshold = 0.5;

        /**
         * Ignore continent fragments smaller than this area [m^2].
         *
         * Useful to remove tiny numerical speckles from the block count.
         */
        double minimum_block_area = 0.0;

        /**
         * Whether to print a one-line summary in the screen output.
         */
        bool output_verbose_screen_line = true;

        /**
         * Whether to write a per-face debug file.
         */
        bool write_surface_map = true;

        /**
         * Whether to write a per-block summary file.
         */
        bool write_block_summary = true;

        /**
         * Prefix for the debug output files.
         */
        std::string output_file_prefix = "continent_statistics";

        /**
         * Convert field names to field indices.
         */
        void resolve_field_indices();

        /**
         * Check whether a face belongs to the top boundary.
         *
         * Check whether a face belongs to the top boundary
         * by comparing its boundary id to the geometry model's
         * symbolic "top" boundary id.
         */
        bool is_top_boundary_face(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  const unsigned int face_no) const;

        /**
         * Write one line per top surface face.
         */
        void write_surface_file(const std::vector<FaceRecord<dim>> &all_faces) const;

        /**
         * Write one line per connected continent block.
         */
        void write_block_file(const std::vector<BlockSummary<dim>> &blocks,
                              const std::vector<Point<dim>> &block_centroids) const;
    };



    template <int dim>
    void ContinentStatistics<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Continent statistics");
        {
          prm.declare_entry("Continent field names",
                            "",
                            Patterns::List(Patterns::Anything()),
                            "Comma-separated list of compositional field names "
                            "that should be interpreted as continental material.");

          prm.declare_entry("Continent threshold",
                            "0.5",
                            Patterns::Double(0.0),
                            "A top-surface face is considered continental if the "
                            "sum of the selected compositional fields exceeds this threshold.");

          prm.declare_entry("Minimum block area",
                            "0.0",
                            Patterns::Double(0.0),
                            "Ignore connected continental blocks with area smaller "
                            "than this value [m^2].");

          prm.declare_entry("Output verbose screen line",
                            "true",
                            Patterns::Bool(),
                            "Whether to print a compact summary line to the screen.");

          prm.declare_entry("Write surface map",
                            "true",
                            Patterns::Bool(),
                            "Write a per-face diagnostic file for the top surface.");

          prm.declare_entry("Write block summary",
                            "true",
                            Patterns::Bool(),
                            "Write a per-block diagnostic summary file.");

          prm.declare_entry("Output file prefix",
                            "continent_statistics",
                            Patterns::Anything(),
                            "Prefix used for the diagnostic output file names.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void ContinentStatistics<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Continent statistics");
        {
          continent_field_names =
            Utilities::split_string_list(prm.get("Continent field names"));

          continent_threshold = prm.get_double("Continent threshold");
          minimum_block_area = prm.get_double("Minimum block area");
          output_verbose_screen_line = prm.get_bool("Output verbose screen line");

          write_surface_map = prm.get_bool("Write surface map");
          write_block_summary = prm.get_bool("Write block summary");
          output_file_prefix = prm.get("Output file prefix");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Resolve immediately so bad field names fail early.
      resolve_field_indices();
    }



    template <int dim>
    void ContinentStatistics<dim>::resolve_field_indices()
    {
      continent_field_indices.clear();

      const std::vector<std::string> &all_names =
        this->introspection().chemical_composition_field_names();

      for (const std::string &requested_name : continent_field_names)
        {
          bool found = false;

          for (unsigned int i = 0; i < all_names.size(); ++i)
            {
              if (all_names[i] == requested_name)
                {
                  continent_field_indices.push_back(i);
                  found = true;
                  break;
                }
            }

          AssertThrow(found,
                      ExcMessage("Continent statistics: field <" +
                                 requested_name +
                                 "> was not found among the compositional fields."));
        }

      AssertThrow(!continent_field_indices.empty(),
                  ExcMessage("Continent statistics: at least one continent field "
                             "must be provided."));
    }



    template <int dim>
    bool ContinentStatistics<dim>::is_top_boundary_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int face_no) const
    {
    // First make sure this face is actually on a boundary.
    if (!cell->face(face_no)->at_boundary())
        return false;

    // Ask the geometry model which boundary id corresponds to the symbolic
    // name "top", then compare it to the face boundary id.
    const types::boundary_id top_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

    return (cell->face(face_no)->boundary_id() == top_id);
    }


    template <int dim>
    void ContinentStatistics<dim>::write_surface_file(const std::vector<FaceRecord<dim>> &all_faces) const
    {
      std::ostringstream filename;
      filename << this->get_output_directory()
               << output_file_prefix
               << "_surface."
               << std::setw(5) << std::setfill('0')
               << this->get_timestep_number();

      std::ofstream out(filename.str().c_str());
      AssertThrow(out,
                  ExcMessage("Could not open continent surface output file: " + filename.str()));

      // Header line documenting the columns.
      out << "# 1:x 2:y ";
      if constexpr (dim == 3)
        out << "3:z 4:radius 5:lon_deg 6:lat_deg 7:area 8:continent_value 9:is_continent 10:speed 11:block_id\n";
      else
        out << "3:area 4:continent_value 5:is_continent 6:speed 7:block_id\n";

      out << std::setprecision(16);

      // One line per top-surface face.
      for (const auto &f : all_faces)
        {
          if constexpr (dim == 3)
            {
              const auto lon_lat = cartesian_to_lon_lat_deg(f.center);

              out << f.center[0] << " "
                  << f.center[1] << " "
                  << f.center[2] << " "
                  << f.center.norm() << " "
                  << lon_lat.first << " "
                  << lon_lat.second << " "
                  << f.area << " "
                  << f.continent_value << " "
                  << (f.is_continent ? 1 : 0) << " "
                  << f.speed << " "
                  << f.block_id << "\n";
            }
          else
            {
              out << f.center[0] << " "
                  << f.center[1] << " "
                  << f.area << " "
                  << f.continent_value << " "
                  << (f.is_continent ? 1 : 0) << " "
                  << f.speed << " "
                  << f.block_id << "\n";
            }
        }
    }



    template <int dim>
    void ContinentStatistics<dim>::write_block_file(const std::vector<BlockSummary<dim>> &blocks,
                                                    const std::vector<Point<dim>> &block_centroids) const
    {
      std::ostringstream filename;
      filename << this->get_output_directory()
               << output_file_prefix
               << "_blocks."
               << std::setw(5) << std::setfill('0')
               << this->get_timestep_number();

      std::ofstream out(filename.str().c_str());
      AssertThrow(out,
                  ExcMessage("Could not open continent block output file: " + filename.str()));

      out << "# 1:block_id ";
      if constexpr (dim == 3)
        out << "2:centroid_x 3:centroid_y 4:centroid_z 5:radius 6:lon_deg 7:lat_deg 8:area 9:mean_speed 10:n_faces\n";
      else
        out << "2:centroid_x 3:centroid_y 4:area 5:mean_speed 6:n_faces\n";

      out << std::setprecision(16);

      // One line per connected continental block.
      for (unsigned int i = 0; i < blocks.size(); ++i)
        {
          const double mean_speed =
            (blocks[i].area > 0.0 ? blocks[i].speed_area_integral / blocks[i].area : 0.0);

          if constexpr (dim == 3)
            {
              const auto lon_lat = cartesian_to_lon_lat_deg(block_centroids[i]);

              out << blocks[i].id << " "
                  << block_centroids[i][0] << " "
                  << block_centroids[i][1] << " "
                  << block_centroids[i][2] << " "
                  << block_centroids[i].norm() << " "
                  << lon_lat.first << " "
                  << lon_lat.second << " "
                  << blocks[i].area << " "
                  << mean_speed << " "
                  << blocks[i].n_faces << "\n";
            }
          else
            {
              out << blocks[i].id << " "
                  << block_centroids[i][0] << " "
                  << block_centroids[i][1] << " "
                  << blocks[i].area << " "
                  << mean_speed << " "
                  << blocks[i].n_faces << "\n";
            }
        }
    }



    template <int dim>
    std::pair<std::string,std::string>
    ContinentStatistics<dim>::execute(TableHandler &statistics)
    {
      // Re-resolve field names every call to be safe.
      resolve_field_indices();

      /**
       * We integrate only on boundary faces, because this is a surface diagnostic.
       *
       * face_quadrature: quadrature rule on faces
       * face_values:     helper object to evaluate FE fields on faces
       */
      const QGauss<dim-1> face_quadrature(this->get_fe().base_element(0).degree + 1);

      FEFaceValues<dim> face_values(this->get_mapping(),
                                    this->get_fe(),
                                    face_quadrature,
                                    update_values |
                                    update_quadrature_points |
                                    update_JxW_values);

      // Velocity values at face quadrature points.
      std::vector<Tensor<1,dim>> velocity_values(face_quadrature.size());

      /**
       * For each selected continent field, store values at the face quadrature points.
       *
       * composition_values_per_field[k][q]
       *   k = selected field index in our local list
       *   q = quadrature point on the current face
       */
      std::vector<std::vector<double>> composition_values_per_field(
        continent_field_indices.size(),
        std::vector<double>(face_quadrature.size()));

      /**
       * Local list of summarized top-surface faces on this MPI rank.
       * We will later gather them onto rank 0.
       */
      std::vector<FaceRecord<dim>> local_faces;

      const auto &velocity_extractor = this->introspection().extractors.velocities;

      // Loop over all active cells.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (!cell->is_locally_owned())
            continue;

          // Loop over the faces of the current cell.
          for (const unsigned int f : cell->face_indices())
            {
              // Keep only top-boundary faces.
              if (!is_top_boundary_face(cell, f))
                continue;

              // Reinitialize FEValues on this face.
              face_values.reinit(cell, f);

              // Evaluate velocity at face quadrature points.
              face_values[velocity_extractor].get_function_values(this->get_solution(),
                                                                  velocity_values);

              // Evaluate each selected continent field at face quadrature points.
              for (unsigned int k = 0; k < continent_field_indices.size(); ++k)
                {
                  const FEValuesExtractors::Scalar comp_extractor =
                    this->introspection().extractors.compositional_fields[continent_field_indices[k]];

                  face_values[comp_extractor].get_function_values(this->get_solution(),
                                                                  composition_values_per_field[k]);
                }

              double face_area = 0.0;
              double face_speed_integral = 0.0;
              double face_continent_integral = 0.0;

              /**
               * Integrate on the current face:
               *   - area
               *   - integral of |u|
               *   - integral of summed continent composition
               */
              for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                {
                  double csum = 0.0;
                  for (unsigned int k = 0; k < continent_field_indices.size(); ++k)
                    csum += composition_values_per_field[k][q];

                  face_area += face_values.JxW(q);
                  face_speed_integral += velocity_values[q].norm() * face_values.JxW(q);
                  face_continent_integral += csum * face_values.JxW(q);
                }

              FaceRecord<dim> rec;
              rec.center = cell->face(f)->center();
              rec.area = face_area;
              rec.continent_value = face_continent_integral / std::max(face_area, 1e-30);
              rec.speed = face_speed_integral / std::max(face_area, 1e-30* year_in_seconds);

              // Binary continent classification using the user threshold.
              rec.is_continent = (rec.continent_value > continent_threshold);

              // Store face vertices for later connectivity/perimeter analysis.
              rec.vertices.reserve(GeometryInfo<dim>::vertices_per_face);
              for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
                rec.vertices.push_back(cell->face(f)->vertex(v));

              local_faces.push_back(rec);
            }
        }

      /**
       * Serialize each local face into a flat vector<double>.
       *
       * This is simple and easy to reconstruct later on rank 0.
       */
      std::vector<double> packed_local;
      packed_local.reserve(local_faces.size() *
                           (4 + dim + dim * GeometryInfo<dim>::vertices_per_face));

      for (const auto &rec : local_faces)
        {
          // 4 scalar values
          packed_local.push_back(rec.area);
          packed_local.push_back(rec.continent_value);
          packed_local.push_back(rec.speed);
          packed_local.push_back(rec.is_continent ? 1.0 : 0.0);

          // face center
          for (unsigned int d = 0; d < dim; ++d)
            packed_local.push_back(rec.center[d]);

          // vertices
          for (const auto &v : rec.vertices)
            for (unsigned int d = 0; d < dim; ++d)
              packed_local.push_back(v[d]);
        }

      /**
       * Gather all packed vectors on rank 0.
       *
       * This is a root-centric algorithm:
       *   - each rank computes its local top faces
       *   - rank 0 reconstructs the global surface map
       *   - rank 0 computes blocks/perimeter and writes files
       *
       * It is simple and good for diagnostics, even if not the most scalable approach.
       */
      const std::vector<std::vector<double>> gathered =
        Utilities::MPI::gather(this->get_mpi_communicator(), packed_local, 0);

      // Final scalar outputs.
      double global_cont_area = 0.0;
      double global_cont_speed_integral = 0.0;
      double global_perimeter = 0.0;
      double global_largest_block = 0.0;
      unsigned int global_n_blocks = 0;
      double fragmentation_index = 0.0;

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          /**
           * Reconstruct all face records on rank 0.
           */
          std::vector<FaceRecord<dim>> all_faces;

          for (const auto &buffer : gathered)
            {
              const unsigned int stride =
                4 + dim + dim * GeometryInfo<dim>::vertices_per_face;

              AssertThrow(buffer.size() % stride == 0,
                          ExcMessage("Continent statistics: invalid gathered buffer size."));

              unsigned int p = 0;
              const unsigned int n_faces = buffer.size() / stride;

              for (unsigned int i = 0; i < n_faces; ++i)
                {
                  FaceRecord<dim> rec;
                  rec.area = buffer[p++];
                  rec.continent_value = buffer[p++];
                  rec.speed = buffer[p++];
                  rec.is_continent = (buffer[p++] > 0.5);

                  for (unsigned int d = 0; d < dim; ++d)
                    rec.center[d] = buffer[p++];

                  rec.vertices.resize(GeometryInfo<dim>::vertices_per_face);
                  for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
                    for (unsigned int d = 0; d < dim; ++d)
                      rec.vertices[v][d] = buffer[p++];

                  all_faces.push_back(rec);
                }
            }

          const unsigned int N = all_faces.size();

          /**
           * 1) Continental area and area-weighted speed integral.
           *
           * Only continent faces contribute.
           */
          for (const auto &f : all_faces)
            {
              if (f.is_continent)
                {
                  global_cont_area += f.area;
                  global_cont_speed_integral += f.speed * f.area;
                }
            }

          /**
           * 2) Build connectivity graph via shared vertices.
           *
           * Two continent faces are merged into the same block
           * if they share at least one vertex.
           *
           * This is slightly permissive, but very practical.
           * Later, you could make it stricter by requiring a shared edge.
           */
          std::map<std::string, std::vector<unsigned int>> vertex_to_faces;
          for (unsigned int i = 0; i < N; ++i)
            for (const auto &v : all_faces[i].vertices)
              vertex_to_faces[point_key(v)].push_back(i);

          DSU dsu(N);

          for (const auto &it : vertex_to_faces)
            {
              const std::vector<unsigned int> &faces = it.second;

              for (unsigned int a = 0; a < faces.size(); ++a)
                for (unsigned int b = a + 1; b < faces.size(); ++b)
                  if (all_faces[faces[a]].is_continent &&
                      all_faces[faces[b]].is_continent)
                    dsu.unite(faces[a], faces[b]);
            }

          /**
           * 3) Compute area of each raw connected component.
           */
          std::map<unsigned int, double> root_area;
          for (unsigned int i = 0; i < N; ++i)
            if (all_faces[i].is_continent)
              root_area[dsu.find(i)] += all_faces[i].area;

          /**
           * 4) Keep only blocks larger than the minimum area threshold.
           * Assign compact block IDs: 0, 1, 2, ...
           */
          std::map<unsigned int, unsigned int> root_to_block_id;
          unsigned int next_block_id = 0;

          for (const auto &it : root_area)
            {
              if (it.second >= minimum_block_area)
                {
                  root_to_block_id[it.first] = next_block_id++;
                  ++global_n_blocks;
                  global_largest_block = std::max(global_largest_block, it.second);
                }
            }

          /**
           * 5) Assign block IDs back to faces.
           *
           * Non-continent faces keep block_id = -1.
           * Continent speckles smaller than minimum_block_area also get -1.
           */
          for (unsigned int i = 0; i < N; ++i)
            {
              if (!all_faces[i].is_continent)
                {
                  all_faces[i].block_id = -1;
                  continue;
                }

              const unsigned int root = dsu.find(i);
              auto it = root_to_block_id.find(root);

              if (it != root_to_block_id.end())
                all_faces[i].block_id = static_cast<int>(it->second);
              else
                all_faces[i].block_id = -1;
            }

          /**
           * 6) Fragmentation index.
           *
           * Here we define:
           *   F = 1 - largest_block_area / total_continent_area
           *
           * Interpretation:
           *   F = 0     -> one dominant coherent continent
           *   F closer to 1 -> more fragmented continent distribution
           */
          if (global_cont_area > 0.0)
            fragmentation_index =
              1.0 - global_largest_block / global_cont_area;

          /**
           * 7) Perimeter estimate (only meaningful in 3D here).
           *
           * We approximate the continent perimeter by summing surface edges
           * that separate at least one continent face from at least one
           * non-continent face.
           */
          if constexpr (dim == 3)
            {
              std::map<std::pair<std::string,std::string>, std::vector<unsigned int>> edge_to_faces;

              // Build edge -> incident faces map.
              for (unsigned int i = 0; i < N; ++i)
                {
                  const unsigned int nv = all_faces[i].vertices.size();

                  for (unsigned int e = 0; e < nv; ++e)
                    {
                      const unsigned int v0 = e;
                      const unsigned int v1 = (e + 1) % nv;

                      std::string k0 = point_key(all_faces[i].vertices[v0]);
                      std::string k1 = point_key(all_faces[i].vertices[v1]);
                      if (k1 < k0)
                        std::swap(k0, k1);

                      edge_to_faces[{k0, k1}].push_back(i);
                    }
                }

              // Sum the lengths of continent/non-continent transition edges.
              for (const auto &it : edge_to_faces)
                {
                  const std::vector<unsigned int> &faces = it.second;

                  bool has_cont = false;
                  bool has_noncont = false;

                  for (const unsigned int idx : faces)
                    {
                      has_cont    = has_cont    || all_faces[idx].is_continent;
                      has_noncont = has_noncont || !all_faces[idx].is_continent;
                    }

                  if (has_cont && has_noncont)
                    {
                      double edge_length = 0.0;
                      bool found = false;

                      // Recover one geometric realization of the edge and compute its length.
                      for (const unsigned int idx : faces)
                        {
                          const auto &verts = all_faces[idx].vertices;

                          for (unsigned int a = 0; a < verts.size(); ++a)
                            {
                              const unsigned int b = (a + 1) % verts.size();

                              std::string ka = point_key(verts[a]);
                              std::string kb = point_key(verts[b]);
                              if (kb < ka)
                                std::swap(ka, kb);

                              if (ka == it.first.first && kb == it.first.second)
                                {
                                  edge_length = point_distance(verts[a], verts[b]);
                                  found = true;
                                  break;
                                }
                            }

                          if (found)
                            break;
                        }

                      global_perimeter += edge_length;
                    }
                }
            }

          /**
           * 8) Build per-block summaries.
           */
          std::vector<BlockSummary<dim>> blocks(global_n_blocks);
          std::vector<Point<dim>> block_centroids(global_n_blocks);

          for (unsigned int b = 0; b < global_n_blocks; ++b)
            blocks[b].id = b;

          for (const auto &f : all_faces)
            {
              if (f.block_id < 0)
                continue;

              BlockSummary<dim> &block = blocks[static_cast<unsigned int>(f.block_id)];
              block.area += f.area;
              block.speed_area_integral += f.speed * f.area;
              block.n_faces += 1;

              for (unsigned int d = 0; d < dim; ++d)
                block.centroid_numerator[d] += f.center[d] * f.area;
            }

          // Compute area-weighted block centroids.
          for (unsigned int b = 0; b < global_n_blocks; ++b)
            {
              if (blocks[b].area > 0.0)
                {
                  for (unsigned int d = 0; d < dim; ++d)
                    block_centroids[b][d] = blocks[b].centroid_numerator[d] / blocks[b].area;
                }
            }

          /**
           * 9) Optional debug outputs.
           */
          if (write_surface_map)
            write_surface_file(all_faces);

          if (write_block_summary)
            write_block_file(blocks, block_centroids);
        }

      /**
       * Broadcast final scalar values back to all ranks,
       * because the TableHandler is used collectively.
       */
      global_cont_area =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_cont_area, 0);
      global_cont_speed_integral =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_cont_speed_integral, 0);
      global_perimeter =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_perimeter, 0);
      global_largest_block =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_largest_block, 0);
      global_n_blocks =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), global_n_blocks, 0);
      fragmentation_index =
        Utilities::MPI::broadcast(this->get_mpi_communicator(), fragmentation_index, 0);

        const double mean_cont_speed_m_per_s =
        (global_cont_area > 0.0 ? global_cont_speed_integral / global_cont_area : 0.0);

        const double mean_cont_speed =
        mean_cont_speed_m_per_s * year_in_seconds;

      /**
       * Add final values to ASPECT's statistics table.
       */
      statistics.add_value("Continental area", global_cont_area);
      statistics.set_precision("Continental area", 8);
      statistics.set_scientific("Continental area", true);

      statistics.add_value("Continental perimeter", global_perimeter);
      statistics.set_precision("Continental perimeter", 8);
      statistics.set_scientific("Continental perimeter", true);

      statistics.add_value("Number of continental blocks", global_n_blocks);
      statistics.set_precision("Number of continental blocks", 0);

      statistics.add_value("Largest continental block area", global_largest_block);
      statistics.set_precision("Largest continental block area", 8);
      statistics.set_scientific("Largest continental block area", true);

      statistics.add_value("Fragmentation index", fragmentation_index);
      statistics.set_precision("Fragmentation index", 6);

      statistics.add_value("Average continental drift speed", mean_cont_speed);
      statistics.set_precision("Average continental drift speed", 8);
      statistics.set_scientific("Average continental drift speed", true);

      /**
       * Compact one-line summary shown in the screen output.
       */
      std::ostringstream out;
      if (output_verbose_screen_line)
        {
          out << "A_cont=" << global_cont_area
              << " m^2, P_cont=" << global_perimeter
              << " m, N_blocks=" << global_n_blocks
              << ", Frag=" << fragmentation_index
              << ", U_cont=" << mean_cont_speed << " m/yr";
        }

      return std::make_pair("Continent statistics", out.str());
    }



    ASPECT_REGISTER_POSTPROCESSOR(
      ContinentStatistics,
      "continent statistics",
      "Computes continent-specific diagnostics on the top surface from a "
      "user-defined set of compositional fields. Outputs continental area, "
      "continent perimeter, number of connected continental blocks, largest "
      "continental block area, fragmentation index, and average continental "
      "drift speed. Can also write per-face and per-block debug files.")
  }
}