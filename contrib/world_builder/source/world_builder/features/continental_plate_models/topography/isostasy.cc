/*
  Copyright (C) 2018-2026 by the authors of the World Builder code.

  This file is part of the World Builder.

  World Builder is free software: you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 2 of the License, or (at your
  option) any later version.
*/

#include "world_builder/features/continental_plate_models/topography/isostasy.h"

#include "world_builder/features/continental_plate_models/topography/interface.h"
#include "world_builder/gravity_model/interface.h"
#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/value_at_points.h"
#include "world_builder/world.h"

#include <cmath>

namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace ContinentalPlateModels
    {
      namespace Topography
      {
        Isostasy::Isostasy(WorldBuilder::World *world_)
          : min_depth(NaN::DSNAN),
            max_depth(NaN::DSNAN),
            maximum_topography_value(NaN::DSNAN),
            operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "isostasy";
        }

        Isostasy::~Isostasy()
          = default;

        void
        Isostasy::declare_entries(Parameters &prm, const std::string & /*parent_name*/)
        {
          prm.declare_entry("", Types::Object(),
                            "Isostatic topography based on vertical reference-density integration.");
          prm.declare_entry("min depth", Types::OneOf(Types::Double(0),
                                                       Types::Array(Types::ValueAtPoints(0.,2)),
                                                       Types::String("")),
                            "Minimum depth over which this topography model applies.");
          prm.declare_entry("max depth", Types::OneOf(Types::Double(std::numeric_limits<double>::max()),
                                                       Types::Array(Types::ValueAtPoints(std::numeric_limits<double>::max(),2)),
                                                       Types::String("")),
                            "Maximum depth over which this topography model applies.");
          prm.declare_entry("maximum topography", Types::Double(20e3),
                            "Guaranteed upper bound for positive topography in meters. ASPECT uses this value to construct a safe spherical manifold.");
        }

        void
        Isostasy::parse_entries(Parameters &prm,
                                const std::vector<Point<2>> &coordinates)
        {
          min_depth_surface = Objects::Surface(prm.get("min depth",coordinates));
          min_depth = min_depth_surface.minimum;
          max_depth_surface = Objects::Surface(prm.get("max depth",coordinates));
          max_depth = max_depth_surface.maximum;
          maximum_topography_value = prm.get<double>("maximum topography");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));

          WBAssertThrow(maximum_topography_value >= 0.0,
                        "The maximum isostatic topography must be nonnegative.");
        }

        double
        Isostasy::get_topography(const Point<3> &position_in_cartesian_coordinates,
                                 const Objects::NaturalCoordinate & /*position_in_natural_coordinates*/,
                                 const double topography) const
        {
          WBAssertThrow(std::isfinite(world->compensation_pressure),
                        "The isostasy topography model requires a valid top-level 'reference profile point'.");
          WBAssertThrow(world->number_integration_points > 1,
                        "The number of isostasy integration points must be greater than one.");

          const double gravity =
            world->parameters.gravity_model->gravity_norm(position_in_cartesian_coordinates);
          WBAssertThrow(std::fabs(gravity) > 0.0,
                        "Gravity must be nonzero for isostatic topography.");

          const double dz = world->compensation_depth /
                            static_cast<double>(world->number_integration_points - 1);
          double local_pressure = 0.0;
          double previous_density =
            world->density(position_in_cartesian_coordinates.get_array(), 0.0);

          for (unsigned int i = 1; i < world->number_integration_points; ++i)
            {
              const double current_density =
                world->density(position_in_cartesian_coordinates.get_array(), i * dz);
              local_pressure +=
                0.5 * (previous_density + current_density) * dz * gravity;
              previous_density = current_density;
            }

          const double isostatic_topography =
            (world->compensation_pressure - local_pressure) /
            (world->background_density * gravity);

          WBAssertThrow(isostatic_topography <= maximum_topography_value,
                        "Computed isostatic topography (" << isostatic_topography
                        << " m) exceeds the configured maximum topography ("
                        << maximum_topography_value << " m).");

          return apply_operation(operation, topography, isostatic_topography);
        }

        double
        Isostasy::maximum_topography() const
        {
          return maximum_topography_value;
        }

        WB_REGISTER_FEATURE_CONTINENTAL_PLATE_TOPOGRAPHY_MODEL(Isostasy, isostasy)
      }
    }
  }
}
