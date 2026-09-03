/*
  Copyright (C) 2018-2026 by the authors of the World Builder code.

  This file is part of the World Builder.

  World Builder is free software: you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 2 of the License, or (at your
  option) any later version.
*/

#ifndef WORLD_BUILDER_FEATURES_CONTINENTAL_PLATE_MODELS_TOPOGRAPHY_ISOSTASY_H
#define WORLD_BUILDER_FEATURES_CONTINENTAL_PLATE_MODELS_TOPOGRAPHY_ISOSTASY_H

#include "world_builder/features/continental_plate_models/topography/interface.h"
#include "world_builder/features/feature_utilities.h"
#include "world_builder/objects/surface.h"

namespace WorldBuilder
{
  namespace Features
  {
    using namespace FeatureUtilities;

    namespace ContinentalPlateModels
    {
      namespace Topography
      {
        /** Isostatic topography from a density-column pressure difference. */
        class Isostasy final: public Interface
        {
          public:
            Isostasy(WorldBuilder::World *world);
            ~Isostasy() override final;

            static void declare_entries(Parameters &prm,
                                        const std::string &parent_name = "");
            void parse_entries(Parameters &prm,
                               const std::vector<Point<2>> &coordinates) override final;

            double get_topography(const Point<3> &position,
                                  const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                  double topography) const override final;

            double maximum_topography() const override final;

          private:
            double min_depth;
            Objects::Surface min_depth_surface;
            double max_depth;
            Objects::Surface max_depth_surface;
            double maximum_topography_value;
            Operations operation;
        };
      }
    }
  }
}

#endif
