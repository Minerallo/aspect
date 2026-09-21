/*
  Copyright (C) 2018-2026 by the authors of the World Builder code.

  This file is part of the World Builder.

  World Builder is free software: you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 2 of the License, or (at your
  option) any later version.
*/

#ifndef WORLD_BUILDER_TYPES_COMPOSITION_PROPERTY_H
#define WORLD_BUILDER_TYPES_COMPOSITION_PROPERTY_H

#include "world_builder/types/interface.h"

#include <string>
#include <vector>

namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {
    /** Schema type for a composition index and its reference density. */
    class CompositionProperty final: public Interface
    {
      public:
        CompositionProperty();
        CompositionProperty(const CompositionProperty &) = default;
        ~CompositionProperty() override final = default;

        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;

        static constexpr double get_default_reference_density()
        {
          return 3300.0;
        }

      protected:
        CompositionProperty *clone_impl() const override final
        {
          return new CompositionProperty(*this);
        }

      private:
        std::vector<std::string> required;
    };
  }
}

#endif
