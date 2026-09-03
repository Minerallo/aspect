/*
  Copyright (C) 2018-2026 by the authors of the World Builder code.

  This file is part of the World Builder.

  World Builder is free software: you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 2 of the License, or (at your
  option) any later version.
*/

#include "world_builder/types/composition_property.h"
#include "world_builder/parameters.h"

namespace WorldBuilder
{
  namespace Types
  {
    CompositionProperty::CompositionProperty()
      : required({"index"})
    {
      this->type_name = Types::type::Object;
    }

    void
    CompositionProperty::write_schema(Parameters &prm,
                                      const std::string &name,
                                      const std::string &documentation) const
    {
      using namespace rapidjson;
      Document &declarations = prm.declarations;
      const std::string base = prm.get_full_json_path() + "/" + name;

      Pointer((base + "/type").c_str()).Set(declarations,"object");
      Pointer((base + "/description").c_str()).Set(declarations,documentation.c_str());
      Pointer((base + "/additionalProperties").c_str()).Set(declarations,false);

      for (unsigned int i = 0; i < required.size(); ++i)
        Pointer((base + "/required/" + std::to_string(i)).c_str()).Set(declarations,
                                                                        required[i].c_str());

      Pointer((base + "/properties/index/type").c_str()).Set(declarations,"integer");
      Pointer((base + "/properties/index/minimum").c_str()).Set(declarations,0);
      Pointer((base + "/properties/index/description").c_str()).Set(declarations,
                                                                      "Composition index used in composition lookups.");

      Pointer((base + "/properties/name/type").c_str()).Set(declarations,"string");
      Pointer((base + "/properties/name/description").c_str()).Set(declarations,
                                                                     "Optional human-readable composition name.");

      Pointer((base + "/properties/reference density/type").c_str()).Set(declarations,"number");
      Pointer((base + "/properties/reference density/default value").c_str()).Set(
        declarations, CompositionProperty::get_default_reference_density());
      Pointer((base + "/properties/reference density/description").c_str()).Set(
        declarations, "Reference density in kg/m^3 used for isostatic column integration.");
    }
  }
}
