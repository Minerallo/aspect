# Regional ice loading with Fastscape C++

This 100 km by 100 km example combines landscape evolution with a reduced
regional solid-Earth response. A 1,000 m ice block loads the central area for
5,000 years. The run is then restarted with the ice removed and rebounds for
another 5,000 years. A separate control keeps the response disabled.

The local equilibrium displacement is

```text
displacement = -ice density * ice thickness / compensation density.
```

An optional immediate fraction is applied first. The remainder approaches
equilibrium exponentially over the chosen relaxation time. The resulting
vertical velocity is added to ASPECT's tectonic surface velocity before
Fastscape advances. Consequently, river incision, glacial erosion, sediment
transport, subsidence, and rebound all use the same evolving elevation.

Run the loading, unloading, and control stages from this directory:

```bash
aspect fastscape_regional_ice_load_loading.prm
aspect fastscape_regional_ice_load_unloading.prm
aspect fastscape_regional_ice_load_no_response.prm
python3 check_regional_ice_load.py
python3 plot_regional_ice_load.py
```

The delayed displacement is stored in ASPECT checkpoints, so a later climate
window can replace the ice-thickness file without losing the loading history.
The option is disabled by default and cannot be combined with the global
degree-two self-gravity correction, which prevents counting the solid-Earth
load response twice.

This is a local isostatic approximation, not a flexural plate or complete
self-gravitating sea-level model. It currently responds to imported ice. A
future sediment-load response should use Fastscape's explicitly conserved
erosion and deposition budget rather than infer loading from elevation, which
would incorrectly include tectonic uplift.
