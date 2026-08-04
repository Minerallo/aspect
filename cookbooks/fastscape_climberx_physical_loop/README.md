# Closed physical climate–surface–solid-Earth loop

This experiment runs one complete feedback cycle using the local installations:

1. CLIMBER-X evolves atmosphere, ocean, sea ice, land, surface and basal ice
   mass balance, and the Yelmo ice-sheet model for one year.
2. A compact binary exchange passes precipitation, surface temperature,
   topography, ice thickness, grounded-ice fraction, and basal sliding speed
   directly to ASPECT and Fastscape. Large climate restart files are not used
   as the coupling interface.
3. ASPECT solves mantle temperature and velocity. Fastscape applies river and
   glacial erosion, hillslope transport, sediment routing, marine sediment
   transport and deposition, and tangential advection of its surface state
   while deforming ASPECT's outer surface.
4. ASPECT's changing density and surface loads update the principal inertia
   axis. An Earth-scale rotational bulge limits true polar wander, and the
   degree-two load response represents the surface deformation that is absent
   from the short, coarse mantle calculation.
5. The surface increment is interpolated back into the current spin-axis
   coordinate system and passed to a second CLIMBER-X year. Its atmosphere,
   ocean, sea ice, land, and dynamic ice then respond to the changed geography.

Run from any directory:

```bash
python3 /Users/ponsm/Desktop/software/aspect_install/aspect_fatscapecc/aspect/cookbooks/fastscape_climberx_physical_loop/run_physical_loop.py
```

The driver runs only one process at a time, fixes the climate model to one
thread, and stops a process if system-wide free memory falls below 15 percent.
Results are written to `output-physical-loop`, which Git ignores. The main
visual result is `physical-loop-summary.png`; `summary.json` records the
surface, sediment, exchange, and polar-wander files.

This short run proves that the feedback path is live. Its one-year climate
windows, 20-year surface window, coarse grids, erosion coefficients,
degree-two load parameters, and component restart strategy require calibration
before the magnitudes can be interpreted as an Earth reconstruction. A longer
experiment should continue component restarts between windows instead of
starting the feedback and control climate years from the same initial state.
