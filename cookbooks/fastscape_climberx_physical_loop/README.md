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

## Continued sequence and runtime comparison

`run_coupling_sequence.py` continues every climate and ice component from its
previous restart, continues ASPECT and Fastscape from their checkpoint, and
returns only the new surface increment at each boundary. It also runs the same
ASPECT temperature and mantle-flow problem without surface coupling for a
measured timing comparison:

```bash
python3 run_coupling_sequence.py --windows 3
```

The output includes `timing-summary.json` and `runtime-comparison.png`.

On the local three-window test, 60 ASPECT/Fastscape years and four consecutive
climate years completed without a model instability. Measured wall times were
501.1 seconds for the full sequence, 16.1 seconds for its three coupled ASPECT
windows, 10.1 seconds for ASPECT without surface coupling but with the same
three restart boundaries, and 4.0 seconds for one continuous uncoupled ASPECT
run. Thus Fastscape and mesh deformation made the matched, windowed ASPECT
part 1.6 times slower. The complete loop was 124.5 times slower than continuous
ASPECT because CLIMBER-X and Yelmo used 483.4 seconds while this deliberately
coarse ASPECT problem used only seconds. That full ratio is not transferable
to a refined production model; it depends mainly on ASPECT resolution and how
often the climate model is called.
