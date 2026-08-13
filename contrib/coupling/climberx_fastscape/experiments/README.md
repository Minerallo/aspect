# Local uplift and glacial erosion comparison

This paired experiment uses a 100 km by 100 km surface above a 50 km deep
ASPECT box. A 30 km by 30 km block at its center rises at 5 mm per year. Both runs
use the same FastScape river erosion and prescribed ice fields; only the
glacial erosion coefficient differs.

The landscape grid contains 64 by 64 cells (1.5625 km spacing), while the
ASPECT volume mesh contains only eight active cells. ASPECT advances in
1,000-year steps and FastScape uses two 500-year steps per ASPECT step. The
eight-point conservative transfer preserves separate area integrals of
positive and negative surface motion when restricting the fine landscape
result to the coarse ASPECT surface. The experiment ends after 10,000 years.
The bottom is an open boundary so that the model can balance the imposed
central uplift without an incompatible prescribed velocity.

Run the pair from this directory with a built ASPECT executable:

```sh
OMP_NUM_THREADS=1 /path/to/aspect local_uplift_no_glacial.prm
OMP_NUM_THREADS=1 /path/to/aspect local_uplift_with_glacial.prm
```

Extract the final surface fields using a Python installation with VTK, then
plot them using a Python installation with Matplotlib:

```sh
python3 compare_local_uplift.py extract \
  output-local-uplift-no-glacial-high-resolution/fastscape_surface_evolution/fastscape-00010.vtu \
  output-local-uplift-with-glacial-high-resolution/fastscape_surface_evolution/fastscape-00010.vtu \
  local-uplift-comparison.npz
python3 compare_local_uplift.py plot \
  local-uplift-comparison.npz local-uplift-comparison.png
```

In the reference run, glacial erosion is active in cells where ice is at least
500 m thick. After 10,000 years those cells are 4.9 to 10 m lower than in the
run without glacial erosion. Peak sediment flux is 3.59 times larger and the
mean flux over the glacier-active cells is 6.98 times larger. The final
500-year landscape step removes at most 0.5 m of rock in one cell, and the
cumulative glacial erosion is 46.35 cubic kilometers. The one-process reference
run used about 187 MB of peak resident memory. These values check that the
coupling responds consistently; they are not a calibration to a particular
glacier or landscape.

The cyan contour in `local-uplift-comparison.png` marks the 500 m ice-thickness
threshold.
