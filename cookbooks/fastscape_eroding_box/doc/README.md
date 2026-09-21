```{tags}
category:cookbook
feature:3d
feature:cartesian
feature:mesh-deformation
```

(sec:cookbooks:fastscape_eroding_box)=
# Fastscape Eroding Box

*This section was contributed by Derek Neuharth, Esther Heckenbach, and Anne Glerum.*

This is a simple model of an eroding central block. The model is primarily used to test that the coupling installation is working properly, or for simple tests when making changes to the plugin to ensure everything is still working properly.

## FastScape C++ comparison

`fastscape_cpp_eroding_box.prm` reproduces the main physical parameters of
the original Fortran-coupled model with the C++ FastScape library: a 2.5 km
central block, a river-incision coefficient of $10^{-4}$, drainage-area
exponent 0.4, slope exponent 1, hillslope diffusivity $10^{-2}$ square meters
per year, five 10,000-year landscape steps per 50,000-year ASPECT step, and a
total duration of 250,000 years.

The C++ model uses an independent 160 by 160 landscape grid with 25,600 cells
and approximately 0.7 km spacing. The ASPECT volume mesh has only 160 cells.
An eight-cell conservative stencil restricts the fine surface motion to the
coarse ASPECT surface. A small deterministic multiscale perturbation replaces
the random seed used by the Fortran coupling, making the drainage network
reproducible.

This is a close physical counterpart rather than a bit-for-bit reproduction.
The current C++ coupling uses single-direction flow routing, whereas the
Fortran cookbook requests variable multi-direction flow, and it does not have
an equivalent of the Fortran bedrock deposition coefficient. Those differences
mainly affect routing across the initially flat lowland and the treatment of
eroded material; the raised block nevertheless develops the expected branching
valleys.

`fastscape_cpp_eroding_box_glacial.prm` includes the same file and adds
glacial erosion over the interior of the raised block. Ice slides at 25 meters
per year, the glacial erosion coefficient is $2\times10^{-5}$, and the minimum
active ice thickness is 500 meters. The two cases therefore differ only by
the glacial erosion coefficient.

Run both cases from the cookbook directory:

```sh
OMP_NUM_THREADS=1 /path/to/aspect fastscape_cpp_eroding_box.prm
OMP_NUM_THREADS=1 /path/to/aspect fastscape_cpp_eroding_box_glacial.prm
```

The reference runs completed all five coupled steps. The normal case develops
a branching river network and retains a maximum elevation of 2486 m. The
glacial case reaches 2361 m and produces 5 m of glacial erosion in the final
10,000-year landscape step. Its peak resident memory was 274 MB with one
process and one thread. Glacial erosion reorganizes several drainage divides,
so local elevation differences include both direct ice erosion and subsequent
river capture. This is a physics and coupling test, not a calibrated glacier.

Create the comparison figure with:

```sh
python3 compare_fastscape_cpp_eroding_box.py extract \
  output-fastscape-cpp-eroding-box/fastscape_surface_evolution/fastscape-00005.vtu \
  output-fastscape-cpp-eroding-box-glacial/fastscape_surface_evolution/fastscape-00005.vtu \
  fastscape-cpp-eroding-box-comparison.npz
python3 compare_fastscape_cpp_eroding_box.py plot \
  fastscape-cpp-eroding-box-comparison.npz \
  fastscape-cpp-eroding-box-comparison.png
```

The river lines in the upper panels mark drainage areas of $10^8$,
$10^{8.5}$, and $10^9$ square meters. The cyan outline marks the 500 m ice
threshold. The generated figure is intentionally not stored in the source
tree.

## Routed-glacier approximation

`fastscape_cpp_eroding_box_routed_glacier.prm` provides an intermediate model
between static prescribed ice fields and a continuum ice-dynamics solver. It
computes elevation-dependent accumulation above an equilibrium-line altitude
(ELA), ablation below it, and routes the remaining ice volume over the
evolving FastScape flow graph. Empirical power laws convert ice discharge
$Q_i$ to glacier width $W$ and thickness $H$. Basal speed follows from

$$u_b = \frac{Q_i}{W H},$$

and the existing glacial erosion law uses this speed. An optional compact
Gaussian graph-distance kernel distributes the centerline erosive volume over
the empirical glacier width. The distribution conserves eroded volume; it
approximates valley widening without solving ice stresses and should therefore
be treated as a geomorphic parameterization rather than an ice-sheet model.

Select the mode with `Glacial erosion mode = routed glacier`. All routed
glacier parameters are optional and the default `prescribed fields` mode is
unchanged. `Terminate routed glacier at sea level = true` treats ice reaching
the coast as exported or calved, preventing a large discharge from following
a flat numerical ocean boundary. Disable it for grounded marine-ice tests.

The spatial CSV and VTU output adds `routed_ice_discharge`,
`routed_glacier_width`, and `routed_ice_mass_balance`. The existing
`ice_thickness`, `basal_ice_velocity`, `glacial_erosion`, and sediment fields
contain the corresponding modeled values and landscape response.

Run the example from this cookbook directory with:

```sh
OMP_NUM_THREADS=1 /path/to/aspect \
  fastscape_cpp_eroding_box_routed_glacier.prm
```

## Lithology, provenance, and depositional layers

`fastscape_cpp_lithology_provenance.prm` assigns 70 percent of the landscape
cells to granite and 30 percent to limestone using a reproducible random seed.
The limestone river-incision factor is two, while granite is the reference
rock with factor one. This represents distinct rocks within the same upper
crust; it does not require separate ASPECT compositional fields.

Eroded sediment carries its source-rock class through river routing, coastal
delivery, marine transport, deposition, and later re-erosion. The surface CSV
files and visualization files contain the flux and deposited fraction of every
class. Each `stratigraphy-*.csv` file is one dated depositional increment, so
the sequence can be stacked to inspect source changes through a sedimentary
basin. The files record net sediment thickness preserved since the preceding
output, so material deposited and re-eroded within that interval is not kept
as a layer. Compaction, chemical transformation, and erosion surfaces are not
yet represented as separate layer objects.

The probabilities are used only when the initial bedrock map is made. A rock
class is therefore stable between time steps. For a mapped geological model,
the next extension should replace this initializer with a structured bedrock
class file while keeping the same transport and output representation.

Run and plot the demonstration with:

```sh
OMP_NUM_THREADS=1 /path/to/aspect fastscape_cpp_lithology_provenance.prm
python3 plot_lithology_provenance.py \
  output-fastscape-cpp-lithology-provenance/fastscape_surface_evolution/fastscape-00004.vtu \
  fastscape-cpp-lithology-provenance.png
```

The generated source-rock and deposited-sediment comparison is intentionally
not stored in the source tree.
