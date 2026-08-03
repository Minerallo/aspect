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
threshold.

```{figure} ../fastscape-cpp-eroding-box-comparison.png
:width: 100%

Fine-grid elevation, drainage area, and glacial erosion after 250,000 years.
```
