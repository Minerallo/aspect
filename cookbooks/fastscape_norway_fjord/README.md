# Norway-like fjord incision and unloading

This inexpensive regional process test combines Fastscape C++ glacial erosion,
river and marine sediment routing, and the ASPECT regional ice-load response.
It represents a synthetic western Norwegian mountain-to-fjord transect rather
than a reconstruction of a named fjord.

The 160 km by 100 km domain has a 1,200 m pre-existing overdeepened trough and
an approximately 1.7 km landscape grid over a
coarser ASPECT volume mesh. A branching valley glacier reaches 1,400 m thickness
and 60 m/yr basal speed. With the linear erosion law and a coefficient of
`3e-5`, the fastest sliding corresponds to 1.8 mm/yr erosion. The speed is close
to the 10--15 cm/day winter basal sliding measured beneath Engabreen; the ice
geometry itself is an idealized large-glaciation forcing. The ice-load
response uses ice and compensation densities of 917 and 3,300 kg/m3 and a
5,000-year regional relaxation time. These are transparent process-test values;
they have not been calibrated against one Norwegian fjord or a full mantle
viscosity model.

Landscape evolution is limited to 100-year substeps. This prevents a long
geodynamic step from placing a large pulse in one coastal cell before the
marine diffusion operator can distribute it offshore.

Run from this directory:

```sh
python3 generate_ice_fields.py
../../../builts/fastscape-release/aspect fastscape_norway_glaciation.prm
../../../builts/fastscape-release/aspect fastscape_norway_unloading.prm
../../../builts/fastscape-release/aspect fastscape_norway_no_glacial_erosion.prm
python3 check_norway_fjord.py
python3 plot_norway_fjord_pyvista.py
python3 render_cinematic_fjord.py
```

The loading run evolves the glacier-covered landscape for 100,000 years and
writes a checkpoint. The unloading run restarts it, removes both ice thickness
and basal sliding, and follows 10,000 years of rebound. The control retains the
same ice load but disables glacial erosion, isolating the incision signal.

The visualization script constructs and triangulates all surfaces with PyVista,
writes `norway_fjord_results.vtm` for interactive inspection, and produces the
four-panel `norway_fjord_pyvista.png`. In a headless terminal it renders the
PyVista triangles with Matplotlib; pass `--native-pyvista` on a machine with a
graphical or OSMesa VTK backend to use PyVista's native renderer.

`render_cinematic_fjord.py` reads all 22 surface time steps, constructs each
terrain with PyVista, and produces a fixed-camera MP4 plus a poster image. It
shows a translucent ice surface at bed elevation plus ice thickness, basal
motion arrows, active glacial erosion, sea level, sediment flux, and unloading
rebound. The default rendering uses every second landscape point and streams
frames to ffmpeg to keep memory use bounded.

Scientific context:

- Steer et al. (2012), *Bimodal Plio-Quaternary glacial erosion of fjords and
  low-relief surfaces in Scandinavia*, https://doi.org/10.1038/ngeo1549.
- Gong et al. (2017), *Basal conditions at Engabreen, Norway, inferred from
  surface measurements and inverse modelling*,
  https://doi.org/10.1017/jog.2017.78.
