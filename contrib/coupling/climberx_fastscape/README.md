# CLIMBER-X–FastScape native surface exchange

This directory defines a small binary coupling boundary for climate,
landscape, and geodynamic fields. It avoids opening and rewriting large
NetCDF restart files during every coupling window. NetCDF remains the archival
and recovery format.

CLIMBER-X writes `precipitation_rate`, `surface_temperature`,
`surface_elevation`, `ice_thickness`, `basal_ice_velocity`, and
`grounded_ice_fraction` when `CLIMBERX_CLIMATE_EXCHANGE_FILE` names an output
file. Dynamic-ice runs obtain basal velocity directly from either Yelmo or
SICOPOLIS. Prescribed-ice runs export ice thickness but use zero basal
velocity, so they do not silently invent ice motion. The adapter converts the
live fields into ASPECT structured inputs:

```bash
python3 surface_exchange.py climate-to-aspect climate.cxe \
  --surface-topography surface-topography.txt \
  --erosion-strength erosion-strength.txt \
  --surface-runoff surface-runoff.txt \
  --ice-thickness ice-thickness.txt \
  --basal-ice-velocity basal-ice-velocity.txt
```

FastScape uses the standard sliding law `E = K u^m` only where ice thickness
exceeds the configured threshold. Both river and glacial erosion enter the
same sediment-routing and deposition calculation, while separate result
fields preserve their individual contributions.

Regional box models can additionally enable `Enable regional ice load
response`. Imported ice then drives local subsidence and unloading drives
rebound. This velocity is added to ASPECT's tectonic surface velocity before
Fastscape advances, so the landscape and sediment-routing calculations use
the displaced bedrock. The delayed displacement is preserved in checkpoints.
See the
[`fastscape_regional_ice_load`](../../../cookbooks/fastscape_regional_ice_load/README.md)
cookbook for a loading, unloading, and disabled-control comparison. This
reduced regional response cannot be enabled with the global degree-two
self-gravity correction.

For prescribed-ice climate runs, the exported basal velocity is zero. A
deliberate sensitivity test can assign a constant velocity only in grounded
ice cells with `--prescribed-basal-ice-velocity VALUE`. This is an explicit
model assumption that needs calibration; it is not used automatically.

After ASPECT/FastScape writes its surface table, create the compact return
field on the climate grid:

```bash
python3 surface_exchange.py surface-to-climate \
  --surface surface-00002.csv --climate climate.cxe \
  --polar-wander-history true_polar_wander.csv \
  --output topography.cxe
```

Set `CLIMBERX_TOPOGRAPHY_EXCHANGE_FILE=topography.cxe` for the next CLIMBER-X
window. CLIMBER-X interpolates the increment onto its high-resolution
geography after loading any geography restart. This preserves both the
continued state and the new surface increment while allowing coastlines and
orography to respond.

When a polar-wander history is supplied, the adapter interprets CLIMBER-X's
longitude and latitude around the latest spin axis and samples the body-fixed
Fastscape surface at those locations. The next climate window therefore sees
the returned Fastscape topography increment in the new spin frame. Omitting
the option uses the geographic north pole and exactly recovers the original
remapping. CLIMBER-X's high-resolution reference geography remains fixed;
rotating that complete reference geography is a separate step required for a
fully self-consistent large polar displacement.
Four nearby Fastscape cells are combined with spherical-distance weights so
small pole movements produce continuous changes instead of nearest-cell
jumps. Set `--interpolation-neighbors 1` to reproduce the original mapping.

The reader accepts either byte order. This matters on systems where
`GFORTRAN_CONVERT_UNIT` overrides the byte order requested on a Fortran stream;
the return file preserves the byte order detected in the climate exchange.
