# CLIMBER-X–FastScape native surface exchange

This directory defines a small binary coupling boundary for climate,
landscape, and geodynamic fields. It avoids opening and rewriting large
NetCDF restart files during every coupling window. NetCDF remains the archival
and recovery format.

CLIMBER-X writes `precipitation_rate`, `surface_temperature`, and
`surface_elevation` when `CLIMBERX_CLIMATE_EXCHANGE_FILE` names an output file.
The adapter converts those live fields into ASPECT structured inputs:

```bash
python3 surface_exchange.py climate-to-aspect climate.cxe \
  --erosion-strength erosion-strength.txt \
  --surface-runoff surface-runoff.txt
```

After ASPECT/FastScape writes its surface table, create the compact return
field on the climate grid:

```bash
python3 surface_exchange.py surface-to-climate \
  --surface surface-00002.csv --climate climate.cxe \
  --output topography.cxe
```

Set `CLIMBERX_TOPOGRAPHY_EXCHANGE_FILE=topography.cxe` for the next CLIMBER-X
window. CLIMBER-X interpolates the increment onto its high-resolution
geography before geography initialization. This preserves the existing
reference geography while allowing coastlines and orography to respond.

## Verified local cycle

The complete path was exercised with the one-year CLIMBER-X climate case and
the 1,536-cell global FastScape surface:

1. CLIMBER-X wrote a 2,592-point, 101 KiB exchange directly from its live
   atmosphere fields.
2. The adapter produced ASPECT's built-in spherical structured-data inputs.
3. ASPECT ran temperature and Stokes solves while FastScape applied spatial
   runoff and erosion, surface-state advection, hillslope diffusion, and marine
   sediment transport for two 1,000-year steps.
4. The adapter remapped FastScape's incremental topography to the 5-degree
   climate grid.
5. A second one-year CLIMBER-X run read this return field and interpolated it
   onto the 2,880 by 1,440 geography before initialization.

The active ASPECT run delivered 417 and 409 million cubic metres of sediment
per year to the coast in the two steps. The deliberately short feedback test
returned a maximum 0.068 m topographic change and produced a maximum 0.012 K
temperature difference. These values validate the plumbing; they are not a
scientific calibration.

The reader accepts either byte order. This matters on systems where
`GFORTRAN_CONVERT_UNIT` overrides the byte order requested on a Fortran stream;
the return file preserves the byte order detected in the climate exchange.
