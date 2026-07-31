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
